// ===========================================================================
//  Mp_p-1_gpu  --  GPU P-1 factoring of Mersenne numbers  M_p = 2^p - 1
// ===========================================================================
//
//  P-1 stage 1: with E the product of the largest power of each prime <= B1,
//  compute x = 3^E mod M_p and take gcd(x-1, M_p). A prime factor q = 2kp+1 of
//  M_p is found exactly when k is B1-smooth.
//
//  The bignum GCD is a self-contained recursive half-GCD -- no GMP dependency.
//
//  LICENCE / PROVENANCE. The modular arithmetic -- the whole OpenCL FFT/NTT
//  engine, its kernels, and the tuner -- is derived from gpuowl / PRPLL by
//  Mihai Preda and George Woltman, ported to MSVC. That project is GPLv3, so
//  this one is too, and the per-file "Copyright (C) Mihai Preda" notices must be
//  kept: GPLv3 section 5 requires preserving them, and removing them would make
//  the result undistributable. See LICENCE.md.
//
//  Build:  build.bat        (needs only Visual Studio; no CUDA/OpenCL SDK)
//  Run:    Mp_p-1_gpu.exe [--config config.txt] [--selftest[=which]]
// ===========================================================================

#include "Args.h"
#include "Background.h"
#include "BigInt.h"
#include "Bounds.h"
#include "Config.h"
#include "Context.h"
#include "FFTConfig.h"
#include "Gcd.h"
#include "Gpu.h"
#include "GpuCommon.h"
#include "PM1.h"
#include "Queue.h"
#include "Stage2Plan.h"
#include "TrigBufCache.h"
#include "clwrap.h"
#include "common.h"
#include "log.h"
#include "timeutil.h"
#include "tune.h"

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>

#include <io.h>

#include <atomic>
#include <filesystem>
#include <cstdio>
#include <ctime>
#include <exception>
#include <string>
#include <vector>

int runBigIntTests();
int runExponentTest();
int runEngineTests(GpuCommon shared, Queue* q, const std::string& fftSpec, bool quick);
int runPM1Tests(GpuCommon shared, Queue* q, const std::string& fftSpec);
int runExtendTests(GpuCommon shared, Queue* q, const std::string& fftSpec);
int runPP1Tests(GpuCommon shared, Queue* q, const std::string& fftSpec);
int runStage2Tests(GpuCommon shared, Queue* q, const std::string& fftSpec);
FFTConfig chooseVerifiedFFT(GpuCommon shared, Queue* q, u32 E,
                            const std::string& forcedSpec, bool verify);
void recommendFFT(GpuCommon shared, Queue* q, u32 E);

std::atomic<bool> gInterrupted{false};

namespace {

BOOL WINAPI ctrlHandler(DWORD type) {
  if (type == CTRL_C_EVENT || type == CTRL_BREAK_EVENT || type == CTRL_CLOSE_EVENT) {
    gInterrupted.store(true);
    printf("\n  interrupt received, finishing current step...\n");
    return TRUE;
  }
  return FALSE;
}

// The stage-2 pairing shape, sized to the GPU's memory. Shared by the real run
// and by --bounds: the shape sets muls-per-prime, which is an input to the cost
// model, so if the two disagreed --bounds would recommend bounds the run would
// not actually choose.
//
// Queries the device but does not create a Gpu, so it is cheap and safe to call
// from a diagnostic.
// The smallest transform that can hold E, with no GPU involved. Only for
// diagnostics: the real run takes the FFT from chooseVerifiedFFT.
FFTConfig smallestFittingFFT(u32 E) {
  for (const FFTShape& shape : FFTShape::allShapes()) {
    FFTConfig f{shape, 101, CARRY_AUTO};
    if (f.maxExp() >= E && double(E) / double(f.size()) >= f.minBpw()) { return f; }
  }
  throw std::runtime_error("no FFT fits this exponent");
}

// `residueBytes` is the size of one T_j buffer: fft.size() * sizeof(Word), taken
// from the transform that will actually be used. It used to be estimated as
// exponent/18 words of 4 bytes, which was wrong twice over -- Word is 64-bit,
// and the NTT transforms at these sizes carry ~39 bits/word, not 18. The two
// errors happened to cancel to within 8%, which is not a property to rely on.
Stage2Shape pickStage2Shape(const Config& cfg, int deviceOverride, u64 residueBytes) {
  if (cfg.stage2D) { return {cfg.stage2D, cfg.stage2W ? cfg.stage2W : 1}; }
  if (cfg.stage2Mode == STAGE2_OFF) { return {210, 1}; }

  u64 budget = 0;
  auto devices = getAllDeviceIDs();
  const int devIx = (deviceOverride >= 0) ? deviceOverride : (cfg.device >= 0 ? cfg.device : 0);
  if (devIx >= 0 && size_t(devIx) < devices.size()) {
    const cl_device_id dev = devices[devIx];
    if (hasFreeMemInfo(dev)) {
      // CL_DEVICE_GLOBAL_FREE_MEMORY_AMD. Leave room for the transform's own
      // working buffers -- running out of VRAM mid-walk throws the walk away.
      const u64 freeMem = getFreeMem(dev);
      budget = freeMem > (1ull << 30) ? freeMem - (1ull << 30) : 0;
    } else {
      // That extension is AMD-only, so on nVidia it always fails. Falling back
      // to a fixed small budget put an 8 GB card on the floor shape (D=210 w=1,
      // 24 buffers) and cost ~25% of stage-2 throughput. Total memory is
      // queryable everywhere, so size from that instead: keep 2 GB for the
      // transform, the display and the driver, and never claim more than 60%.
      const u64 total = u64(getGpuRamGB(dev) * double(1ull << 30));
      if (total > (2ull << 30)) {
        budget = std::min(total - (2ull << 30), u64(total * 0.6));
      }
    }
  }
  if (!budget) { budget = 512ull << 20; }   // no device info at all: stay small

  Stage2Shape s = chooseStage2Shape(budget, residueBytes, cfg.b1 ? cfg.b1 : 100000);
  if (!s.d) { s = {210, 1}; }
  if (!s.w) { s.w = 1; }
  return s;
}

void banner() {
  printf("=======================================================================\n");
  printf(" %s %s  --  GPU P-1 factoring of M_p = 2^p - 1\n",
         PROGRAM_NAME, PROGRAM_VERSION);
  printf("=======================================================================\n\n");
}

void usage() {
  printf(
"usage: Mp_p-1_gpu.exe [mode] [options]\n"
"\n"
"MODES  (default: run the job described in config.txt)\n"
"  --config <file>        read the job from <file> instead of config.txt\n"
"  --list-devices         list the OpenCL GPUs and exit\n"
"  --selftest[=which]     run self-checks; omit =which to run them all\n"
"      =gcd                 BigInt and GCD          (no GPU used)\n"
"      =exponent            stage-1 exponent build  (no GPU used)\n"
"      =stage2plan          stage-2 pairing plan    (no GPU used)\n"
"      =bounds              B1/B2 choice model      (no GPU used)\n"
"      =engine              FFT/squaring engine vs a CPU reference\n"
"      =pm1                 P-1 against known factors\n"
"      =stage2              stage-2 engine vs a CPU reference\n"
"  --bench                time every FFT config that can hold an exponent\n"
"  --tune[=opts]          measure FFT configs and write tune.txt\n"
"\n"
"TYPICAL USE (a new exponent)\n"
"  1. put the exponent in config.txt, leave  fft = auto\n"
"  2. Mp_p-1_gpu.exe --tune noconfig     once -- targets that exponent and\n"
"                                          ends with a VERIFIED recommendation\n"
"  3. Mp_p-1_gpu.exe                     auto picks the best verified config\n"
"\n"
"  You do not need to copy the recommendation into config.txt: `auto` already\n"
"  takes tune.txt in cost order and verifies each entry before use, caching the\n"
"  verdict in fft-verified.txt so later runs pay nothing. Pin a spec only to\n"
"  override that choice.\n"
"\n"
"  'noconfig' is advised: the kernel-option search picks its baseline shape\n"
"  without checking it first, and if that shape is broken every measurement\n"
"  returns the failure sentinel, options get chosen from noise, and those\n"
"  options then corrupt the FFT timings. That has happened on this hardware.\n"
"\n"
"TUNING\n"
"  tune.txt is what makes FFT selection cost-based. Without it the engine\n"
"  falls back to \"first shape that fits\", which ignores the arithmetic type;\n"
"  configs of the SAME size were measured here to differ by up to 6x.\n"
"  Suggested -use kernel settings are appended to Mp_p-1_gpu-tune-config.txt (NOT\n"
"  to config.txt, which is this program's job file and a different format).\n"
"  That file is read back automatically on every normal run.\n"
"\n"
"  A bare --tune targets the exponent in config.txt and times both FP64 and\n"
"  NTT. Add 'noconfig' to skip the slower kernel-option search.\n"
"\n"
"  opts is a comma-separated list. Values must be PLAIN DIGITS -- they are\n"
"  parsed with stoull, so \"80M\" silently becomes 80, not 80000000.\n"
"      fp64               time the FP64 transforms\n"
"      ntt                time the integer-NTT transforms\n"
"                         (if you name NEITHER, both are timed. the upstream\n"
"                          tuner defaults to neither, i.e. it times nothing.)\n"
"      nofp32             skip transforms that use FP32\n"
"      quick=1..10        1 = slowest/most accurate, 10 = fastest. Default 7.\n"
"      minexp=N           only time FFTs usable for exponents >= N\n"
"      maxexp=N           only time FFTs usable for exponents <= N\n"
"      noconfig           skip the -use option search (times FFTs only)\n"
"  A full sweep takes roughly 30-90 min. To tune only what you need:\n"
"      --tune quick=10,minexp=80000000,maxexp=90000000\n"
"\n"
"OPTIONS\n"
"  -d <n>                 use GPU <n> (see --list-devices); default 0\n"
"  -fft <spec>            force an FFT config, e.g. 4:256:16:256\n"
"                         overrides both config.txt and tune.txt\n"
"  -h, --help             this text\n"
"\n"
"  --bench only:\n"
"  -E <exponent>          exponent to time at        (default 86599237)\n"
"  -n <iters>             iterations per config      (default 200)\n"
"  -slack <factor>        also try FFTs up to this multiple of the smallest\n"
"                         viable size                (default 2.1)\n"
"\n"
"STAGES\n"
"  Stage 1 finds q = 2kp+1 when k is entirely B1-smooth. Stage 2 also finds it\n"
"  when k is B1-smooth apart from ONE prime in (B1,B2] -- the most common\n"
"  near-miss by far, so it roughly doubles the yield for comparable work.\n"
"  config.txt:  stages = auto | both | 1     (or stage2 = auto | yes | no)\n"
"  With `stages = 1` B1 is optimised alone and comes out LARGER, since nothing\n"
"  else is left to catch a near-miss.\n"
"\n"
"Long runs write a single progress line to a terminal, or periodic timestamped\n"
"lines when redirected to a file. Stage 2 checkpoints as it goes and resumes\n"
"automatically; a resumed walk reproduces the accumulator bit for bit.\n");
}

// One line per PRIME. A gcd routinely carries several factors multiplied
// together -- every one whose k was smooth comes out of the same gcd -- so
// reporting it raw would be a composite masquerading as a factor.
// One JSON object per line, the shape PrimeNet accepts and the same convention
// Prime95 uses for its own results.txt -- so the file can be uploaded to
// mersenne.org's "Manual Results" page as-is.
//
//   {"status":"F","exponent":81679223,"worktype":"P-1","b1":2000000,
//    "b2":60000000,"factors":["..."],"program":{"name":"...","version":"1.0"},
//    "timestamp":"2026-07-28 17:29:54"}
//
// Only factors that are PRIME and verified to divide M_p are reported as
// factors: a gcd routinely carries several multiplied together, and submitting
// that product would be a composite masquerading as a factor. Anything that
// could not be split is written with status "C" and is not a submittable
// result -- it is recorded so the run is not silently lost.
void writeResultJson(const Config& cfg, const char* worktype, u64 b1, u64 b2,
                     const std::vector<FoundFactor>& factors, u32 seed) {
  FILE* f = fopen(cfg.resultsFile.c_str(), "a");
  if (!f) {
    printf("  WARNING: could not append to %s\n", cfg.resultsFile.c_str());
    return;
  }

  char stamp[32] = "";
  const time_t now = time(nullptr);
  struct tm utc;
  if (gmtime_s(&utc, &now) == 0) { strftime(stamp, sizeof(stamp), "%Y-%m-%d %H:%M:%S", &utc); }

  auto common = [&]() {
    fprintf(f, ",\"exponent\":%u,\"worktype\":\"%s\",\"b1\":%llu", cfg.exponent,
            worktype, (unsigned long long) b1);
    if (b2 > b1) { fprintf(f, ",\"b2\":%llu", (unsigned long long) b2); }
    if (seed) { fprintf(f, ",\"seed\":%u", seed); }
    fprintf(f, ",\"program\":{\"name\":\"%s\",\"version\":\"%s\"}", PROGRAM_NAME, PROGRAM_VERSION);
    fprintf(f, ",\"timestamp\":\"%s\"", stamp);
    if (!cfg.username.empty()) { fprintf(f, ",\"user\":\"%s\"", cfg.username.c_str()); }
    if (!cfg.computerName.empty()) { fprintf(f, ",\"computer\":\"%s\"", cfg.computerName.c_str()); }
  };

  std::vector<const FoundFactor*> good, bad;
  for (const FoundFactor& ff : factors) {
    ((ff.prime && ff.dividesMp) ? good : bad).push_back(&ff);
  }

  if (!good.empty()) {
    fprintf(f, "{\"status\":\"F\"");
    common();
    fprintf(f, ",\"factors\":[");
    for (size_t i = 0; i < good.size(); ++i) {
      fprintf(f, "%s\"%s\"", i ? "," : "", good[i]->value.dec().c_str());
    }
    fprintf(f, "]}\n");
  }
  for (const FoundFactor* ff : bad) {
    fprintf(f, "{\"status\":\"C\"");
    common();
    fprintf(f, ",\"composite\":\"%s\",\"note\":\"%s\"}\n", ff->value.dec().c_str(),
            ff->dividesMp ? "could not be split into primes"
                          : "DOES NOT DIVIDE M_p -- please report");
  }
  if (good.empty() && bad.empty()) {
    fprintf(f, "{\"status\":\"NF\"");
    common();
    fprintf(f, "}\n");
  }
  fclose(f);
}

// Console reporting. Separate from the results file on purpose: the file is
// machine-readable for submission, the console is for the person watching.
void reportFactors(const Config& cfg, const std::vector<FoundFactor>& factors,
                   u64 b1, u64 b2, const char* worktype, u32 seed) {
  for (const FoundFactor& ff : factors) {
    if (ff.prime && ff.dividesMp) {
      printf("  *** M%u has a factor: %s ***\n", cfg.exponent, ff.value.dec().c_str());
      printf("      %zu bits, k = %llu, verified 2^p == 1 (mod q)\n",
             ff.value.bits(), (unsigned long long) ff.k);
    } else if (ff.dividesMp) {
      printf("  *** M%u has a COMPOSITE divisor: %s ***\n",
             cfg.exponent, ff.value.dec().c_str());
      printf("      %zu bits; could not be split -- its factors have k larger\n"
             "      than the trial-division limit\n", ff.value.bits());
    } else {
      printf("  !!! %s does NOT divide M%u -- this is a bug, please report\n",
             ff.value.dec().c_str(), cfg.exponent);
    }
  }
  writeResultJson(cfg, worktype, b1, b2, factors, seed);
}

} // namespace

// Set from config once it is read; consulted by main() on every exit path,
// including the ones that return early or throw.
int gPauseMode = PAUSE_AUTO;

namespace {

// True when this process owns the console by itself, i.e. it was double-clicked
// rather than started from an existing shell. Closing such a window destroys
// the output before it can be read, which is the case worth pausing for.
bool launchedByDoubleClick() {
  DWORD pids[4];
  const DWORD n = GetConsoleProcessList(pids, 4);
  return n <= 1 && _isatty(_fileno(stdin)) != 0;
}

} // namespace

static int runMain(int argc, char** argv) {
  // Unbuffered throughout. Redirected to a file, stdout is fully buffered, so
  // anything after the last explicit fflush stays invisible -- which makes a
  // stall look like it happened at the last line that happened to be flushed.
  // The in-place progress line needs unbuffered output anyway.
  setvbuf(stdout, nullptr, _IONBF, 0);

  SetConsoleCtrlHandler(ctrlHandler, TRUE);

  try {
    std::string configPath = "config.txt";
    std::string fftSpec, selftest, tuneOpts;
    bool doSelftest = false, doBench = false, doTune = false, doList = false;
    bool doBounds = false;
    int deviceOverride = -1;
    u32 benchE = 86599237, benchIters = 200;
    double sizeSlack = 2.1;

    for (int i = 1; i < argc; ++i) {
      const std::string a = argv[i];
      if (a == "--config" && i + 1 < argc) { configPath = argv[++i]; }
      else if (a == "--bounds") { doBounds = true; }
      else if (a == "--selftest") { doSelftest = true; selftest = "all"; }
      else if (a.rfind("--selftest=", 0) == 0) { doSelftest = true; selftest = a.substr(11); }
      else if (a == "--list-devices") { doList = true; }
      else if (a == "--bench") { doBench = true; }
      // Tuning options are passed through to gpuowl's Tune, which reads them
      // from Args::tune as a comma-separated list. Accepts both
      //   --tune quick=10,minexp=80000000,maxexp=90000000
      //   --tune=quick=10,...
      // Recognised keys: quick=1..10 (10 is fastest/least accurate),
      // minexp=, maxexp= to restrict the exponent range that gets timed, and
      // "noconfig" to skip the -use option search.
      else if (a == "--tune") {
        doTune = true;
        if (i + 1 < argc && argv[i + 1][0] != '-') { tuneOpts = argv[++i]; }
      }
      else if (a.rfind("--tune=", 0) == 0) { doTune = true; tuneOpts = a.substr(7); }
      else if (a == "-d" && i + 1 < argc) { deviceOverride = atoi(argv[++i]); }
      else if (a == "-fft" && i + 1 < argc) { fftSpec = argv[++i]; }
      else if (a == "-E" && i + 1 < argc) { benchE = strtoul(argv[++i], nullptr, 10); }
      else if (a == "-n" && i + 1 < argc) { benchIters = strtoul(argv[++i], nullptr, 10); }
      else if (a == "-slack" && i + 1 < argc) { sizeSlack = atof(argv[++i]); }
      else if (a == "-h" || a == "--help") { banner(); usage(); return 0; }
      else { banner(); printf("unrecognised option '%s'\n\n", a.c_str()); usage(); return 2; }
    }

    banner();

    // ---- load the job -----------------------------------------------------
    // Read before ANY early return, so pause_on_exit governs every exit path --
    // including --list-devices, the self-tests and a config error itself, which
    // are exactly the short paths a double-clicked window would flash past.
    // Also read for --tune: a tune that does not cover the exponent you
    // actually run is wasted effort, so the exponent is worth knowing there.
    Config cfg;
    bool haveConfig = false;
    {
      std::string err;
      haveConfig = loadConfig(configPath, cfg, err);
      if (haveConfig) { gPauseMode = cfg.pauseMode; }
      if (!haveConfig && !doSelftest && !doBench && !doTune && !doList) {
        printf("config error: %s\n", err.c_str());
        return 2;
      }
      if (haveConfig) {
        if (!fftSpec.empty()) { cfg.fftSpec = fftSpec; }
        if (deviceOverride >= 0) { cfg.device = deviceOverride; }
      }
    }

    if (doList) {
      auto devices = getAllDeviceIDs();
      printf("OpenCL GPU devices:\n");
      for (u32 i = 0; i < devices.size(); ++i) {
        printf("  device %u: %s  [driver %s]\n", i,
               getDeviceName(devices[i]).c_str(), getDriverVersion(devices[i]).c_str());
      }
      return 0;
    }

    // These self-tests need no GPU, so handle them before touching a device.
    if (doSelftest && selftest == "gcd") { return runBigIntTests(); }
    if (doSelftest && selftest == "exponent") { return runExponentTest(); }
    if (doSelftest && selftest == "stage2plan") { return runStage2PlanTests(); }
    if (doSelftest && selftest == "bounds") { return runBoundsTests(); }

    // --bounds: show the surface around the automatic choice, for the job in
    // config.txt. No GPU, so the pairing shape comes from the configured values
    // or a mid-range default rather than from a device query.
    if (doBounds) {
      if (!haveConfig) { printf("--bounds needs a valid config.txt\n"); return 2; }
      CostModel cm;
      cm.gcdIters = gcdIterCost(cfg.exponent);
      cm.tolerance = cfg.boundsTolerance;
      // The SAME shape the real run would pick, queried from the device, so the
      // bounds shown here are the bounds you will actually get. The one
      // difference: with no GPU up, the transform size comes from the smallest
      // shape that fits rather than from tune.txt's choice.
      const FFTConfig f = smallestFittingFFT(cfg.exponent);
      const u64 residueBytes = u64(f.size()) * sizeof(Word);
      const Stage2Shape s = pickStage2Shape(cfg, deviceOverride, residueBytes);
      cm.mulsPerPrime = buildStage2Plan(1000000, 4000000, s.d, s.w).mulsPerPrime();
      printf("stage-2 shape D=%u w=%u (%u T-buffers of %.1f MB = %.2f GB)\n",
             s.d, s.w, stage2NumJ(s.d, s.w), double(residueBytes) / (1 << 20),
             double(stage2NumJ(s.d, s.w)) * double(residueBytes) / (1u << 30));
      printBoundsSurface(cfg.exponent, cfg.factoredTo, cfg.bias, cm);
      return 0;
    }

    // ---- bring up the device ---------------------------------------------
    Args args;
    args.device = (deviceOverride >= 0) ? deviceOverride
                                        : (cfg.device >= 0 ? cfg.device : 0);
    // Command line only. config.txt's fft= describes the JOB's exponent and
    // must not leak into --selftest (whose exponents are ~1e6, where a pinned
    // 2M-word transform is far too large) or into --tune.
    args.fftSpec = fftSpec;

    // gpuowl's Tune defaults time_FFTs and time_NTTs BOTH to 0, so a bare
    // "--tune" skips every shape and writes no tune.txt at all -- it looks like
    // a silent failure. Neither selector is a sensible default for a user who
    // just asked to tune, so select both unless one was named explicitly.
    // Apply the kernel options a previous --tune discovered.
    //
    // The tune's option search measures things like MODM31 and TABMUL_CHAIN61
    // and writes the winners as "-use KEY=VAL" lines. Args::readConfig() parses
    // exactly that format -- but upstream nothing ever calls it here, so those
    // settings were being computed and then discarded. Skipped during --tune
    // itself, so a tune always measures from the stock baseline and stays
    // reproducible.
    // The old name is still read if present: the file holds measured kernel
    // settings that took 30-90 minutes to produce, and silently ignoring one
    // left over from before the rename would quietly lose that tuning.
    if (!doTune) {
      const char* tuneCfg = fs::exists("Mp_p-1_gpu-tune-config.txt") ? "Mp_p-1_gpu-tune-config.txt"
                          : fs::exists("gpuowl-tune-config.txt")    ? "gpuowl-tune-config.txt"
                          : nullptr;
      if (tuneCfg) {
        args.readConfig(tuneCfg);
        if (!args.flags.empty()) {
          printf("Applied %u tuned kernel option(s) from %s\n", u32(args.flags.size()), tuneCfg);
        }
      }
    }

    args.tune = tuneOpts;
    if (doTune) {
      if (tuneOpts.find("fp64") == std::string::npos
          && tuneOpts.find("ntt") == std::string::npos) {
        args.tune = args.tune.empty() ? "fp64,ntt" : args.tune + ",fp64,ntt";
        printf("Tuning both FP64 and NTT transforms (pass 'fp64' or 'ntt' to pick one).\n");
      }
      // gpuowl's default tune window is 75M-350M, which has nothing to do with
      // the job being run: tuning outside the exponent's range yields a
      // tune.txt that bestFit then cannot use for it. Default to the configured
      // exponent so a bare "--tune" produces something actually applicable.
      if (haveConfig && cfg.exponent
          && tuneOpts.find("minexp") == std::string::npos
          && tuneOpts.find("maxexp") == std::string::npos) {
        args.tune += ",minexp=" + std::to_string(cfg.exponent)
                   + ",maxexp=" + std::to_string(cfg.exponent);
        printf("Tuning for M%u (from %s); override with minexp=/maxexp=.\n",
               cfg.exponent, configPath.c_str());
      }
      printf("\n");
    }
    args.setDefaults();

    cl_device_id devId = getDevice(args.device);
    Context context{devId};
    Queue queue{context, false};
    TrigBufCache bufCache{&context};
    Background background;
    GpuCommon shared{&args, &bufCache, &background};

    printf("Device: %s\n\n", getDeviceName(devId).c_str());

    if (doTune) {
      Tune tune{&queue, shared};
      tune.tune();
      // Never recommend tune's winner unverified -- see recommendFFT().
      if (haveConfig && cfg.exponent) { recommendFFT(shared, &queue, cfg.exponent); }
      else { printf("\n  (no exponent in config.txt, so tune.txt was not verified)\n"); }
      return 0;
    }

    if (doBench) {
      extern void benchmarkFFT(GpuCommon, Queue*, u32, u32, double, const std::string&);
      benchmarkFFT(shared, &queue, benchE, benchIters, sizeSlack, args.fftSpec);
      return 0;
    }

    if (doSelftest) {
      int rc = 0;
      if (selftest == "all")                          { rc |= runBigIntTests(); }
      if (selftest == "all")                          { rc |= runStage2PlanTests(); }
      if (selftest == "all")                          { rc |= runBoundsTests(); }
      if (selftest == "all" || selftest == "engine")  { rc |= runEngineTests(shared, &queue, args.fftSpec, selftest != "all"); }
      if (selftest == "all" || selftest == "pm1")     { rc |= runPM1Tests(shared, &queue, args.fftSpec); }
      if (selftest == "all" || selftest == "extend")  { rc |= runExtendTests(shared, &queue, args.fftSpec); }
      if (selftest == "all" || selftest == "pp1")     { rc |= runPP1Tests(shared, &queue, args.fftSpec); }
      if (selftest == "all" || selftest == "stage2")  { rc |= runStage2Tests(shared, &queue, args.fftSpec); }
      if (selftest != "all" && selftest != "engine" && selftest != "pm1"
          && selftest != "extend" && selftest != "pp1" && selftest != "stage2") {
        printf("unknown selftest '%s'\n"
               "  no GPU: gcd, exponent, stage2plan, bounds\n"
               "  GPU:    engine, pm1, extend, pp1, stage2\n"
               "  omit =which to run them all\n",
               selftest.c_str());
        return 2;
      }
      return rc;
    }

    // Phase 3 (the gcd) is CPU-only; this is the knob that speeds it up.
    setGcdThreads(cfg.gcdThreads);
    printf("  gcd worker threads: %s\n",
           cfg.gcdThreads ? std::to_string(cfg.gcdThreads).c_str() : "auto (all cores)");

    // ---- FFT, then bounds -------------------------------------------------
    // Order matters. The transform is chosen first because it fixes the size of
    // one residue, which fixes how many T_j buffers fit, which fixes the stage-2
    // cost per prime -- and that is an input to the bounds model. None of it
    // depends on B1/B2, so there is no circularity, but guessing the residue
    // size instead (as an earlier version did) was wrong by 2x in each direction
    // at once.
    const std::string jobSpec = !fftSpec.empty() ? fftSpec : cfg.fftSpec;
    FFTConfig fft = chooseVerifiedFFT(shared, &queue, cfg.exponent,
                                      jobSpec, cfg.verifyFft);
    printf("  FFT %s (%llu words, %.2f bits/word)\n", fft.spec().c_str(),
           (unsigned long long) fft.size(), double(cfg.exponent) / fft.size());

    auto gpu = Gpu::make(&queue, cfg.exponent, shared, fft, {}, false);

    CostModel cost;
    cost.gcdIters = gcdIterCost(cfg.exponent);
    cost.tolerance = cfg.boundsTolerance;

    // Ask the Gpu, do not compute it here: stage 2 holds its T_j table in
    // middle-transform form, which is 1.5x a Word buffer for these NTTs and also
    // depends on the -use INPLACE/PAD options. Budgeting with the Word size
    // would over-commit by half and could exhaust VRAM mid-walk.
    const u64 residueBytes = gpu->stage2BufferBytes();
    Stage2Shape shape = pickStage2Shape(cfg, deviceOverride, residueBytes);
    cost.mulsPerPrime = buildStage2Plan(1000000, 4000000, shape.d, shape.w).mulsPerPrime();

    const bool wantStage2 = cfg.doPM1 && cfg.stage2Mode != STAGE2_OFF;
    Bounds bounds = chooseBounds(cfg.exponent, cfg.factoredTo, cfg.bias, cost,
                                 cfg.b1, cfg.b2, wantStage2);
    // stages = both forces a stage 2 even where the model would decline one.
    if (wantStage2 && cfg.stage2Mode == STAGE2_ON && bounds.b2 <= bounds.b1) {
      bounds = chooseBounds(cfg.exponent, cfg.factoredTo, cfg.bias, cost,
                            bounds.b1, cfg.b2 ? cfg.b2 : bounds.b1 * 30, true);
    }

    u64 b1 = bounds.b1;
    const bool runStage2 = wantStage2 && bounds.b2 > bounds.b1;
    gPhaseTotal = runStage2 ? 5 : 3;

    printf("M%u, trial-factored to %u bits, bias %.1f\n", cfg.exponent, cfg.factoredTo, cfg.bias);
    printf("  B1 = %llu%s\n", (unsigned long long) b1, cfg.b1 ? "" : " (auto)");
    if (runStage2) {
      printf("  B2 = %llu%s   pairing D=%u w=%u, %.3f muls/prime\n",
             (unsigned long long) bounds.b2, cfg.b2 ? "" : " (auto)",
             shape.d, shape.w, cost.mulsPerPrime);
      printf("  stage-2 table: %u buffers x %.1f MB = %.2f GB of GPU memory\n",
             stage2NumJ(shape.d, shape.w), double(residueBytes) / (1 << 20),
             double(stage2NumJ(shape.d, shape.w)) * double(residueBytes) / (1u << 30));
    } else if (cfg.stage2Mode == STAGE2_OFF) {
      printf("  stage 2 disabled (stages = 1); B1 is optimised alone, which makes\n"
             "  it larger than it would be with a stage 2 to catch near-misses\n");
    } else {
      printf("  stage 2 not worth running at this exponent and bias -- raise `bias`\n"
             "  if a factor is worth more to you than one PRP test\n");
    }
    if (runStage2) {
      printf("  estimated success %.3f%%  (stage 1 %.3f%% + stage 2 %.3f%%)\n\n",
             bounds.prob() * 100, bounds.probStage1 * 100, bounds.probStage2 * 100);
    } else {
      printf("  estimated success %.3f%%\n\n", bounds.probStage1 * 100);
    }

    // ---- P+1, if asked for -------------------------------------------------
    // Run before P-1 only when P-1 is not also requested; otherwise P-1 first,
    // since it is half the cost per bound and more likely to succeed.
    bool anyFactor = false;
    if (cfg.doPP1) {
      u32 seedIx = 0;
      for (u32 seed : cfg.pp1Seeds) {
        if (gInterrupted.load()) { break; }
        // Index AND value: the default seeds are 3, 5, 7 and there are 3 of
        // them, so "seed 3 (of 3)" read as though it were the last one.
        printf("\n  P+1 attempt %u of %u (seed %u)\n", ++seedIx,
               u32(cfg.pp1Seeds.size()), seed);
        PP1Result pr = runPP1Stage1(*gpu, cfg, b1, seed, true);
        if (pr.interrupted) {
          printf("\n  interrupted; P+1 progress for seed %u is checkpointed.\n", seed);
          return 1;
        }
        if (pr.foundFactor) {
          anyFactor = true;
          reportFactors(cfg, pr.factors, b1, b1, "P+1", seed);
          break;                      // no point trying further seeds
        }
        printf("  P+1 seed %u: no factor\n", seed);
      }
      if (!anyFactor) { writeResultJson(cfg, "P+1", b1, b1, {}, 0); }
      if (!cfg.doPM1) {
        printf("\n  appended to %s\n", cfg.resultsFile.c_str());
        return 0;
      }
    }

    PM1Result r = runPM1Stage1(*gpu, cfg, b1, true);

    if (r.interrupted) {
      printf("\n  interrupted before stage 1 completed; nothing written.\n");
      return 1;
    }

    printf("\n");
    if (r.foundFactor) {
      if (r.factors.size() > 1) {
        printf("  the gcd was a product of %zu factors (every factor with a\n"
               "  B1-smooth k comes out of the same gcd):\n\n", r.factors.size());
      }
      for (const FoundFactor& ff : r.factors) {
        if (ff.prime && ff.dividesMp) {
          printf("  *** M%u has a factor: %s ***\n", cfg.exponent, ff.value.dec().c_str());
          printf("      %zu bits, k = %llu, verified 2^p == 1 (mod q)\n",
                 ff.value.bits(), (unsigned long long) ff.k);
        } else if (ff.dividesMp) {
          printf("  *** M%u has a COMPOSITE divisor: %s ***\n",
                 cfg.exponent, ff.value.dec().c_str());
          printf("      %zu bits; could not be split into primes -- its factors\n"
                 "      have k larger than the trial-division limit\n", ff.value.bits());
        } else {
          printf("  !!! %s does NOT divide M%u -- this is a bug, please report\n",
                 ff.value.dec().c_str(), cfg.exponent);
        }
      }
    } else {
      printf("  M%u: no factor found with B1 = %llu\n",
             cfg.exponent, (unsigned long long) r.b1Used);
    }
    // ONE result per job, not one per stage. If stage 2 is going to run, the
    // result is written after it with both bounds; emitting a stage-1-only "NF"
    // here as well would double-report the exponent to PrimeNet and understate
    // the work actually done. A factor found in stage 1 ends the job, so that
    // case is reported here and stage 2 never runs.
    if (r.foundFactor || !runStage2) {
      writeResultJson(cfg, "P-1", r.b1Used, r.b1Used, r.factors, 0);
    }
    printf("  appended to %s\n", cfg.resultsFile.c_str());

    // ---- stage 2 ----------------------------------------------------------
    // Only when stage 1 came up empty. A factor already in hand makes the whole
    // second stage wasted work.
    if (runStage2 && !r.foundFactor && !gInterrupted.load()) {
      printf("\n");
      const Stage2Plan plan = buildStage2Plan(b1, bounds.b2, shape.d, shape.w);
      PM1Stage2Result s2 = runPM1Stage2(*gpu, cfg, r.residue, plan, true);

      printf("\n");
      if (s2.interrupted) {
        printf("  interrupted during stage 2; resume by running again.\n");
        return 1;
      }
      if (s2.foundFactor) {
        reportFactors(cfg, s2.factors, b1, bounds.b2, "P-1", 0);
      } else {
        printf("  M%u: no factor found in stage 2 either (B1=%llu, B2=%llu)\n",
               cfg.exponent, (unsigned long long) b1, (unsigned long long) bounds.b2);
        writeResultJson(cfg, "P-1", b1, bounds.b2, {}, 0);
      }
      printf("  appended to %s\n", cfg.resultsFile.c_str());
    }
    return 0;

  } catch (const char* s) {
    printf("\nFAILED: %s\n", s);
    return 2;
  } catch (const std::exception& e) {
    printf("\nFAILED: %s\n", e.what());
    return 2;
  }
}

int main(int argc, char** argv) {
  const int rc = runMain(argc, argv);

  // Wrapping runMain() rather than pausing inside it means every exit path is
  // covered -- early returns, config errors and exceptions included. Those are
  // exactly the cases a double-clicked window would otherwise flash and close.
  const bool hold = (gPauseMode == PAUSE_ALWAYS) ||
                    (gPauseMode == PAUSE_AUTO && launchedByDoubleClick());
  if (hold) {
    printf("\n  Press Enter to close this window . . . ");
    fflush(stdout);
    int c;
    while ((c = getchar()) != '\n' && c != EOF) {}
  }
  return rc;
}
