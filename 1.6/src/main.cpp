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
#include "Worktodo.h"
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
#include <future>
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
int runB2ExtendTests(GpuCommon shared, Queue* q, const std::string& fftSpec);
int runPp1Stage2Tests(GpuCommon shared, Queue* q, const std::string& fftSpec);
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
// b1Hint sizes the pairing window (w*D/2 must not exceed whatever B1 stage 2
// actually runs against, see Stage2Plan.cpp) -- it is a HINT because the real
// B1 isn't known yet the first time this is called (chooseBounds needs the
// shape's mulsPerPrime first). Callers that later learn the real B1 undershot
// the hint call this again with the real value; see main.cpp's job driver.
Stage2Shape pickStage2Shape(const Config& cfg, int deviceOverride, u64 residueBytes, u64 b1Hint) {
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

  Stage2Shape s = chooseStage2Shape(budget, residueBytes, b1Hint);
  if (!s.d) { s = {210, 1}; }
  if (!s.w) { s.w = 1; }
  return s;
}

void banner() {
  printf("=======================================================================\n");
  printf(" %s %s  --  GPU P-1 / P+1 factoring of M_p = 2^p - 1\n",
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
"      =worktodo            worktodo.txt parsing    (no GPU used)\n"
"      =engine              FFT/squaring engine vs a CPU reference\n"
"      =pm1                 P-1 against known factors\n"
"      =extend              raising B1 on a completed stage 1\n"
"      =pp1                 P+1 Lucas ladder vs a CPU reference\n"
"      =stage2              stage-2 engine vs a CPU reference\n"
"      =b2extend            raising B2 on a completed stage 2\n"
"      =pp1stage2           P+1 stage 2\n"
"  --bench                time FFT configs near the smallest that can hold an\n"
"                         exponent (see -slack below)\n"
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
"  'noconfig' is advised for a first run: the kernel-option search verifies its\n"
"  baseline shape before trusting ~50 measurements against it (substituting a\n"
"  working shape, or skipping the search cleanly, if the baseline is broken --\n"
"  this has happened on this hardware) but it is still the slow part of\n"
"  --tune and rarely changes the FFT choice that matters most.\n"
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
"      minexp=N           skip FFTs too small to hold exponent N\n"
"      maxexp=N           skip FFTs much larger than exponent N (roughly\n"
"                         1.2x-2x, not an exact cutoff; if minexp==maxexp\n"
"                         finds nothing, it retries with a 10x ceiling)\n"
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
"  P+1 (method = pp1 | both) targets a different quantity, q+1 rather than\n"
"  q-1, so it catches factors P-1 cannot and vice versa; it has its own B1\n"
"  model but still shares P-1's B2 and stage-2 pairing shape.\n"
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

// Console AND log reporting -- this is the single most important line an
// unattended run produces, so it goes through log() rather than printf(),
// unlike most of this file's routine progress output. Separate from the
// results file on purpose: that one is machine-readable for submission,
// this is for the person (or the log) watching.
void reportFactors(const Config& cfg, const std::vector<FoundFactor>& factors,
                   u64 b1, u64 b2, const char* worktype, u32 seed) {
  if (factors.size() > 1) {
    log("  the gcd was a product of %zu factors (every factor with a\n"
        "  B1-smooth k comes out of the same gcd):\n\n", factors.size());
  }
  for (const FoundFactor& ff : factors) {
    if (ff.prime && ff.dividesMp) {
      log("  *** M%u has a factor: %s ***\n", cfg.exponent, ff.value.dec().c_str());
      log("      %zu bits, k = %llu, verified 2^p == 1 (mod q)\n",
          ff.value.bits(), (unsigned long long) ff.k);
    } else if (ff.dividesMp) {
      log("  *** M%u has a COMPOSITE divisor: %s ***\n",
          cfg.exponent, ff.value.dec().c_str());
      log("      %zu bits; could not be split -- its factors have k larger\n"
          "      than the trial-division limit\n", ff.value.bits());
    } else {
      log("  !!! %s does NOT divide M%u -- this is a bug, please report\n",
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

// Runs one exponent's full job -- FFT/bounds selection through P+1/P-1
// stage 1 and (if warranted) stage 2 -- and returns 0 (done, whether or not
// a factor turned up), 1 (interrupted; nothing more to do until relaunched),
// or throws on failure (caught by runMain's own try/catch, which turns it
// into the usual "FAILED: ..." message and return 2). Called once per
// worktodo.txt entry by runMain's queue loop; cfg.exponent/cfg.factoredTo
// must already be set for the entry being run.
static int runOneJob(Config& cfg, GpuCommon shared, Queue& queue,
                     const std::string& fftSpec, int deviceOverride) {
  // Phase 3 (the gcd) is CPU-only; this is the knob that speeds it up.
  setGcdThreads(cfg.gcdThreads);
  printf("  gcd worker threads: %s\n",
         cfg.gcdThreads ? std::to_string(cfg.gcdThreads).c_str() : "auto (all cores)");

  // ---- FFT, then bounds ---------------------------------------------------
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
  Stage2Shape shape = pickStage2Shape(cfg, deviceOverride, residueBytes, cfg.b1 ? cfg.b1 : 100000);
  cost.mulsPerPrime = buildStage2Plan(1000000, 4000000, shape.d, shape.w).mulsPerPrime();

  // Whether ANY active method wants a stage 2 considered at all -- this has
  // to be method-agnostic, not just P-1's own wantStage2 below: a pure
  // `method = pp1` run has cfg.doPM1 == false, and chooseBounds's
  // allowStage2 argument decides whether bounds.b2 gets a non-trivial value
  // in the first place. Gating that on P-1 alone would silently keep
  // bounds.b2 == bounds.b1 forever for a P+1-only job, no matter what `b2`
  // or `stages` say.
  const bool wantAnyStage2 = cfg.stage2Mode != STAGE2_OFF && (cfg.doPM1 || cfg.doPP1);
  const bool wantStage2 = cfg.doPM1 && cfg.stage2Mode != STAGE2_OFF;
  Bounds bounds = chooseBounds(cfg.exponent, cfg.factoredTo, cfg.bias, cost,
                               cfg.b1, cfg.b2, wantAnyStage2);
  // stages = both forces a stage 2 even where the model would decline one.
  if (wantAnyStage2 && cfg.stage2Mode == STAGE2_ON && bounds.b2 <= bounds.b1) {
    bounds = chooseBounds(cfg.exponent, cfg.factoredTo, cfg.bias, cost,
                          bounds.b1, cfg.b2 ? cfg.b2 : bounds.b1 * 30, true);
  }

  u64 b1 = bounds.b1;
  const bool runStage2 = wantStage2 && bounds.b2 > bounds.b1;

  // P+1's own bounds: q+1's smoothness target is different from P-1's q-1
  // (see pp1Prob's own comment in Bounds.cpp), so reusing P-1's B1 or B2
  // here would optimise the wrong quantity. The pairing SHAPE (D, w,
  // T-table sizing) is still shared -- that's a GPU-memory budget decision,
  // not a cost-model gap, and splitting it would double VRAM use for both
  // methods running at once. `stages = 1` still means "no stage 2 for
  // either method," a single consistent switch (wantAnyStage2 above).
  const Bounds pp1Bounds = choosePP1Bounds(cfg.exponent, cfg.factoredTo, cfg.bias, cost,
                                           u32(cfg.pp1Seeds.size()), cfg.b1, cfg.b2,
                                           wantAnyStage2);
  const u64 pp1B1 = pp1Bounds.b1;
  const bool wantPp1Stage2 = cfg.doPP1 && cfg.stage2Mode != STAGE2_OFF && pp1Bounds.b2 > pp1B1;

  // The shape above was sized against a GUESS (cfg.b1, or 100000 when b1 is
  // auto) because the real B1s aren't known until chooseBounds/choosePP1Bounds
  // run, and those need the shape's mulsPerPrime first. A small exponent or
  // a low bias/factoredTo can auto-pick a real B1 well under that guess --
  // and buildStage2Plan refuses to run stage 2 when w*D/2 exceeds B1 (every
  // prime in range must exceed it, see Stage2Plan.cpp). Pre-existing bug,
  // not new here: confirmed the untouched v1.2 binary crashes the same way
  // at M859433 TF20 bias1 b1=auto. Shrink the shape to fit whichever real B1
  // is smallest among the methods that will actually run stage 2 against it
  // -- unless the shape was pinned explicitly (stage2_d in config.txt),
  // which is the user's own choice, not this function's to override.
  if (!cfg.stage2D && (runStage2 || wantPp1Stage2)) {
    u64 minStageB1 = ~0ull;
    if (runStage2) { minStageB1 = std::min(minStageB1, b1); }
    if (wantPp1Stage2) { minStageB1 = std::min(minStageB1, pp1B1); }
    if (u64(shape.d) * shape.w / 2 > minStageB1) {
      shape = pickStage2Shape(cfg, deviceOverride, residueBytes, minStageB1);
      cost.mulsPerPrime = buildStage2Plan(1000000, 4000000, shape.d, shape.w).mulsPerPrime();
    }
  }

  const bool anyRunStage2 = runStage2 || wantPp1Stage2;
  printf("  M%u, trial-factored to %u bits, bias %.1f\n", cfg.exponent, cfg.factoredTo, cfg.bias);
  if (cfg.doPM1 && cfg.doPP1) {
    // Two different B1s now -- label which is which, the same way the
    // per-method sections below are already labelled.
    printf("  P-1: B1 = %llu%s\n", (unsigned long long) b1, cfg.b1 ? "" : " (auto)");
    printf("  P+1: B1 = %llu%s (%u seeds)\n", (unsigned long long) pp1B1,
           cfg.b1 ? "" : " (auto)", u32(cfg.pp1Seeds.size()));
  } else if (cfg.doPM1) {
    printf("  B1 = %llu%s\n", (unsigned long long) b1, cfg.b1 ? "" : " (auto)");
  } else {
    printf("  B1 = %llu%s (%u seeds)\n", (unsigned long long) pp1B1,
           cfg.b1 ? "" : " (auto)", u32(cfg.pp1Seeds.size()));
  }
  if (anyRunStage2) {
    if (runStage2 && wantPp1Stage2) {
      // Two independently-chosen B2s now -- label which is which, same
      // convention as the B1 split above. The pairing shape (D, w, table
      // size) is a shared GPU-memory decision, so it gets its own line
      // rather than repeating per method.
      printf("  P-1: B2 = %llu%s\n", (unsigned long long) bounds.b2, cfg.b2 ? "" : " (auto)");
      printf("  P+1: B2 = %llu%s\n", (unsigned long long) pp1Bounds.b2, cfg.b2 ? "" : " (auto)");
      printf("  pairing D=%u w=%u, %.3f muls/prime\n", shape.d, shape.w, cost.mulsPerPrime);
    } else {
      const u64 b2Shown = runStage2 ? bounds.b2 : pp1Bounds.b2;
      printf("  B2 = %llu%s   pairing D=%u w=%u, %.3f muls/prime\n",
             (unsigned long long) b2Shown, cfg.b2 ? "" : " (auto)",
             shape.d, shape.w, cost.mulsPerPrime);
    }
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
  if (cfg.doPM1) {
    if (runStage2) {
      printf("  estimated P-1 success %.3f%%  (stage 1 %.3f%% + stage 2 %.3f%%)\n",
             bounds.prob() * 100, bounds.probStage1 * 100, bounds.probStage2 * 100);
    } else {
      printf("  estimated P-1 success %.3f%%\n", bounds.probStage1 * 100);
    }
  }
  if (cfg.doPP1) {
    if (wantPp1Stage2) {
      printf("  estimated P+1 success %.3f%%  (stage 1 %.3f%% + stage 2 %.3f%%, %u seeds)\n",
             pp1Bounds.prob() * 100, pp1Bounds.probStage1 * 100, pp1Bounds.probStage2 * 100,
             u32(cfg.pp1Seeds.size()));
    } else {
      printf("  estimated P+1 success %.3f%%  (%u seeds)\n",
             pp1Bounds.probStage1 * 100, u32(cfg.pp1Seeds.size()));
    }
  }
  printf("\n");

  // ---- P+1, if asked for ---------------------------------------------------
  // Run before P-1 only when P-1 is not also requested; otherwise P-1 first,
  // since it is half the cost per bound and more likely to succeed.
  bool anyFactor = false;
  if (cfg.doPP1) {
    gPhaseTotal = wantPp1Stage2 ? 5 : 3;
    u32 seedIx = 0;
    for (u32 seed : cfg.pp1Seeds) {
      if (gInterrupted.load()) { break; }
      // Index AND value: the default seeds are 3, 5, 7 and there are 3 of
      // them, so "seed 3 (of 3)" read as though it were the last one.
      printf("\n  P+1 attempt %u of %u (seed %u)\n", ++seedIx,
             u32(cfg.pp1Seeds.size()), seed);
      PP1Result pr = runPP1Stage1(*gpu, cfg, pp1B1, seed, true);
      if (pr.interrupted) {
        log("\n  interrupted; P+1 progress for seed %u is checkpointed.\n", seed);
        return 1;
      }
      if (pr.foundFactor) {
        anyFactor = true;
        reportFactors(cfg, pr.factors, pp1B1, pp1B1, "P+1", seed);
        break;                      // no point trying further seeds
      }
      log("  P+1 seed %u: no factor\n", seed);

      // Only when stage 1 came up empty for this seed -- a factor already in
      // hand makes the second stage wasted work, exactly like P-1's own rule.
      if (wantPp1Stage2 && !gInterrupted.load()) {
        PP1Stage2Result s2 = runPP1Stage2(*gpu, cfg, pr.residue, seed, pp1B1,
                                          pp1Bounds.b2, shape.d, shape.w, true);
        if (s2.interrupted) {
          log("\n  interrupted during P+1 stage 2 for seed %u; resume by running again.\n", seed);
          return 1;
        }
        if (s2.foundFactor) {
          anyFactor = true;
          reportFactors(cfg, s2.factors, pp1B1, pp1Bounds.b2, "P+1", seed);
          break;
        }
        log("  P+1 seed %u: no factor in stage 2 either (B2=%llu)\n",
            seed, (unsigned long long) pp1Bounds.b2);
      }
    }
    if (!anyFactor) {
      writeResultJson(cfg, "P+1", pp1B1, wantPp1Stage2 ? pp1Bounds.b2 : pp1B1, {}, 0);
    }
    log("\n  appended to %s\n", cfg.resultsFile.c_str());
    if (!cfg.doPM1) {
      return 0;
    }
  }

  // Only label this when it follows a P+1 section above -- a pure `method =
  // pm1` job has nothing before it to separate from.
  if (cfg.doPP1) {
    printf("\n  P-1 attempt\n");
  }
  gPhaseTotal = runStage2 ? 5 : 3;

  // The stage-1 gcd is CPU-only and used to run with the GPU sitting idle,
  // which at 100M-digit exponents costs more wall time than the stage-1
  // squarings that produced its input. Stage 2 does not actually need the
  // gcd's ANSWER -- it continues from the stage-1 residue, which is ready the
  // moment stage 1 ends -- it only needs to know whether a factor turned up,
  // because a factor in hand makes stage 2 pointless. So when a stage 2 is
  // coming, run the gcd on another thread and start stage 2 immediately,
  // rather than doing them one after the other.
  //
  // That is a bet, and the odds are on the record: the bounds model prints
  // its own stage-1 success estimate, typically a few percent. In the ~97% of
  // jobs where the gcd comes up empty the overlap is pure profit; in the rest
  // the stage-2 work done so far is discarded -- but a factor was found,
  // which ends the job successfully anyway. Losing that race is the good
  // outcome.
  const bool overlapGcd = runStage2;
  PM1Result r = runPM1Stage1(*gpu, cfg, b1, true, !overlapGcd);

  if (r.interrupted) {
    log("\n  interrupted before stage 1 completed; nothing written.\n");
    return 1;
  }

  // Started here, joined after stage 2 (or immediately below when there is no
  // stage 2 to overlap with). Silent while it runs: its progress line and
  // stage 2's would otherwise interleave on the same terminal.
  std::future<void> gcdTask;
  if (overlapGcd && !gInterrupted.load()) {
    printf("\n  [%u/5 gcd CPU] gcd(x-1, M_p) running alongside stage 2 --"
           " reported when both finish\n", 3);
    fflush(stdout);
    gcdTask = std::async(std::launch::async,
                         [&r, &cfg] { finishStage1Gcd(r, cfg.exponent, false, /*announce=*/false); });
  }

  // Reported only once the gcd has actually answered. When it is overlapped
  // its result is not in yet, so this waits until after stage 2 -- announcing
  // "no factor found with B1" here would be a guess, not a result.
  if (!overlapGcd) {
    printf("\n");
    if (r.foundFactor) {
      reportFactors(cfg, r.factors, r.b1Used, r.b1Used, "P-1", 0);
      log("  appended to %s\n", cfg.resultsFile.c_str());
    } else {
      log("  M%u: no factor found with B1 = %llu\n",
          cfg.exponent, (unsigned long long) r.b1Used);
      // ONE result per job, not one per stage. If stage 2 is going to run, the
      // result is written after it with both bounds; emitting a stage-1-only
      // "NF" here as well would double-report the exponent to PrimeNet and
      // understate the work actually done. A factor found in stage 1 ends the
      // job, so that case is always reported above regardless of runStage2 --
      // and the "appended" line only prints when a write actually happened.
      if (!runStage2) {
        writeResultJson(cfg, "P-1", r.b1Used, r.b1Used, {}, 0);
        log("  appended to %s\n", cfg.resultsFile.c_str());
      }
    }
  }

  // ---- stage 2 --------------------------------------------------------------
  // Only when stage 1 came up empty -- a factor already in hand makes the whole
  // second stage wasted work. When the gcd is overlapped, "came up empty" is
  // not known yet and this runs speculatively; see the bet described above.
  bool ranStage2 = false;
  PM1Stage2Result s2;
  if (runStage2 && !r.foundFactor && !gInterrupted.load()) {
    printf("\n");
    // The plan is built inside: which range actually has to be walked depends
    // on whether a completed stage 2 for a smaller B2 is on disk.
    s2 = runPM1Stage2(*gpu, cfg, r.residue, b1, bounds.b2,
                      shape.d, shape.w, true);
    ranStage2 = true;
  }

  // Collect the overlapped gcd. Its answer decides which of the two results
  // above is the one worth reporting, so nothing is announced until it lands.
  if (gcdTask.valid()) {
    Timer waitTimer;
    gcdTask.get();
    printf("\n  [%u/5 gcd CPU] done in %s", 3, fmtDuration(r.gcdSecs).c_str());
    // Only interesting when the gcd outlasted stage 2; the other way round it
    // finished during the walk and cost nothing at all.
    if (waitTimer.at() > 1.0) {
      printf(" (%s of it after stage 2 finished)", fmtDuration(waitTimer.at()).c_str());
    }
    printf("\n");

    printf("\n");
    if (r.foundFactor) {
      // Stage 1 had it all along. Whatever stage 2 computed meanwhile is
      // discarded -- that is the losing side of the bet, and it still ends
      // the job with a factor.
      if (ranStage2) {
        printf("  stage 1's gcd found a factor; the overlapped stage 2 is discarded\n");
      }
      reportFactors(cfg, r.factors, r.b1Used, r.b1Used, "P-1", 0);
      log("  appended to %s\n", cfg.resultsFile.c_str());
      return 0;
    }
    log("  M%u: no factor found with B1 = %llu\n",
        cfg.exponent, (unsigned long long) r.b1Used);
    if (!ranStage2) {
      writeResultJson(cfg, "P-1", r.b1Used, r.b1Used, {}, 0);
      log("  appended to %s\n", cfg.resultsFile.c_str());
    }
  }

  if (ranStage2) {
    printf("\n");
    if (s2.interrupted) {
      log("  interrupted during stage 2; resume by running again.\n");
      return 1;
    }
    if (s2.foundFactor) {
      reportFactors(cfg, s2.factors, b1, bounds.b2, "P-1", 0);
    } else {
      log("  M%u: no factor found in stage 2 either (B1=%llu, B2=%llu)\n",
          cfg.exponent, (unsigned long long) b1, (unsigned long long) bounds.b2);
      writeResultJson(cfg, "P-1", b1, bounds.b2, {}, 0);
    }
    log("  appended to %s\n", cfg.resultsFile.c_str());
  }
  return 0;
}

static int runMain(int argc, char** argv) {
  // Unbuffered throughout. Redirected to a file, stdout is fully buffered, so
  // anything after the last explicit fflush stays invisible -- which makes a
  // stall look like it happened at the last line that happened to be flushed.
  // The in-place progress line needs unbuffered output anyway.
  setvbuf(stdout, nullptr, _IONBF, 0);

  // log() (used throughout tune.cpp, and for job outcomes below) has always
  // written to stdout only: initLog() -- the call that opens its file half --
  // was declared and defined but never actually invoked anywhere in this
  // program, upstream included. Appended, like every other file this program
  // accumulates (tune.txt, Mp_p-1_gpu-tune-config.txt, ...), so nothing here
  // needs rotation or size management; delete it to reset.
  initLog("Mp_p-1_gpu.log");

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

    // factoredTo's auto-default (0 == pick from exponent) can no longer be
    // resolved inside loadConfig -- it needs an exponent, which now comes
    // from worktodo.txt, not config.txt. Captured once here: an explicit
    // factored_to = N in config.txt applies to every queue entry; "auto" (0)
    // is re-resolved per entry, below and in the queue loop, so a later
    // exponent does not inherit an earlier one's default.
    const u32 configuredFactoredTo = cfg.factoredTo;

    // Peek at the next queued exponent WITHOUT consuming it, purely so
    // --bounds and --tune (which read cfg.exponent same as always) have
    // something to scope themselves to. The actual job loop below re-reads
    // worktodo.txt itself and consumes entries as their jobs complete.
    if (haveConfig) {
      std::vector<WorktodoEntry> peek;
      std::string werr;
      if (loadWorktodo(cfg.worktodoFile, peek, werr) && !peek.empty()) {
        cfg.exponent = peek.front().exponent;
        cfg.factoredTo = configuredFactoredTo ? configuredFactoredTo
                                              : defaultFactoredTo(cfg.exponent);
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
    if (doSelftest && selftest == "worktodo") { return runWorktodoTests(); }

    // --bounds: show the surface around the automatic choice, for the job in
    // config.txt. No GPU, so the pairing shape comes from the configured values
    // or a mid-range default rather than from a device query.
    if (doBounds) {
      if (!haveConfig) { printf("--bounds needs a valid config.txt\n"); return 2; }
      if (!cfg.exponent) { printf("--bounds needs an exponent in worktodo.txt\n"); return 2; }
      CostModel cm;
      cm.gcdIters = gcdIterCost(cfg.exponent);
      cm.tolerance = cfg.boundsTolerance;
      // The SAME shape the real run would pick, queried from the device, so the
      // bounds shown here are the bounds you will actually get. The one
      // difference: with no GPU up, the transform size comes from the smallest
      // shape that fits rather than from tune.txt's choice.
      const FFTConfig f = smallestFittingFFT(cfg.exponent);
      const u64 residueBytes = u64(f.size()) * sizeof(Word);
      const Stage2Shape s = pickStage2Shape(cfg, deviceOverride, residueBytes, cfg.b1 ? cfg.b1 : 100000);
      cm.mulsPerPrime = buildStage2Plan(1000000, 4000000, s.d, s.w).mulsPerPrime();
      printf("stage-2 shape D=%u w=%u (%u T-buffers of %.1f MB = %.2f GB)\n",
             s.d, s.w, stage2NumJ(s.d, s.w), double(residueBytes) / (1 << 20),
             double(stage2NumJ(s.d, s.w)) * double(residueBytes) / (1u << 30));
      printBoundsSurface(cfg.exponent, cfg.factoredTo, cfg.bias, cm);
      if (cfg.doPP1) {
        // Not a full surface -- just enough to make P+1's own model visible
        // through this diagnostic rather than only at job-run time. B1 and
        // B2 are both genuinely chosen for P+1 here; the pairing SHAPE (D,
        // w, T-table sizing) is still shared with P-1, a GPU-memory budget
        // decision rather than a cost-model gap.
        const Bounds pp1b = choosePP1Bounds(cfg.exponent, cfg.factoredTo, cfg.bias, cm,
                                            u32(cfg.pp1Seeds.size()), cfg.b1, cfg.b2);
        if (pp1b.b2 > pp1b.b1) {
          printf("\nP+1 (%u seeds, own B1/B2 model):\n"
                 "  B1=%llu B2=%llu  P(factor)=%.3f%% (%.3f%% + %.3f%%)  work %.1f%% of a PRP test\n",
                 u32(cfg.pp1Seeds.size()), (unsigned long long) pp1b.b1,
                 (unsigned long long) pp1b.b2, pp1b.prob() * 100,
                 pp1b.probStage1 * 100, pp1b.probStage2 * 100,
                 (pp1b.workStage1 + pp1b.workStage2) / cfg.exponent * 100);
        } else {
          printf("\nP+1 (%u seeds, own B1/B2 model):\n"
                 "  B1=%llu (no stage 2)  P(factor)=%.3f%%  work %.1f%% of a PRP test\n",
                 u32(cfg.pp1Seeds.size()), (unsigned long long) pp1b.b1, pp1b.prob() * 100,
                 pp1b.workStage1 / cfg.exponent * 100);
        }
      }
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
               cfg.exponent, cfg.worktodoFile.c_str());
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
      else { printf("\n  (no exponent in worktodo.txt, so tune.txt was not verified)\n"); }
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
      if (selftest == "all")                          { rc |= runWorktodoTests(); }
      if (selftest == "all" || selftest == "engine")  { rc |= runEngineTests(shared, &queue, args.fftSpec, selftest != "all"); }
      if (selftest == "all" || selftest == "pm1")     { rc |= runPM1Tests(shared, &queue, args.fftSpec); }
      if (selftest == "all" || selftest == "extend")  { rc |= runExtendTests(shared, &queue, args.fftSpec); }
      if (selftest == "all" || selftest == "pp1")     { rc |= runPP1Tests(shared, &queue, args.fftSpec); }
      if (selftest == "all" || selftest == "stage2")  { rc |= runStage2Tests(shared, &queue, args.fftSpec); }
      if (selftest == "all" || selftest == "b2extend") { rc |= runB2ExtendTests(shared, &queue, args.fftSpec); }
      if (selftest == "all" || selftest == "pp1stage2") { rc |= runPp1Stage2Tests(shared, &queue, args.fftSpec); }
      if (selftest != "all" && selftest != "engine" && selftest != "pm1"
          && selftest != "extend" && selftest != "pp1" && selftest != "stage2"
          && selftest != "b2extend" && selftest != "pp1stage2") {
        printf("unknown selftest '%s'\n"
               "  no GPU: gcd, exponent, stage2plan, bounds, worktodo\n"
               "  GPU:    engine, pm1, extend, pp1, stage2, b2extend, pp1stage2\n"
               "  omit =which to run them all\n",
               selftest.c_str());
        return 2;
      }
      return rc;
    }

    // ---- process worktodo.txt as a queue -----------------------------------
    // Each entry runs under this SAME cfg (method/b1/b2/bias/... are
    // config.txt-wide); only cfg.exponent (and cfg.factoredTo's auto default)
    // varies per entry. A job that throws is NOT caught here -- it propagates
    // to this function's own try/catch below, which stops the whole queue
    // (the entry stays queued, same as a single-exponent run failing today).
    for (;;) {
      std::vector<WorktodoEntry> entries;
      std::string werr;
      if (!loadWorktodo(cfg.worktodoFile, entries, werr)) {
        printf("worktodo error: %s\n", werr.c_str());
        return 2;
      }
      if (entries.empty()) {
        printf("worktodo.txt is empty -- nothing to do. "
               "Add an exponent, one per line, and run again.\n");
        return 0;
      }

      const WorktodoEntry job = entries.front();
      cfg.exponent = job.exponent;
      cfg.factoredTo = configuredFactoredTo ? configuredFactoredTo
                                            : defaultFactoredTo(cfg.exponent);
      printf("\n=== M%u (%zu queued) ===\n\n", job.exponent, entries.size());

      const int rc = runOneJob(cfg, shared, queue, fftSpec, deviceOverride);
      if (rc == 1) { return 1; }   // interrupted mid-job; entry stays queued

      std::string cerr;
      if (!consumeWorktodoEntry(cfg.worktodoFile, job, cerr)) {
        log("WARNING: could not remove completed exponent %u from worktodo.txt: %s\n",
            job.exponent, cerr.c_str());
        return 2;   // do not risk silently reprocessing forever
      }
    }
  } catch (const char* s) {
    log("\nFAILED: %s\n", s);
    return 2;
  } catch (const std::exception& e) {
    log("\nFAILED: %s\n", e.what());
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
