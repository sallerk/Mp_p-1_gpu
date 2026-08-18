// Copyright (C) Mp_p-1_gpu
//
// Gates G1 (engine) and G3 (P-1), plus the FFT benchmark. G2 lives in
// testBigInt.cpp because it needs no GPU.

#include "BigInt.h"
#include "Bounds.h"
#include "Config.h"
#include "FFTConfig.h"
#include "Gcd.h"
#include "Gpu.h"
#include "GpuCommon.h"
#include "log.h"
#include "PM1.h"
#include "Queue.h"
#include "Stage2Plan.h"
#include "Stage2Save.h"
#include "Pp1Stage2Save.h"
#include "TuneEntry.h"
#include "clwrap.h"
#include "Context.h"
#include "Args.h"
#include "common.h"
#include "timeutil.h"

#include <cstdio>
#include <string>
#include <algorithm>
#include <filesystem>
#include <map>
#include <random>
#include <vector>

using namespace std;

namespace {

int failures = 0, checks = 0;

void check(bool ok, const string& what) {
  ++checks;
  if (!ok) { ++failures; printf("   FAIL  %s\n", what.c_str()); }
}

// Pick a size-appropriate FFT for a test exponent.
//
// FFTConfig::bestFit() consults tune.txt and returns the cheapest entry whose
// maxExp merely EXCEEDS E. Once tune.txt is tuned for ~8e7, that hands a ~1e6
// test a 2-million-word transform -- 0.4 bits/word, which the engine rejects
// outright. Tests want the smallest transform that fits, so choose directly.
FFTConfig smallestFFT(u32 E, const string& spec) {
  if (!spec.empty()) { return FFTConfig{spec}; }
  for (const FFTShape& shape : FFTShape::allShapes()) {
    FFTConfig f{shape, 101, CARRY_AUTO};
    if (f.maxExp() >= E) { return f; }
  }
  throw "no FFT fits this exponent";
}

bool isExactlyNine(const Words& w) {
  if (w.empty() || w[0] != 9) { return false; }
  for (u32 i = 1; i < w.size(); ++i) { if (w[i]) { return false; } }
  return true;
}

} // namespace

// ---------------------------------------------------------------------------
// G1 -- the FFT/squaring engine.
//
// A. exactness at scale: 3^(2^n) mod (2^p-1) for n = 1000 against exact
//    integer arithmetic on the CPU (tools/prpref.py).
// B. end to end: the full p squarings for a known Mersenne PRIME must land on
//    exactly 9.
// ---------------------------------------------------------------------------
int runEngineTests(GpuCommon shared, Queue* q, const string& fftSpec, bool quick) {
  printf("G1: FFT engine\n\n");
  const int before = failures;

  struct Partial { u32 exponent, iters; u64 want; };
  static const Partial PARTIAL[] = {
    { 859433, 1000, 0x5aef4a8ea6c14b3eull},
    {1257787, 1000, 0x91d0e6e562cb2541ull},
  };

  printf("  A. exact residue vs CPU reference\n");
  for (const Partial& c : PARTIAL) {
    FFTConfig fft = smallestFFT(c.exponent, fftSpec);
    auto gpu = Gpu::make(q, c.exponent, shared, fft, {}, false);
    Timer t;
    Words r = gpu->expExp2(makeWords(c.exponent, 3), c.iters);
    const u64 got = res64(r);
    check(got == c.want, "engine residue at E=" + to_string(c.exponent));
    printf("     %s  E=%-8u n=%-5u res64=%016llx  fft=%-18s %6.2fs\n",
           got == c.want ? "PASS" : "FAIL", c.exponent, c.iters,
           (unsigned long long) got, fft.spec().c_str(), t.at());
  }

  if (quick) {
    printf("\n  (skipping the full-PRP check; run --selftest for it)\n");
  } else {
    printf("\n  B. full PRP of a known Mersenne prime, residue must be exactly 9\n");
    for (u32 E : {u32(859433), u32(1257787)}) {
      FFTConfig fft = smallestFFT(E, fftSpec);
      auto gpu = Gpu::make(q, E, shared, fft, {}, false);
      Timer t;
      Words r = gpu->expExp2(makeWords(E, 3), E);
      const bool ok = isExactlyNine(r);
      check(ok, "full PRP == 9 at E=" + to_string(E));
      printf("     %s  E=%-8u res64=%016llx  %6.1fs (%.1f us/it)\n",
             ok ? "PASS" : "FAIL", E, (unsigned long long) res64(r),
             t.at(), t.at() * 1e6 / E);
    }
  }

  printf("\nG1: %d checks, %d failed.\n\n", checks, failures - before);
  return failures == before ? 0 : 1;
}

// ---------------------------------------------------------------------------
// G3 -- P-1 stage 1 against factors that are known to be findable.
//
// From gpuowl's test-pm1/pm1.txt, whose rows are exponent,factor,bits,b1,b2.
// A row is found by stage 1 alone when B1 >= its b2 column, so these are cheap:
// B1 in the hundreds means only a few thousand squarings, yet the FFT is
// full target size.
// ---------------------------------------------------------------------------
int runPM1Tests(GpuCommon shared, Queue* q, const string& fftSpec) {
  printf("G3: P-1 stage 1 against known factors\n\n");
  const int before = failures;

  struct Vector { u32 exponent; u64 b1; const char* factor; };
  static const Vector VECTORS[] = {
    { 86599237,  659, "801888655371890025340351"   },
    { 86275669, 1297, "36963258258267123886417"    },
    { 86121953, 1873, "56111791165910924949551"    },
    { 86237119, 2749, "7349026679767614782449"     },
  };

  // The factor is known here, so there is no need to run the full GCD: q
  // divides gcd(x-1, M_p) exactly when q divides x-1, and reducing an 86-million
  // bit residue modulo an 80-bit q is a 2-limb division that takes well under a
  // second. The full GCD is exercised once at the end to prove the extraction
  // path, rather than four times at ~14 minutes each.
  printf("  A. residue check against the known factor (no GCD)\n");
  Nat firstResidue;
  u32 firstExponent = 0;

  for (const Vector& v : VECTORS) {
    Config cfg;
    cfg.exponent = v.exponent;
    cfg.fftSpec = fftSpec;
    cfg.reportEvery = 0;              // no progress line inside a test

    FFTConfig fft = smallestFFT(v.exponent, fftSpec);
    auto gpu = Gpu::make(q, v.exponent, shared, fft, {}, false);

    Timer t;
    PM1Result r = runPM1Stage1(*gpu, cfg, v.b1, false, /*doGcd=*/false);

    Nat want;
    const bool parsed = fromDecimal(v.factor, want);
    check(parsed, "parse expected factor");

    const bool ok = parsed && !want.isZero() && !r.xMinusOne.isZero() &&
                    mod(r.xMinusOne, want).isZero();
    check(ok, "P-1 residue is divisible by the known factor of M" + to_string(v.exponent));
    printf("     %s  M%-9u B1=%-5llu  q | (x-1): %-3s  (%s)\n",
           ok ? "PASS" : "FAIL", v.exponent, (unsigned long long) v.b1,
           ok ? "yes" : "NO", fmtDuration(t.at()).c_str());

    if (!firstExponent) { firstResidue = r.xMinusOne; firstExponent = v.exponent; }
  }

  // One full GCD, to prove the extraction path end to end.
  if (firstExponent && !firstResidue.isZero()) {
    printf("\n  B. full gcd(x-1, M_p) on the first vector\n");
    Nat want;
    fromDecimal(VECTORS[0].factor, want);
    Timer t;
    Nat g = gcd(firstResidue, mersenne(firstExponent));
    const bool ok = !g.isOne() && mod(g, want).isZero();
    check(ok, "full GCD extracts the known factor");
    printf("     %s  M%-9u factor %s  (%s)\n", ok ? "PASS" : "FAIL",
           firstExponent, g.hex().c_str(), fmtDuration(t.at()).c_str());
  }

  printf("\n  C. interrupt and resume reproduce an uninterrupted ladder exactly\n");
  {
    // Real bug this guards: the interrupt path used to have nothing valid to
    // checkpoint -- Ctrl-C between two scheduled saves silently discarded the
    // squarings since the last one. reportEvery=1 makes the stop point exact
    // (not rounded to a report boundary); saveEvery is set far beyond the
    // whole walk so the ONLY save that fires is the interrupt-triggered one,
    // proving that path specifically rather than a coincidental scheduled one.
    const u32 p = 86599237;
    const u64 b1 = 2000;
    FFTConfig fft = smallestFFT(p, fftSpec);
    auto gpu = Gpu::make(q, p, shared, fft, {}, false);
    const Nat E = stage1Exponent(b1, p);
    auto always = [](u64, u64) { return true; };

    const Words full = gpu->powBase3(E.toVector(), 0, always);

    const u64 total = E.bits() ? E.bits() - 1 : 0;
    const u64 stopAfter = std::max<u64>(1, total / 2);
    u64 seen = 0;
    Words savedResidue;
    u64 savedBit = ~0ull;
    bool saveCalled = false;
    auto stopHalfway = [&](u64 done, u64) { seen = done; return done < stopAfter; };
    auto captureSave = [&](const Words& w, u64 bit) {
      savedResidue = w; savedBit = bit; saveCalled = true;
    };
    gpu->powBase3(E.toVector(), 1, stopHalfway, nullptr, 0, 1'000'000'000u, captureSave);

    const bool stoppedEarly = saveCalled && seen >= stopAfter;
    check(stoppedEarly, "the ladder actually stopped before completion, and checkpointed");

    const Words resumed = gpu->powBase3(E.toVector(), 0, always, &savedResidue, savedBit);
    const bool ok = stoppedEarly && (resumed == full);
    check(ok, "resumed ladder matches an uninterrupted run");
    printf("     %s  M%u B1=%llu  stopped at %llu/%llu, resumed  full res64=%016llx  resumed res64=%016llx\n",
           ok ? "PASS" : "FAIL", p, (unsigned long long) b1,
           (unsigned long long) seen, (unsigned long long) total,
           (unsigned long long) res64(full), (unsigned long long) res64(resumed));
  }

  printf("\nG3: %d failed.\n\n", failures - before);
  return failures == before ? 0 : 1;
}

// ---------------------------------------------------------------------------
// M4 -- P+1 stage 1 (Lucas ladder).
//
// Checked against exact CPU arithmetic (tools/pp1ref.py), which itself is
// cross-checked against the plain V_{n+1} = P*V_n - V_{n-1} recurrence. The
// GPU ladder is only trustworthy if it reproduces those residues exactly.
// ---------------------------------------------------------------------------
int runPP1Tests(GpuCommon shared, Queue* q, const string& fftSpec) {
  printf("M4: P+1 Lucas ladder\n\n");
  const int before = failures;

  struct Case { u32 exponent; u64 b1; u32 seed; u64 want; };
  static const Case CASES[] = {
    {859433, 100, 3, 0x5a38d19823252e67ull},
    {859433, 500, 5, 0x9524944e21ca339eull},
  };

  printf("  V_E(seed,1) mod M_p vs exact CPU reference\n");
  for (const Case& c : CASES) {
    FFTConfig fft = smallestFFT(c.exponent, fftSpec);
    auto gpu = Gpu::make(q, c.exponent, shared, fft, {}, false);
    const Nat E = stage1Exponent(c.b1, c.exponent);
    auto quiet = [](u64, u64) { return true; };

    Timer t;
    Words v = gpu->lucasV(c.seed, E.toVector(), 0, quiet);
    const u64 got = res64(v);
    const bool ok = got == c.want;
    check(ok, "P+1 V_E at E=" + to_string(c.exponent) + " B1=" + to_string(c.b1));
    printf("     %s  M%-8u B1=%-5llu seed=%u  res64=%016llx  (%zu bits, %.2fs)\n",
           ok ? "PASS" : "FAIL", c.exponent, (unsigned long long) c.b1, c.seed,
           (unsigned long long) got, E.bits(), t.at());
    if (!ok) { printf("        expected %016llx\n", (unsigned long long) c.want); }
  }

  printf("\n  interrupt and resume reproduce an uninterrupted ladder exactly\n");
  {
    // Both A and B in lucasV are always fully carried (see Gpu.cpp), so unlike
    // powBase3 this path needed no extra carry work -- just actually calling
    // save() on interrupt, which it did not before. reportEvery=1 makes the
    // stop point exact; saveEvery is set far beyond the walk so the only save
    // that fires is the interrupt-triggered one.
    const u32 p = 859433;
    const u64 b1 = 500;
    const u32 seed = 3;
    FFTConfig fft = smallestFFT(p, fftSpec);
    auto gpu = Gpu::make(q, p, shared, fft, {}, false);
    const Nat E = stage1Exponent(b1, p);
    auto always = [](u64, u64) { return true; };

    const Words full = gpu->lucasV(seed, E.toVector(), 0, always);

    const u64 total = E.bits() ? E.bits() - 1 : 0;
    const u64 stopAfter = std::max<u64>(1, total / 2);
    u64 seen = 0;
    Words savedA, savedB;
    u64 savedBit = ~0ull;
    bool saveCalled = false;
    auto stopHalfway = [&](u64 done, u64) { seen = done; return done < stopAfter; };
    auto captureSave = [&](const Words& a, const Words& b, u64 bit) {
      savedA = a; savedB = b; savedBit = bit; saveCalled = true;
    };
    gpu->lucasV(seed, E.toVector(), 1, stopHalfway, nullptr, nullptr, 0,
               1'000'000'000u, captureSave);

    const bool stoppedEarly = saveCalled && seen >= stopAfter;
    check(stoppedEarly, "the ladder actually stopped before completion, and checkpointed");

    const Words resumed = gpu->lucasV(seed, E.toVector(), 0, always, &savedA, &savedB, savedBit);
    const bool ok = stoppedEarly && (resumed == full);
    check(ok, "resumed ladder matches an uninterrupted run");
    printf("     %s  M%u B1=%llu seed=%u  stopped at %llu/%llu  full res64=%016llx  resumed res64=%016llx\n",
           ok ? "PASS" : "FAIL", p, (unsigned long long) b1, seed,
           (unsigned long long) seen, (unsigned long long) total,
           (unsigned long long) res64(full), (unsigned long long) res64(resumed));
  }

  printf("\nM4: %d failed.\n\n", failures - before);
  return failures == before ? 0 : 1;
}

// ---------------------------------------------------------------------------
// M5b -- B1 extension.
//
// The whole claim is that extending a completed B1 gives BIT-IDENTICAL results
// to running the larger B1 from scratch. Anything less and an extended run is
// quietly not the computation it says it is, so both halves are checked:
//   A. on the CPU, E(from) * R == E(to) exactly;
//   B. on the GPU, x_from ^ R == x_to, against a from-scratch run.
// ---------------------------------------------------------------------------
int runExtendTests(GpuCommon shared, Queue* q, const string& fftSpec) {
  printf("M5b: B1 extension\n\n");
  const int before = failures;

  printf("  A. E(from) * R == E(to)   (CPU)\n");
  const u32 p = 86599237;
  for (auto pr : {pair<u64,u64>{659, 1873}, {1000, 5000}, {5000, 50000},
                  {50000, 200000}}) {
    const Nat eFrom = stage1Exponent(pr.first, p);
    const Nat eTo = stage1Exponent(pr.second, p);
    const Nat r = stage1ExponentDelta(pr.first, pr.second, p);
    const bool ok = mul(eFrom, r) == eTo;
    check(ok, "delta " + to_string(pr.first) + " -> " + to_string(pr.second));
    printf("     %s  %6llu -> %-7llu  E %zu + R %zu = %zu bits\n",
           ok ? "PASS" : "FAIL",
           (unsigned long long) pr.first, (unsigned long long) pr.second,
           eFrom.bits(), r.bits(), eTo.bits());
  }

  printf("\n  B. x_from ^ R == x_to     (GPU, vs from scratch)\n");
  {
    const u64 b1From = 659, b1To = 1873;
    FFTConfig fft = smallestFFT(p, fftSpec);
    auto gpu = Gpu::make(q, p, shared, fft, {}, false);
    auto quiet = [](u64, u64) { return true; };

    const Nat eFrom = stage1Exponent(b1From, p);
    const Nat eTo = stage1Exponent(b1To, p);
    const Nat r = stage1ExponentDelta(b1From, b1To, p);

    Timer t1;
    Words xFrom = gpu->powBase3(eFrom.toVector(), 0, quiet);
    Words xTo = gpu->powBase3(eTo.toVector(), 0, quiet);
    const double freshSecs = t1.at();

    Timer t2;
    Words extended = gpu->powResidue(xFrom, r.toVector(), 0, quiet);
    const double extSecs = t2.at();

    const bool ok = extended == xTo;
    check(ok, "extended residue == from-scratch residue");
    printf("     %s  B1 %llu -> %llu   res64 extended=%016llx scratch=%016llx\n",
           ok ? "PASS" : "FAIL",
           (unsigned long long) b1From, (unsigned long long) b1To,
           (unsigned long long) res64(extended), (unsigned long long) res64(xTo));
    printf("     extension %zu squarings vs %zu from scratch (%.0f%% saved), %.1fs vs %.1fs\n",
           r.bits() - 1, eTo.bits() - 1,
           100.0 * (1.0 - double(r.bits()) / double(eTo.bits())),
           extSecs, freshSecs);
  }

  printf("\nM5b: %d failed.\n\n", failures - before);
  return failures == before ? 0 : 1;
}

// ---------------------------------------------------------------------------
// Verified FFT selection.
//
// Tuning this GPU showed 5 of ~20 configs computing WRONG results -- all of
// them width=512/height=512 NTT shapes. Tune catches that with a Gerbicz check
// and scores them unusable, so tune.txt is clean. But FFTConfig::bestFit's
// no-tune fallback ("first shape that fits") verifies nothing, so a run without
// tune.txt can silently commit hours to a broken transform.
//
// So: before a long run, actually check the config. Gpu::timePRP() runs a few
// hundred iterations with the same check tune uses, and returns a sentinel cost
// on failure. A couple of seconds against a multi-hour stage 1.
// ---------------------------------------------------------------------------
namespace {

// timePRP marks a failed correctness check by returning 0.1 s/iteration.
const double TIME_PRP_FAILED = 99999.0;

bool verifyOne(GpuCommon shared, Queue* q, u32 E, const FFTConfig& fft, double* cost) {
  try {
    const double c = Gpu::make(q, E, shared, fft, {}, false)->timePRP(10);
    if (cost) { *cost = c; }
    return c < TIME_PRP_FAILED;
  } catch (const char* s) {
    printf("      %s: threw \"%s\"\n", fft.spec().c_str(), s);
    return false;
  } catch (const std::exception& e) {
    printf("      %s: threw \"%s\"\n", fft.spec().c_str(), e.what());
    return false;
  }
}

} // namespace

// After --tune: verify what it wrote and recommend a line for config.txt.
//
// Recommending tune's winner UNVERIFIED would be actively dangerous. A real
// tune on this machine produced a tune.txt whose every entry computes wrong
// results: its kernel-option search happened to pick a broken shape as the
// baseline, so all 50 option measurements returned the failure sentinel, "best"
// was chosen among identical failures, and those meaningless options then
// corrupted the FFT timings that followed. Auto-pinning that winner would have
// committed hours to garbage.
void recommendFFT(GpuCommon shared, Queue* q, u32 E) {
  vector<TuneEntry> entries = TuneEntry::readTuneFile(*shared.args);
  printf("\n=======================================================================\n");
  printf(" Verifying what --tune wrote, for M%u\n", E);
  printf("=======================================================================\n\n");

  if (entries.empty()) {
    printf("  tune.txt is empty -- nothing was timed.\n");
    return;
  }

  string bestSpec;
  double bestCost = 0;
  int usable = 0, broken = 0;

  for (const TuneEntry& e : entries) {
    // maxExp() alone is only an upper bound -- a transform can equally be too
    // LARGE for this exponent (too few bits/word), which throws rather than
    // mis-computing. Skip those here too (chooseVerifiedFFT's own `usable`
    // applies the same minBpw check) so this report only verifies genuinely
    // plausible candidates instead of ones known in advance to be unusable.
    if (E > e.fft.maxExp()) { continue; }        // cannot hold this exponent
    if (double(E) / double(e.fft.size()) < e.fft.minBpw()) { continue; }  // too few bits/word
    printf("  %-22s %8.1f us/it ... ", e.fft.spec().c_str(), e.cost);
    fflush(stdout);
    if (verifyOne(shared, q, E, e.fft, nullptr)) {
      ++usable;
      printf("verified\n");
      if (bestSpec.empty()) { bestSpec = e.fft.spec(); bestCost = e.cost; }
    } else {
      ++broken;
      printf("COMPUTES WRONG RESULTS -- do not use\n");
    }
  }

  printf("\n");
  if (!bestSpec.empty()) {
    printf("  Verified: fft = %s (%.1f us/it, correct at M%u).\n",
           bestSpec.c_str(), bestCost, E);
    printf("  fft = auto (config.txt's default) already picks and caches this same\n"
           "  entry automatically -- pin it explicitly only to skip that one-time\n"
           "  re-verify:\n\n      fft = %s\n\n", bestSpec.c_str());
  }
  if (broken) {
    printf("\n  Note: %d of %d entries large enough for M%u compute WRONG results on\n"
           "  this hardware. That's usually just a config OpenCL mis-compiles here --\n"
           "  no action needed: this check and the automatic fft = auto path both skip\n"
           "  broken entries and use the next working one, same as just happened above.\n"
           "  To stop tune.txt listing them, delete tune.txt (not\n"
           "  Mp_p-1_gpu-tune-config.txt, which holds unrelated kernel options) and\n"
           "  re-run --tune.\n",
           broken, usable + broken, E);
  }
  if (!usable) {
    printf("\n  NOTHING in tune.txt is usable. Do not set fft= from it.\n"
           "  Find a working config with:  Mp_p-1_gpu.exe --bench -E %u\n", E);
  }
}

// Cache of verification verdicts, so `fft = auto` costs nothing after the first
// run. Keyed by exponent and spec, and scoped to the GPU + driver: a driver
// update can change kernel codegen, so a verdict from a different driver is not
// evidence about this one.
namespace {

const char* VERIFY_CACHE = "fft-verified.txt";

string deviceTag(Queue* q) {
  const cl_device_id id = q->context->deviceId();
  return getDeviceName(id) + " / " + getDriverVersion(id);
}

// spec -> verified?   (only entries for this exponent and this device)
std::map<string, bool> loadVerifyCache(u32 E, const string& tag) {
  std::map<string, bool> out;
  FILE* f = fopen(VERIFY_CACHE, "r");
  if (!f) { return out; }
  char line[512];
  bool tagOk = false;
  while (fgets(line, sizeof(line), f)) {
    string s = line;
    while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) { s.pop_back(); }
    if (s.rfind("# device ", 0) == 0) { tagOk = (s.substr(9) == tag); continue; }
    if (s.empty() || s[0] == '#' || !tagOk) { continue; }
    char spec[128], verdict[16];
    unsigned long long exp = 0;
    if (sscanf(s.c_str(), "%llu %127s %15s", &exp, spec, verdict) == 3 && exp == E) {
      out[spec] = (string(verdict) == "ok");
    }
  }
  fclose(f);
  return out;
}

void appendVerifyCache(u32 E, const string& tag, const string& spec, bool ok,
                       double cost) {
  const bool fresh = !std::filesystem::exists(VERIFY_CACHE);
  FILE* f = fopen(VERIFY_CACHE, "a");
  if (!f) { return; }
  if (fresh) {
    fprintf(f, "# Verification verdicts for FFT configs. Delete to re-verify.\n");
  }
  static bool wroteTag = false;
  if (!wroteTag) { fprintf(f, "# device %s\n", tag.c_str()); wroteTag = true; }
  fprintf(f, "%u %s %s %.1f\n", E, spec.c_str(), ok ? "ok" : "bad", cost);
  fclose(f);
}

} // namespace

FFTConfig chooseVerifiedFFT(GpuCommon shared, Queue* q, u32 E,
                            const string& forcedSpec, bool verify) {
  // Candidate order: an explicit -fft wins; otherwise tune.txt by measured
  // cost, then bestFit's own fallback as a last resort.
  vector<FFTConfig> candidates;
  auto addUnique = [&](const FFTConfig& f) {
    for (const FFTConfig& c : candidates) { if (c.spec() == f.spec()) { return; } }
    candidates.push_back(f);
  };

  // maxExp() alone is only an UPPER bound. A transform can equally be too LARGE
  // for an exponent: below fft.minBpw() bits per word the engine refuses to run
  // at all. With a tune.txt built for ~8e7, every entry passes maxExp for an
  // 8e5 exponent and the cheapest is a 2-million-word transform at 0.41
  // bits/word, so the job died on the first config it tried. Apply the engine's
  // own test here instead of offering it configs it will reject.
  auto usable = [E](const FFTConfig& f) {
    return E <= f.maxExp() && double(E) / double(f.size()) >= f.minBpw();
  };

  if (!forcedSpec.empty()) {
    candidates.push_back(FFTConfig{forcedSpec});
  } else {
    // 1. tune.txt, cheapest first.
    for (const TuneEntry& e : TuneEntry::readTuneFile(*shared.args)) {
      if (usable(e.fft)) { addUnique(e.fft); }
    }
    const size_t nTuned = candidates.size();
    if (!nTuned) {
      log("  no tune.txt entry covers M%u.\n"
          "  Consider:  Mp_p-1_gpu.exe --tune quick=10,minexp=%u,maxexp=%u\n",
          E, E, E);
    }

    // 2. An INDEPENDENT fallback. Note FFTConfig::bestFit() also reads
    //    tune.txt, so calling it here would just re-offer the same entries --
    //    useless when every tuned entry is broken, which is exactly the case
    //    that needs a fallback. Enumerate the shapes directly instead.
    //
    //    Observed on this GPU: every width=512/height=512 NTT shape fails its
    //    correctness check deterministically (same wrong residue every time),
    //    and tune.txt was populated entirely with those. Everything that has
    //    ever worked here has height 256, so ordering by size and verifying is
    //    what actually finds a usable config.
    for (u32 variant : {101u, 202u}) {
      for (const FFTShape& shape : FFTShape::allShapes()) {
        FFTConfig f{shape, variant, CARRY_AUTO};
        if (usable(f)) { addUnique(f); }
      }
    }
    if (nTuned) {
      printf("  %u tuned candidate(s), then %u untuned as fallback\n",
             u32(nTuned), u32(candidates.size() - nTuned));
    }
  }

  if (!verify) { return candidates.front(); }

  // Verifying costs a few seconds each; do not crawl through hundreds of them.
  const size_t MAX_TRIES = 8;
  const size_t tries = std::min(candidates.size(), MAX_TRIES);

  const string tag = deviceTag(q);
  const std::map<string, bool> cache = loadVerifyCache(E, tag);

  for (size_t i = 0; i < tries; ++i) {
    const FFTConfig& fft = candidates[i];
    const string spec = fft.spec();

    // A remembered verdict, from this exponent on this GPU and driver.
    if (auto it = cache.find(spec); it != cache.end()) {
      if (it->second) {
        printf("  FFT %s (verified earlier)\n", spec.c_str());
        return fft;
      }
      printf("  skipping %s (failed verification earlier)\n", spec.c_str());
      continue;
    }

    printf("  verifying FFT %s ... ", spec.c_str());
    fflush(stdout);
    double cost = 0;
    Timer t;
    const bool ok = verifyOne(shared, q, E, fft, &cost);
    appendVerifyCache(E, tag, spec, ok, cost);
    if (ok) {
      printf("OK (%.0f us/it, %.1fs)\n", cost, t.at());
      return fft;
    }
    printf("FAILED its correctness check -- NOT USING IT\n");
    if (i + 1 < tries) { printf("  trying the next candidate\n"); }
  }

  log("\n  No FFT config passed verification for M%u.\n"
      "  Run --tune, or pin a known-good one with -fft <spec>.\n", E);
  throw "no usable FFT config";
}

// ---------------------------------------------------------------------------
// Exponent-build isolation: pure CPU, no GPU, no threads. Used to tell a slow
// product tree apart from a stall elsewhere.
// ---------------------------------------------------------------------------
int runExponentTest() {
  setvbuf(stdout, nullptr, _IONBF, 0);
  printf("stage1Exponent build timing (CPU only)\n\n");
  for (u64 b1 : {u64(10000), u64(50000), u64(200000), u64(1000000), u64(3000000)}) {
    printf("  B1 = %-9llu ... ", (unsigned long long) b1);
    Timer t;
    Nat e = stage1Exponent(b1, 82589959);
    printf("%zu bits in %.2fs\n", e.bits(), t.at());
  }
  return 0;
}

// ---------------------------------------------------------------------------
// FFT benchmark (moved out of main).
// ---------------------------------------------------------------------------
void benchmarkFFT(GpuCommon shared, Queue* q, u32 E, u32 iters, double sizeSlack,
                  const string& onlySpec) {
  setvbuf(stdout, nullptr, _IONBF, 0);
  printf("Benchmark: E=%u, %u iterations per config, size slack %.2fx\n\n",
         E, iters, sizeSlack);

  vector<FFTConfig> todo;
  u64 smallest = 0;
  if (!onlySpec.empty()) {
    todo.push_back(FFTConfig{onlySpec});
    smallest = todo.front().size();
  } else {
    for (const FFTShape& shape : FFTShape::allShapes()) {
      FFTConfig fft{shape, 101, CARRY_AUTO};
      if (fft.maxExp() < E) { continue; }
      if (!smallest) { smallest = fft.size(); }
      if (fft.size() > smallest * sizeSlack) { break; }
      todo.push_back(fft);
    }
  }
  printf("%u candidate configs (smallest viable size %llu)\n\n",
         u32(todo.size()), (unsigned long long) smallest);

  double best = 1e30;
  string bestSpec;
  for (u32 i = 0; i < todo.size(); ++i) {
    FFTConfig fft = todo[i];
    // Announced before building, so a stall names the config that caused it.
    printf("  [%2u/%2u] %-22s size %-9llu bpw %5.2f ... ",
           i + 1, u32(todo.size()), fft.spec().c_str(),
           (unsigned long long) fft.size(), double(E) / fft.size());
    try {
      Timer stage;
      printf("make ");
      auto gpu = Gpu::make(q, E, shared, fft, {}, false);
      printf("(%.1fs) warmup ", stage.reset());
      gpu->expExp2(makeWords(E, 3), 4);
      printf("(%.1fs) run ", stage.reset());
      Timer timer;
      gpu->expExp2(makeWords(E, 3), iters);
      const double us = timer.at() * 1e6 / iters;
      if (us < best) { best = us; bestSpec = fft.spec(); }
      printf("-> %8.1f us/it\n", us);
    } catch (const char* s) {
      printf("SKIP (%s)\n", s);
    } catch (const std::exception& e) {
      printf("SKIP (%s)\n", e.what());
    }
  }
  if (!bestSpec.empty()) { printf("\n  fastest: %s at %.1f us/it\n", bestSpec.c_str(), best); }
}

// ---------------------------------------------------------------------------
// M6b -- P-1 stage 2 engine.
//
// The plan itself is audited separately (--selftest=stage2plan); here the plan
// is INPUT, and what is under test is the GPU: the T_j table, the A/S
// recurrence, the un-normalised subWords difference, and the accumulator.
//
// The CPU reference is deliberately built differently from the GPU. Every A_m
// is a fresh modular exponentiation instead of being carried forward by
// A *= S, so an error in that recurrence cannot be mirrored on both sides and
// cancel out. The arithmetic paths share nothing either: Karatsuba and Knuth-D
// division on one side, an FFT with carry propagation on the other.
// ---------------------------------------------------------------------------
namespace {

// x mod (2^p - 1) by folding the high half down, which is linear. BigInt's
// generic mod() is Knuth-D at O(limbs^2), and the smallest FFT this engine
// supports needs p ~ 860000, i.e. 13400 limbs -- generic division there would
// make this reference take minutes per case instead of seconds.
Nat modMersenneFast(Nat x, u32 p) {
  const Nat M = mersenne(p);
  while (x.bits() > p) {
    Nat hi = shr(x, p);
    x = add(hi, sub(x, shl(hi, p)));
  }
  if (cmp(x, M) >= 0) { x = sub(x, M); }
  return x;
}

Nat mulModM(const Nat& a, const Nat& b, u32 p) { return modMersenneFast(mul(a, b), p); }

Nat stage2Reference(u32 p, u32 base, const Stage2Plan& plan) {
  const Nat M = mersenne(p);
  const Nat x{base};

  // T_j = x^(j^2), by the second-difference recurrence
  // x^((j+1)^2) = x^(j^2) * x^(2j+1). The GPU builds the same table by direct
  // square-and-multiply on each j^2, so the two share no method.
  const u64 jmax = u64(plan.w) * plan.d / 2;
  vector<int> idxOf(size_t(jmax), -1);
  for (size_t i = 0; i < plan.jset.size(); ++i) { idxOf[plan.jset[i]] = int(i); }

  vector<Nat> T(plan.jset.size());
  Nat u{1};                        // x^(j^2) at j = 0
  Nat step = x;                    // x^(2j+1) at j = 0
  const Nat xx = mulModM(x, x, p);
  for (u64 j = 0; j + 1 < jmax; ++j) {
    u = mulModM(u, step, p);       // u = x^((j+1)^2)
    step = mulModM(step, xx, p);   // step = x^(2(j+1)+1)
    const int i = idxOf[size_t(j + 1)];
    if (i >= 0) { T[size_t(i)] = u; }
  }

  // A_m by a fresh exponentiation per block, NOT by the A *= S recurrence the
  // GPU uses. That is the whole point of this reference.
  Nat acc{1};
  for (u64 m = plan.mFirst; m <= plan.mLast; ++m) {
    bool any = false;
    for (size_t i = 0; i < plan.jset.size() && !any; ++i) { any = plan.slot(m, i); }
    if (!any) { continue; }

    const u64 e = (m * plan.d) * (m * plan.d);
    Nat A{1};
    for (int b = 63 - std::countl_zero(e); b >= 0; --b) {
      A = mulModM(A, A, p);
      if ((e >> b) & 1) { A = mulModM(A, x, p); }
    }

    for (size_t i = 0; i < plan.jset.size(); ++i) {
      if (!plan.slot(m, i)) { continue; }
      // A + M - T is always non-negative, and the reduction handles the rest,
      // so there is no need to branch on A >= T.
      acc = mulModM(acc, sub(add(A, M), T[i]), p);
    }
  }
  return modMersenneFast(acc, p);
}

} // namespace

int runStage2Tests(GpuCommon shared, Queue* q, const string& fftSpec) {
  printf("M6b: P-1 stage 2 engine vs exact CPU arithmetic\n\n");
  const int before = failures;

  struct Case { u32 exponent; u64 b1; u64 b2; u32 d; u32 w; };
  // The exponent is bounded below by the smallest FFT this engine has -- around
  // 860000 -- so the CPU side works on 13400-limb numbers and the ranges have to
  // stay small to keep the reference to seconds rather than minutes.
  static const Case CASES[] = {
    {859433, 1000, 1800, 210, 1},   // no window: exactly one slot per prime
    {859433, 1000, 1800, 210, 3},   // window: several candidate slots per prime
    {859433, 2000, 3000, 420, 1},
  };

  printf("  accumulated product mod M_p, GPU vs CPU\n");
  for (const Case& c : CASES) {
    Timer t;
    const Stage2Plan plan = buildStage2Plan(c.b1, c.b2, c.d, c.w);

    FFTConfig fft = smallestFFT(c.exponent, fftSpec);
    auto gpu = Gpu::make(q, c.exponent, shared, fft, {}, false);

    const Words got = gpu->stage2(makeWords(c.exponent, 3), plan, 0,
                                  [](u64, u64) { return true; });
    const Nat want = stage2Reference(c.exponent, 3, plan);
    const Nat have = mod(fromWords(got), mersenne(c.exponent));

    const bool ok = (have == want) && !want.isZero();
    check(ok, "stage-2 accumulator matches the CPU reference");
    printf("     %s  M%-6u B1=%-5llu B2=%-5llu D=%-4u w=%u  %5llu slots  %s  (%s)\n",
           ok ? "PASS" : "FAIL", c.exponent, (unsigned long long) c.b1,
           (unsigned long long) c.b2, c.d, c.w,
           (unsigned long long) plan.nSlots,
           ok ? "match" : ("MISMATCH got " + have.hex() + " want " + want.hex()).c_str(),
           fmtDuration(t.at()).c_str());
  }

  // The one real risk in this design: A_m - T_j is not carry normalised, so its
  // words are about one bit wider than usual and every accumulator multiply sees
  // a wider operand than a squaring would. At a realistic exponent, where
  // bits/word is closest to the limit, that is where it would break.
  //
  // Testing this by round-off statistics turned out to be useless: the
  // transforms chosen at this size are largely NTT (exact modular arithmetic),
  // so the reported ROE describes only the float component and cannot be read as
  // headroom. The direct test instead: run the SAME plan twice, once with the
  // wide difference and once with the difference carried back into the balanced
  // range first. Both multiply the identical mathematical value, so the two
  // accumulators must agree bit for bit -- and if the extra bit ever overflowed
  // the transform, they would not.
  printf("\n  un-normalised difference at a realistic exponent\n");
  {
    const u32 p = 86599237;
    try {
      FFTConfig fft = smallestFFT(p, fftSpec);
      auto gpu = Gpu::make(q, p, shared, fft, {}, false);

      // Small: the normalised control pays a GPU->CPU->GPU round trip of ~10 MB
      // per slot.
      const Stage2Plan plan = buildStage2Plan(1000000, 1005000, 210, 3);
      auto never = [](u64, u64) { return true; };

      const Words wide = gpu->stage2(makeWords(p, 3), plan, 0, never);
      const Words norm = gpu->stage2(makeWords(p, 3), plan, 0, never, nullptr, true);

      const bool ok = (wide == norm);
      check(ok, "wide difference gives the same accumulator as a normalised one");
      printf("     %s  M%u  %s  %llu slots  %s\n", ok ? "PASS" : "FAIL", p,
             fft.spec().c_str(), (unsigned long long) plan.nSlots,
             ok ? "identical" : "DIFFER -- the extra bit is overflowing the transform");

      // Marginal cost per multiply, by differencing two plans that share the
      // same D and w. Setup -- the T_j table and the A/S/C seeds -- is identical
      // in both and cancels exactly, which a single run cannot separate out: on
      // a plan this small setup is most of the wall time, and dividing total
      // time by loop multiplies would overstate the cost roughly tenfold.
      const Stage2Plan small = buildStage2Plan(1000000, 1010000, 210, 1);
      const Stage2Plan big   = buildStage2Plan(1000000, 1100000, 210, 1);
      auto mulsOf = [](const Stage2Plan& q2) { return q2.nSlots + 2 * (q2.nBlocks() - 1); };

      Timer t;
      gpu->stage2(makeWords(p, 3), small, 0, never);
      const double tSmall = t.at();
      gpu->stage2(makeWords(p, 3), big, 0, never);
      const double tBig = t.at();

      const double usPerMul = (tBig - tSmall) * 1e6 / double(mulsOf(big) - mulsOf(small));

      // The bounds model costs stage 2 in units of stage-1 squarings, so the
      // ratio has to be measured here, on this transform, in this run -- not
      // quoted from a stage-1 timing taken elsewhere.
      const u32 nSq = 2000;
      Timer ts;
      gpu->expExp2(makeWords(p, 3), nSq);
      const double usPerSquaring = ts.at() * 1e6 / nSq;

      printf("     %.0f us per accumulator multiply (marginal, %llu vs %llu muls)\n",
             usPerMul, (unsigned long long) mulsOf(big), (unsigned long long) mulsOf(small));
      printf("     %.0f us per stage-1 squaring on the same transform\n", usPerSquaring);
      printf("     ratio %.2fx -- a multiply needs two forward transforms plus an\n"
             "     inverse, where the fused-carry squaring needs one of each\n",
             usPerMul / usPerSquaring);
    } catch (const char* s) {
      printf("     SKIP (%s)\n", s);
    } catch (const std::exception& e) {
      printf("     SKIP (%s)\n", e.what());
    }
  }

  printf("\nM6b: %d failed.\n\n", failures - before);
  return failures == before ? 0 : 1;
}

// ---------------------------------------------------------------------------
// M6e -- B2 extension.
//
// The claim: walking (b2old, b2new] with the accumulator seeded from a completed
// (b1, b2old] gives the same factoring power as walking (b1, b2new] in one go.
// It is not the same NUMBER -- the two cover the same primes through different
// slots, and each slot carries a harmless x^(j^2) unit -- so a res64 comparison
// against a from-scratch run would be wrong to demand. Three things are checked
// instead, each aimed at a different way this could be broken:
//
//   A. the seeded accumulator is EXACTLY the product of the two ranges'
//      accumulators, against the same independent CPU reference M6b uses;
//   B. a real factor of a real Mersenne number, whose missing prime lies in the
//      gap, is found by the extension -- and is NOT found before it, so the
//      test can fail;
//   C. no stale, foreign or unfinished accumulator can be picked up as a seed.
//
// A and B use DIFFERENT pairing shapes for the two halves on purpose: nothing
// requires them to agree, and this is where that is asserted.
// ---------------------------------------------------------------------------
int runB2ExtendTests(GpuCommon shared, Queue* q, const string& fftSpec) {
  printf("M6e: B2 extension\n\n");
  const int before = failures;

  auto quiet = [](u64, u64) { return true; };

  printf("  A. seeded accumulator == product of the two ranges   (GPU vs exact CPU)\n");
  {
    const u32 p = 859433;
    const u64 b1 = 2000, mid = 3000, b2 = 5000;
    Timer t;
    const Stage2Plan lo = buildStage2Plan(b1, mid, 420, 1);
    const Stage2Plan hi = buildStage2Plan(mid, b2, 210, 3);

    FFTConfig fft = smallestFFT(p, fftSpec);
    auto gpu = Gpu::make(q, p, shared, fft, {}, false);
    const Words x = makeWords(p, 3);

    const Words accLo = gpu->stage2(x, lo, 0, quiet);
    const Words accExt = gpu->stage2(x, hi, 0, quiet, nullptr, false, nullptr, 0, {},
                                     nullptr, &accLo);

    const Nat want = mulModM(stage2Reference(p, 3, lo), stage2Reference(p, 3, hi), p);
    const Nat have = mod(fromWords(accExt), mersenne(p));

    const bool ok = (have == want) && !want.isZero();
    check(ok, "extended accumulator equals acc(low) * acc(gap)");
    printf("     %s  M%u  (%llu,%llu] D=420 w=1  then  (%llu,%llu] D=210 w=3  %s  (%s)\n",
           ok ? "PASS" : "FAIL", p, (unsigned long long) b1, (unsigned long long) mid,
           (unsigned long long) mid, (unsigned long long) b2,
           ok ? "match" : ("MISMATCH got " + have.hex() + " want " + want.hex()).c_str(),
           fmtDuration(t.at()).c_str());

    // A seeded walk over an empty seed must reduce to the ordinary one, or the
    // seeding path is not the same code path production uses.
    const Words one = makeWords(p, 1);
    const Words accSeeded1 = gpu->stage2(x, lo, 0, quiet, nullptr, false, nullptr, 0, {},
                                         nullptr, &one);
    const bool ok1 = (accSeeded1 == accLo);
    check(ok1, "seeding with 1 reproduces the unseeded accumulator");
    printf("     %s  seed = 1 gives the identical accumulator\n", ok1 ? "PASS" : "FAIL");
  }

  printf("\n  B. a real factor whose missing prime is in the gap   (M86255591)\n");
  {
    // From gpuowl's test-pm1/pm1.txt: q = 2kp+1 divides M86255591, and k is
    // 1997-smooth apart from the single prime 30949. So stage 1 to B1=2000 plus
    // stage 2 to B2=20000 must NOT find it, and extending to B2=40000 must.
    const u32 p = 86255591;
    const char* factor = "28428815677762982290409";
    const u64 b1 = 2000, mid = 20000, b2 = 40000;

    Nat qFactor;
    check(fromDecimal(factor, qFactor), "parse expected factor");

    Config cfg;
    cfg.exponent = p;
    cfg.fftSpec = fftSpec;
    cfg.reportEvery = 0;
    cfg.checkpoint = false;             // hermetic: leave no save files behind

    FFTConfig fft = smallestFFT(p, fftSpec);
    auto gpu = Gpu::make(q, p, shared, fft, {}, false);

    Timer t1;
    PM1Result r = runPM1Stage1(*gpu, cfg, b1, false, /*doGcd=*/false);
    printf("     stage 1 to B1=%llu: %llu squarings (%s)\n",
           (unsigned long long) b1, (unsigned long long) r.squarings,
           fmtDuration(t1.at()).c_str());

    auto divisible = [&](const Words& acc) {
      return !qFactor.isZero() && mod(fromWords(acc), qFactor).isZero();
    };

    const Stage2Plan lo = buildStage2Plan(b1, mid, 210, 1);
    const Stage2Plan hi = buildStage2Plan(mid, b2, 420, 3);

    Timer t2;
    const Words accLo = gpu->stage2(r.residue, lo, 0, quiet);
    const bool foundEarly = divisible(accLo);
    check(!foundEarly, "the factor is NOT reachable with B2=20000");
    printf("     %s  (%llu,%llu] alone: q divides acc? %-3s (must be no, or the\n"
           "           test proves nothing)  (%s)\n", foundEarly ? "FAIL" : "PASS",
           (unsigned long long) b1, (unsigned long long) mid,
           foundEarly ? "YES" : "no", fmtDuration(t2.at()).c_str());

    Timer t3;
    const Words accExt = gpu->stage2(r.residue, hi, 0, quiet, nullptr, false, nullptr, 0, {},
                                     nullptr, &accLo);
    const bool foundExt = divisible(accExt);
    check(foundExt, "the extension to B2=40000 finds the factor");
    printf("     %s  extended to (%llu,%llu]: q divides acc? %-3s  (%s)\n",
           foundExt ? "PASS" : "FAIL", (unsigned long long) mid, (unsigned long long) b2,
           foundExt ? "yes" : "NO", fmtDuration(t3.at()).c_str());

    // Control: one walk over the whole range must find it too. If this failed
    // while the extension passed, the extension would be finding it for the
    // wrong reason.
    Timer t4;
    const Stage2Plan full = buildStage2Plan(b1, b2, 210, 1);
    const Words accFull = gpu->stage2(r.residue, full, 0, quiet);
    const bool foundFull = divisible(accFull);
    check(foundFull, "a single walk over (b1, b2] finds the same factor");
    printf("     %s  from scratch (%llu,%llu]: q divides acc? %-3s  (%s)\n",
           foundFull ? "PASS" : "FAIL", (unsigned long long) b1, (unsigned long long) b2,
           foundFull ? "yes" : "NO", fmtDuration(t4.at()).c_str());
  }

  printf("\n  C. only a matching, completed accumulator may be reused   (CPU)\n");
  {
    // Exponent 1 cannot collide with a real save file.
    const u32 p = 1;
    const u64 b1 = 100, b2 = 200;
    const u64 xr = 0x0123456789abcdefull;
    const string path = defaultStage2Path(p, b1, b2);

    Stage2State st;
    st.exponent = p;
    st.b1 = b1;
    st.b2 = b2;
    st.d = 210;
    st.w = 1;
    st.complete = true;
    st.xRes64 = xr;
    st.acc = Words{1, 2, 3, 4};

    string err;
    check(saveStage2(path, st, err), "write a completed record");

    struct Case { const char* what; u32 e; u64 b1; u64 xr; bool want; };
    static const Case CASES[] = {
      {"matching (exponent, B1, x)", 1, 100, 0x0123456789abcdefull, true },
      {"different exponent",         2, 100, 0x0123456789abcdefull, false},
      {"different B1",               1, 101, 0x0123456789abcdefull, false},
      {"different stage-1 residue",  1, 100, 0xdeadbeefdeadbeefull, false},
      {"unknown stage-1 residue",    1, 100, 0,                     false},
    };
    for (const Case& c : CASES) {
      Stage2State got;
      string e2;
      const bool got_ok = loadCompletedStage2(path, c.e, c.b1, c.xr, got, e2);
      const bool ok = (got_ok == c.want) && (!got_ok || got.acc == st.acc);
      check(ok, string("seed acceptance: ") + c.what);
      const string outcome = got_ok ? string("accepted") : "rejected (" + e2 + ")";
      printf("     %s  %-28s -> %s\n", ok ? "PASS" : "FAIL", c.what, outcome.c_str());
    }

    // An UNFINISHED walk is a resume point, never a seed: its accumulator is
    // missing every slot after the position it stopped at.
    Stage2State partial = st;
    partial.complete = false;
    partial.m = 5;
    check(saveStage2(path, partial, err), "overwrite with a partial record");
    Stage2State got;
    string e2;
    const bool tookPartial = loadCompletedStage2(path, p, b1, xr, got, e2);
    check(!tookPartial, "an unfinished accumulator is never used as a seed");
    const string partialOut = tookPartial ? string("ACCEPTED") : "rejected (" + e2 + ")";
    printf("     %s  %-28s -> %s\n", tookPartial ? "FAIL" : "PASS",
           "unfinished walk", partialOut.c_str());

    // fromB2 separates two walks that agree on every other field: one seeded
    // from a completed lower range, one that covered the whole range itself.
    Stage2State want2;
    want2.exponent = p; want2.b1 = b1; want2.b2 = b2; want2.d = 210; want2.w = 1;
    want2.xRes64 = xr; want2.fromB2 = 150;
    const bool crossed = loadStage2(path, want2, got, e2);
    check(!crossed, "a from-scratch checkpoint cannot resume an extension");
    const string crossOut = crossed ? string("ACCEPTED") : "rejected (" + e2 + ")";
    printf("     %s  %-28s -> %s\n", crossed ? "FAIL" : "PASS",
           "fromB2 mismatch", crossOut.c_str());

    // And the directory scan the driver relies on.
    const string other = defaultStage2Path(p, b1, 500);
    Stage2State st2 = st;
    st2.b2 = 500;
    saveStage2(other, st2, err);
    const vector<u64> found = findStage2Saves(p, b1);
    const bool listed = found.size() == 2 && found[0] == 500 && found[1] == 200;
    check(listed, "findStage2Saves returns every B2, largest first");
    printf("     %s  %-28s -> %zu found\n", listed ? "PASS" : "FAIL",
           "directory scan", found.size());

    std::error_code ec;
    filesystem::remove(path, ec);
    filesystem::remove(other, ec);
  }

  printf("\nM6e: %d failed.\n\n", failures - before);
  return failures == before ? 0 : 1;
}

// ---------------------------------------------------------------------------
// M7 -- P+1 stage 2.
//
// P-1's exact pairing trick (A_m=x^((mD)^2), T_j=x^(j^2), difference-of-
// squares factoring (mD-j)(mD+j)) does not transfer: P-1 works with an
// EXPLICIT ring element (x^n, freely divisible); P+1's group element is never
// explicit, only its Lucas-V trace is computable in Z/NZ. Writing
// y1 = V_E(seed,1) (stage 1's result), beta = alpha^E for a hidden factor q's
// P+1 unit alpha, and W_n := V_n(y1,1) (so W_n mod q = beta^n+beta^-n), direct
// expansion gives, for any a, b:
//
//   W_a - W_b = beta^-a * (beta^(a-b) - 1) * (beta^(a+b) - 1)
//
// so W_a == W_b (mod q) iff ord(beta) | (a-b) OR (a+b) -- V's built-in
// evenness (V_-n=V_n) gives the same "catches either candidate" property P-1
// manufactures by squaring, but on LINEAR indices. See Gpu.h's fuller
// derivation above Gpu::pp1Stage2.
// ---------------------------------------------------------------------------
namespace {

// V_n(P,1) mod M_p, by ordinary square-and-multiply exponentiation of r in
// the quadratic ring Z[r]/(r^2-Pr+1) -- structurally unrelated to
// Gpu::lucasVResidue's paired Montgomery ladder, so a bug in that recurrence
// cannot be mirrored here and cancel out. r^n = a+br; V_n = trace(r^n) =
// 2a+Pb, since the two roots of x^2-Px+1 sum to P and multiply to 1 (matching
// Q=1).
Nat vRef(u32 p, const Nat& P, u64 n) {
  if (n == 0) { return Nat{2}; }
  const Nat M = mersenne(p);
  auto subM = [&](const Nat& x, const Nat& y) { return modMersenneFast(sub(add(x, M), y), p); };

  Nat a{1}, b{};   // r^0 = 1
  for (int bit = 63 - std::countl_zero(n); bit >= 0; --bit) {
    // square: (a+br)^2 = (a^2-b^2) + (2ab + b^2 P) r
    Nat a2 = mulModM(a, a, p), b2 = mulModM(b, b, p), ab = mulModM(a, b, p);
    Nat newA = subM(a2, b2);
    Nat newB = modMersenneFast(add(add(ab, ab), mulModM(b2, P, p)), p);
    a = newA; b = newB;
    if ((n >> bit) & 1) {
      // multiply by r: (a+br)*r = -b + (a+bP) r
      Nat na = subM(Nat{}, b);
      Nat nb = modMersenneFast(add(a, mulModM(b, P, p)), p);
      a = na; b = nb;
    }
  }
  return modMersenneFast(add(add(a, a), mulModM(b, P, p)), p);
}

// Full P+1 stage-2 accumulator, independent of Gpu::pp1Stage2's incremental
// recurrences: every T_j and every A_m is a FRESH call to vRef, never built by
// stepping from a neighbour -- matching M6b's own "fresh exponentiation per
// block" standard for A_m, and going a step further for T_j (P-1's own
// reference still builds T_j incrementally; this doesn't need to, since vRef
// is cheap enough per entry at selftest scale).
Nat pp1Stage2Reference(u32 p, const Nat& y1, const Stage2Plan& plan) {
  const Nat M = mersenne(p);
  vector<Nat> T(plan.jset.size());
  for (size_t i = 0; i < plan.jset.size(); ++i) { T[i] = vRef(p, y1, plan.jset[i]); }

  Nat acc{1};
  for (u64 m = plan.mFirst; m <= plan.mLast; ++m) {
    bool any = false;
    for (size_t i = 0; i < plan.jset.size() && !any; ++i) { any = plan.slot(m, i); }
    if (!any) { continue; }
    const Nat A = vRef(p, y1, m * plan.d);
    for (size_t i = 0; i < plan.jset.size(); ++i) {
      if (!plan.slot(m, i)) { continue; }
      acc = mulModM(acc, sub(add(A, M), T[i]), p);
    }
  }
  return modMersenneFast(acc, p);
}

// V_n(P,1) mod a SMALL prime q, plain u64 arithmetic. Used only to check the
// algebraic identity itself (gate A) -- no GPU, no BigInt, and no relation to
// how the GPU or vRef compute it, since this gate isn't cross-validating an
// implementation, it's checking that the underlying theorem is true.
u64 vSeqSmall(u64 P, u64 n, u64 mod) {
  if (n == 0) { return 2 % mod; }
  const size_t nBits = std::bit_width(n);
  u64 A = P % mod;
  u64 B = (A * A % mod + mod - 2 % mod) % mod;
  for (size_t bi = nBits - 1; bi-- > 0; ) {
    const bool bit = (n >> bi) & 1;
    const u64 T = ((A * B) % mod + mod - P % mod) % mod;
    if (bit) {
      B = (B * B % mod + mod - 2 % mod) % mod;
      A = T;
    } else {
      A = (A * A % mod + mod - 2 % mod) % mod;
      B = T;
    }
  }
  return A;
}

} // namespace

int runPp1Stage2Tests(GpuCommon shared, Queue* q, const string& fftSpec) {
  printf("M7: P+1 stage 2\n\n");
  const int before = failures;

  printf("  A. W_a == W_b (mod q)  <=>  ord(beta) | (a-b) or (a+b)   (no GPU)\n");
  {
    int localFails = 0, localChecks = 0;
    std::mt19937_64 rng(12345);
    for (u64 qq : {13ull, 17ull, 19ull, 23ull, 29ull, 31ull, 37ull, 41ull, 43ull, 47ull, 101ull, 997ull, 7919ull}) {
      u64 y1 = 0, ordk = 0;
      for (u64 cand = 3; cand < qq && !ordk; ++cand) {
        u64 disc = (cand * cand % qq + qq - 4 % qq) % qq;
        if (disc == 0) { continue; }
        // Euler's criterion: disc^((q-1)/2) mod q == q-1 means non-residue,
        // which is exactly the condition P+1 needs for this base to apply.
        u64 e = (qq - 1) / 2, base = disc, r = 1;
        while (e) { if (e & 1) { r = r * base % qq; } base = base * base % qq; e >>= 1; }
        if (r != qq - 1) { continue; }
        for (u64 k = 1; k <= qq + 1; ++k) {
          if (vSeqSmall(cand, k, qq) == 2) { y1 = cand; ordk = k; break; }
        }
      }
      if (!ordk) { continue; }
      for (int t = 0; t < 30; ++t) {
        const u64 a = rng() % (3 * ordk), b = rng() % (3 * ordk);
        const u64 Wa = vSeqSmall(y1, a, qq), Wb = vSeqSmall(y1, b, qq);
        const bool lhs = (Wa == Wb);
        const u64 diffab = (a >= b) ? (a - b) : (b - a);
        const bool rhs = (diffab % ordk == 0) || ((a + b) % ordk == 0);
        ++localChecks;
        if (lhs != rhs) {
          ++localFails;
          if (localFails <= 3) {
            printf("     FAIL q=%llu y1=%llu ord=%llu a=%llu b=%llu Wa=%llu Wb=%llu\n",
                   (unsigned long long) qq, (unsigned long long) y1, (unsigned long long) ordk,
                   (unsigned long long) a, (unsigned long long) b,
                   (unsigned long long) Wa, (unsigned long long) Wb);
          }
        }
      }
    }
    check(localFails == 0, "identity holds across every (q,a,b) check");
    printf("     %s  %d checks over several primes, %d failed\n",
           localFails ? "FAIL" : "PASS", localChecks, localFails);
  }

  printf("\n  B. lucasVResidue vs an independent CPU reference (ring exponentiation)\n");
  {
    const u32 p = 859433;
    const Nat P{5};
    const Words base = makeWords(p, 5);
    FFTConfig fft = smallestFFT(p, fftSpec);
    auto gpu = Gpu::make(q, p, shared, fft, {}, false);

    for (u64 n : {0ull, 1ull, 2ull, 3ull, 4ull, 7ull, 8ull, 1000ull, 123457ull}) {
      const Words got = gpu->lucasVResidue(base, n);
      const Nat want = vRef(p, P, n);
      const Nat have = mod(fromWords(got), mersenne(p));
      const bool ok = (have == want);
      check(ok, "lucasVResidue at n=" + to_string(n));
      printf("     %s  n=%-8llu  %s\n", ok ? "PASS" : "FAIL", (unsigned long long) n,
             ok ? "match" : ("MISMATCH got " + have.hex() + " want " + want.hex()).c_str());
    }
  }

  printf("\n  C. full pp1Stage2 walk vs the CPU reference\n");
  {
    struct Case { u32 exponent; u64 b1; u64 b2; u32 d; u32 w; };
    static const Case CASES[] = {
      {859433, 1000, 1800, 210, 1},   // no window: exactly one slot per prime
      {859433, 1000, 1800, 210, 3},   // window: several candidate slots per prime
      {859433, 2000, 3000, 420, 1},
    };
    for (const Case& c : CASES) {
      Timer t;
      const Stage2Plan plan = buildStage2Plan(c.b1, c.b2, c.d, c.w);
      FFTConfig fft = smallestFFT(c.exponent, fftSpec);
      auto gpu = Gpu::make(q, c.exponent, shared, fft, {}, false);
      const Words y1 = makeWords(c.exponent, 5);

      const Words got = gpu->pp1Stage2(y1, plan, 0, [](u64, u64) { return true; });
      const Nat want = pp1Stage2Reference(c.exponent, Nat{5}, plan);
      const Nat have = mod(fromWords(got), mersenne(c.exponent));

      const bool ok = (have == want) && !want.isZero();
      check(ok, "pp1Stage2 accumulator matches the CPU reference");
      printf("     %s  M%-6u B1=%-5llu B2=%-5llu D=%-4u w=%u  %5llu slots  %s  (%s)\n",
             ok ? "PASS" : "FAIL", c.exponent, (unsigned long long) c.b1,
             (unsigned long long) c.b2, c.d, c.w,
             (unsigned long long) plan.nSlots,
             ok ? "match" : ("MISMATCH got " + have.hex() + " want " + want.hex()).c_str(),
             fmtDuration(t.at()).c_str());
    }
  }

  printf("\n  D. renormalize() vs skipping it, at a realistic exponent\n");
  {
    // The naive worry, when this was designed: an uncarried "one extra bit"
    // width would compound across the block-step/T-table's persisted state,
    // since each step's kSubWords result feeds the next step. Measuring says
    // otherwise -- see Gpu::renormalize's comment for why. Confirmed here at
    // a modest scale (fast enough to run every time); confirmed again during
    // development up to ~100000 accumulator slots at this same exponent, with
    // the same result, recorded in LESSONS_LEARNED.md rather than run here
    // every time.
    const u32 p = 86599237;
    try {
      FFTConfig fft = smallestFFT(p, fftSpec);
      auto gpu = Gpu::make(q, p, shared, fft, {}, false);
      const Stage2Plan plan = buildStage2Plan(1000000, 1005000, 210, 3);
      auto never = [](u64, u64) { return true; };
      const Words y1 = makeWords(p, 5);

      const Words correct = gpu->pp1Stage2(y1, plan, 0, never);
      const Words broken = gpu->pp1Stage2(y1, plan, 0, never, nullptr, false, nullptr, 0, {},
                                          nullptr, nullptr, /*skipRenormalize=*/true);

      const bool identical = (correct == broken);
      check(identical, "renormalize() and skipping it produce the same accumulator");
      printf("     %s  M%u  %llu slots  res64=%016llx  %s\n",
             identical ? "PASS" : "FAIL", p, (unsigned long long) plan.nSlots,
             (unsigned long long) residue(correct),
             identical ? "identical, as measured during development"
                       : "DIFFER -- the precaution turns out to matter after all, investigate before removing it");

      // choosePP1Bounds (Bounds.cpp) reuses cost.mulsPerPrime/s2MulCost --
      // measured for P-1's own stage2 -- on the grounds that pp1Stage2 reuses
      // the identical (m,j) pairing geometry and modMul-based accumulator
      // recurrence (see the comment above Gpu::pp1Stage2). Measured the same
      // way M6b measures P-1's own ratio: informational, not asserted, since
      // hardware/driver variance makes an exact match unrealistic to demand --
      // this is evidence the reuse is reasonable, not a correctness gate.
      const Stage2Plan small = buildStage2Plan(1000000, 1010000, 210, 1);
      const Stage2Plan big   = buildStage2Plan(1000000, 1100000, 210, 1);
      auto mulsOf = [](const Stage2Plan& q2) { return q2.nSlots + 2 * (q2.nBlocks() - 1); };

      Timer pt;
      gpu->pp1Stage2(y1, small, 0, never);
      const double ptSmall = pt.at();
      gpu->pp1Stage2(y1, big, 0, never);
      const double ptBig = pt.at();
      const double usPerMul = (ptBig - ptSmall) * 1e6 / double(mulsOf(big) - mulsOf(small));

      const u32 nSq = 2000;
      Timer ts;
      gpu->expExp2(makeWords(p, 3), nSq);
      const double usPerSquaring = ts.at() * 1e6 / nSq;

      printf("     %.0f us per P+1 accumulator multiply (marginal, %llu vs %llu muls)\n",
             usPerMul, (unsigned long long) mulsOf(big), (unsigned long long) mulsOf(small));
      printf("     %.0f us per stage-1 squaring on the same transform\n", usPerSquaring);
      printf("     ratio %.2fx (P-1's own measured ratio is %.2fx, from M6b) -- reusing\n"
             "     P-1's s2MulCost for P+1's own bounds model\n",
             usPerMul / usPerSquaring, CostModel{}.s2MulCost);
    } catch (const char* s) {
      printf("     SKIP (%s)\n", s);
    } catch (const std::exception& e) {
      printf("     SKIP (%s)\n", e.what());
    }
  }

  printf("\n  E. interrupt and resume reproduce an uninterrupted walk exactly\n");
  {
    const u32 p = 859433;
    FFTConfig fft = smallestFFT(p, fftSpec);
    auto gpu = Gpu::make(q, p, shared, fft, {}, false);
    const Stage2Plan plan = buildStage2Plan(2000, 3000, 420, 1);
    const Words y1 = makeWords(p, 5);
    auto always = [](u64, u64) { return true; };

    const Words full = gpu->pp1Stage2(y1, plan, 0, always);

    // Interrupt partway: stop about halfway through the slots (reportEvery=1
    // so the stop point is exact, not rounded to the next report boundary).
    const u64 stopAfter = std::max<u64>(1, plan.nSlots / 2);
    u64 seen = 0;
    Gpu::Stage2Pos stoppedAt;
    auto stopHalfway = [&](u64 done, u64) { seen = done; return done < stopAfter; };
    gpu->pp1Stage2(y1, plan, 1, stopHalfway, nullptr, false, nullptr, 0, {}, &stoppedAt);

    const bool stoppedEarly = !stoppedAt.acc.empty() && seen >= stopAfter;
    check(stoppedEarly, "the walk actually stopped before completion");

    const Words resumed = gpu->pp1Stage2(y1, plan, 0, always, nullptr, false, &stoppedAt);

    const bool ok = stoppedEarly && (resumed == full);
    check(ok, "resumed accumulator matches an uninterrupted run");
    printf("     %s  M%u  stopped at %llu/%llu, resumed  full res64=%016llx  resumed res64=%016llx\n",
           ok ? "PASS" : "FAIL", p, (unsigned long long) seen,
           (unsigned long long) (plan.nSlots + 2 * (plan.nBlocks() - 1)),
           (unsigned long long) residue(full), (unsigned long long) residue(resumed));
  }

  printf("\n  F. checkpoint format: round-trip and rejection rules   (CPU)\n");
  {
    const u32 p = 1;   // cannot collide with a real save file
    const u64 b1 = 100, b2 = 200;
    const u32 seed = 3;
    const u64 yr = 0x0123456789abcdefull;
    const string path = defaultPp1Stage2Path(p, b1, b2, seed);

    Pp1Stage2State st;
    st.exponent = p; st.b1 = b1; st.b2 = b2; st.d = 210; st.w = 1; st.seed = seed;
    st.yRes64 = yr;
    st.acc = Words{1, 2, 3, 4};
    st.a = Words{5, 6};
    st.s = Words{7, 8};

    string err;
    check(savePp1Stage2(path, st, err), "write a checkpoint");

    Pp1Stage2State want = st;
    Pp1Stage2State got;
    const bool ok1 = loadPp1Stage2(path, want, got, err) && got.acc == st.acc;
    check(ok1, "round-trip: matching fields load back identically");
    printf("     %s  round-trip\n", ok1 ? "PASS" : "FAIL");

    struct Case { const char* what; u32 seed; u64 yr; };
    static const Case CASES[] = {
      {"different seed",            5, yr},
      {"different stage-1 residue", 3, 0xdeadbeefdeadbeefull},
    };
    for (const Case& c : CASES) {
      Pp1Stage2State w2 = want;
      w2.seed = c.seed;
      w2.yRes64 = c.yr;
      Pp1Stage2State g2;
      string e2;
      const bool gotOk = loadPp1Stage2(path, w2, g2, e2);
      check(!gotOk, string("checkpoint rejection: ") + c.what);
      const string outcome = gotOk ? string("ACCEPTED") : "rejected (" + e2 + ")";
      printf("     %s  %-28s -> %s\n", gotOk ? "FAIL" : "PASS", c.what, outcome.c_str());
    }

    std::error_code ec;
    filesystem::remove(path, ec);
  }

  printf("\n  G. a real factor whose missing prime is in the gap   (M862907)\n");
  {
    // From an independent search (trial-divide M_p for q=2kp+1, factor q+1),
    // filtered to exclude p with ANY OTHER q'=2k'p+1 | M_p whose k' is itself
    // B1-smooth: since stage1Exponent's shared E bakes in a factor of p, and
    // q'-1 = 2k'p always contains p trivially, such a q' is caught by stage 1
    // regardless of seed the moment k' is B1-smooth -- the first candidate
    // tried (M881917) turned out to have exactly this kind of unrelated,
    // easier factor, caught here before it could look like a bug.
    //
    // q=125984423 = 2*73*862907+1 divides M862907. q+1 = 2^3*3*523*10037, so
    // k=73 is 523-smooth apart from the single prime 10037, and seed=3 gives
    // (3^2-4)=5 a quadratic non-residue mod q (required for P+1 to apply with
    // this seed). So B1=600 must NOT find it, and B2=10100 must.
    const u32 p = 862907;
    const char* factor = "125984423";
    const u32 seed = 3;
    const u64 b1 = 600, b2 = 10100;

    Nat qFactor;
    check(fromDecimal(factor, qFactor), "parse expected factor");

    Config cfg;
    cfg.exponent = p;
    cfg.fftSpec = fftSpec;
    cfg.reportEvery = 0;
    cfg.checkpoint = false;             // hermetic: leave no save files behind

    FFTConfig fft = smallestFFT(p, fftSpec);
    auto gpu = Gpu::make(q, p, shared, fft, {}, false);

    Timer t1;
    PP1Result pr = runPP1Stage1(*gpu, cfg, b1, seed, false);
    const bool foundAtStage1 = pr.foundFactor;
    check(!foundAtStage1, "the factor is NOT reachable by stage 1 alone at B1=600");
    printf("     %s  stage 1 to B1=%llu: found early? %-3s  (%s)\n",
           foundAtStage1 ? "FAIL" : "PASS", (unsigned long long) b1,
           foundAtStage1 ? "YES" : "no", fmtDuration(t1.at()).c_str());

    Timer t2;
    PP1Stage2Result s2 = runPP1Stage2(*gpu, cfg, pr.residue, seed, b1, b2,
                                      /*d=*/210, /*w=*/1, false);
    const bool foundAtStage2 = s2.foundFactor &&
        std::any_of(s2.factors.begin(), s2.factors.end(),
                    [&](const FoundFactor& ff) { return cmp(ff.value, qFactor) == 0; });
    check(foundAtStage2, "stage 2 to B2=10100 finds the factor");
    printf("     %s  stage 2 to B2=%llu: found? %-3s  (%s)\n",
           foundAtStage2 ? "PASS" : "FAIL", (unsigned long long) b2,
           foundAtStage2 ? "yes" : "NO", fmtDuration(t2.at()).c_str());
  }

  printf("\nM7: %d failed.\n\n", failures - before);
  return failures == before ? 0 : 1;
}
