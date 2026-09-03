// Copyright (C) Mp_p-1_gpu

#include "PM1.h"

#include "Gcd.h"
#include "Save.h"
#include "Stage2Plan.h"
#include "Stage2Save.h"
#include "Pp1Stage2Save.h"
#include "Gpu.h"
#include "Primes.h"
#include "timeutil.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <io.h>

#include <atomic>
#include <cstdio>
#include <vector>

using namespace std;

// Set by the Ctrl-C handler installed in main.
extern std::atomic<bool> gInterrupted;

// True when stdout is a console rather than a file or a pipe. The in-place
// progress line is right for a terminal and useless in a log, where it produces
// one enormous line of \r-joined updates, so the two are formatted differently.
static bool stdoutIsTerminal() {
  static const bool tty = _isatty(_fileno(stdout)) != 0;
  return tty;
}

// How many phases this job has: 3 for stage 1 alone (build E, squarings, gcd),
// 5 when a stage 2 follows (plus the stage-2 walk and its own gcd). Set once by
// the driver before anything prints.
u32 gPhaseTotal = 3;

// "2/3" or "2/5". Returns a rotating buffer so several can appear in one
// printf; three live at once is more than any call site needs.
static const char* ph(u32 n) {
  static char bufs[4][12];
  static int k = 0;
  k = (k + 1) & 3;
  snprintf(bufs[k], sizeof(bufs[k]), "%u/%u", n, gPhaseTotal);
  return bufs[k];
}

// gcd(v, M_p), reporting as it goes. Shared by stage 1 and stage 2 -- at these
// sizes it runs for 10-20 minutes and would otherwise print nothing at all,
// which looks exactly like a hang.
// `announce` off means print nothing at all -- not the header, not the closing
// "done in". The overlapped stage-1 gcd needs that: it runs while stage 2 owns
// the terminal, so its lines would land in the middle of stage 2's, and its
// header would claim the GPU is idle when stage 2 is in fact driving it. The
// caller prints the summary once both have finished.
// `interrupted` is set true, and the returned Nat is meaningless (Nat{}, not
// a real gcd), if gInterrupted became set partway through: the ONLY place in
// a run that could previously not be stopped by Ctrl-C, since gcdHalf's
// progress hook used to be void-returning and had nothing to check. Every
// caller must look at `interrupted` before touching the return value -- an
// aborted gcd is not "no factor found," it is "no answer at all."
static Nat gcdWithProgress(const Nat& v, const Nat& Mp, u32 phase, const char* what,
                           bool showProgress, double& secsOut, bool& interrupted,
                           bool announce = true) {
  interrupted = false;
  if (announce) {
    printf("  [%s gcd CPU] gcd(%s, M_p), %zu bits -- CPU only, GPU idle from here\n",
           ph(phase), what, Mp.bits());
    fflush(stdout);
  }
  Timer gcdTimer;

  // gcd() is gcdGmp now (see Gcd.h's file comment) -- a single opaque call
  // with no progress callback, so there is no per-tick ETA line to show
  // (showProgress no longer has anything to gate here; it still reaches this
  // function because callers thread it through for their own other prints)
  // and no way to throw GcdAborted out of a call already in flight. The best
  // available check is immediately before the call: if Ctrl-C landed while
  // this phase was starting up, honor it before spending gcdGmp's ~12s
  // production-scale worst case on an answer that will be thrown away
  // anyway. See gcdGmp's comment in Gcd.cpp for why that worst case is small
  // enough now that running one already-started call to completion is an
  // acceptable trade.
  Nat g;
  if (gInterrupted.load()) {
    interrupted = true;
  } else {
    g = v.isZero() ? Mp : gcd(v, Mp);
  }

  secsOut = gcdTimer.at();
  if (announce) {
    if (interrupted) {
      printf("  [%s gcd CPU] interrupted after %s\n", ph(phase), fmtDuration(secsOut).c_str());
    } else {
      printf("  [%s gcd CPU] done in %s\n", ph(phase), fmtDuration(secsOut).c_str());
    }
  }
  return g;
}

string fmtDuration(double s) {
  if (!(s >= 0.0) || s > 3.1e9) { return "--"; }
  unsigned long long t = (unsigned long long) (s + 0.5);
  const unsigned long long h = t / 3600, m = (t % 3600) / 60, sec = t % 60;
  char buf[64];
  if (h >= 24)  { snprintf(buf, sizeof(buf), "%llud%02lluh", h / 24, h % 24); }
  else if (h)   { snprintf(buf, sizeof(buf), "%lluh%02llum", h, m); }
  else if (m)   { snprintf(buf, sizeof(buf), "%llum%02llus", m, sec); }
  else          { snprintf(buf, sizeof(buf), "%llus", sec); }
  return string(buf);
}

namespace {

// gpuowl's src/Primes.h is a small fixed sieve for prev/next-prime queries and
// cannot enumerate; a plain sieve of Eratosthenes is what is wanted here.
vector<u32> primesUpTo(u32 limit) {
  vector<u32> out;
  if (limit < 2) { return out; }
  vector<bool> composite(size_t(limit) + 1, false);
  for (u32 i = 2; u64(i) * i <= limit; ++i) {
    if (!composite[i]) {
      for (u64 j = u64(i) * i; j <= limit; j += i) { composite[size_t(j)] = true; }
    }
  }
  for (u32 i = 2; i <= limit; ++i) { if (!composite[i]) { out.push_back(i); } }
  return out;
}

} // namespace

namespace {

FoundFactor describe(const Nat& v, u32 exponent) {
  FoundFactor f;
  f.value = v;
  f.prime = isProbablePrime(v);
  // A genuine factor of M_p satisfies 2^p == 1 (mod v). Checking it costs
  // microseconds and catches any error in the gcd or the splitting.
  f.dividesMp = powMod(Nat(2), Nat(exponent), v).isOne();
  // k = (v-1)/(2p), reported because it is what P-1's bounds actually act on.
  const Nat km1 = sub(v, Nat(1));
  const u64 twoP = 2ull * exponent;
  if (modU64(km1, twoP) == 0) {
    Nat q, r;
    divrem(km1, Nat(twoP), q, r);
    if (q.size() <= 1) { f.k = q[0]; }
  }
  return f;
}

} // namespace

std::vector<FoundFactor> splitFactors(const Nat& g, u32 exponent, u64 maxK) {
  std::vector<FoundFactor> out;
  if (g.isZero() || g.isOne()) { return out; }

  Nat rem = g;
  const u64 twoP = 2ull * exponent;

  // Every prime factor of M_p is 2kp+1, so walking k tries only the ~1/(2p)
  // of integers that could possibly divide. Ascending k also means the first
  // divisor found is the smallest, hence prime.
  for (u64 k = 1; k <= maxK && !rem.isOne(); ++k) {
    if (k > (~u64(0) - 1) / twoP) { break; }        // 2kp+1 would overflow
    const u64 cand = twoP * k + 1;
    if (modU64(rem, cand) != 0) { continue; }

    const Nat c(cand);
    while (modU64(rem, cand) == 0) {
      Nat q, r;
      divrem(rem, c, q, r);
      rem = std::move(q);
      out.push_back(describe(c, exponent));
    }
  }

  // Whatever is left: prime (report it) or an unsplit composite (say so).
  if (!rem.isOne() && !rem.isZero()) { out.push_back(describe(rem, exponent)); }
  return out;
}

Nat stage1Exponent(u64 b1, u32 exponent) {
  // E = 2*p * prod over primes q <= B1 of q^floor(log_q B1).
  //
  // The 2*p is essential, not an optimisation. Every factor of M_p has the form
  // q = 2kp+1, so ord_q(a) divides q-1 = 2kp -- and for a generic base it does
  // contain p. (Measured on the M86599237 test vector: 3 is a primitive root
  // mod q, so ord_q(3) = q-1 and p | ord_q(3).) P-1 succeeds only when
  // ord_q(a) | E, so E must carry the p and the 2 as well as the B1-smooth part
  // of k. p is known in advance, so this costs ~log2(p) extra bits -- 27 out of
  // ~1.4 million at realistic bounds -- and without it the algorithm cannot
  // find anything.
  const vector<u32> primes = primesUpTo(u32(b1));

  // Combine with a balanced product tree. Folding one small factor at a time
  // into a growing accumulator would be quadratic; pairing equal-sized pieces
  // keeps Karatsuba doing useful work.
  vector<Nat> stack;
  for (u32 q : primes) {
    u64 power = q;
    while (power <= b1 / q) { power *= q; }
    stack.push_back(Nat(power));
    // Merge while the top two are comparable in size.
    while (stack.size() >= 2 &&
           stack[stack.size() - 1].size() >= stack[stack.size() - 2].size()) {
      Nat a = std::move(stack.back()); stack.pop_back();
      Nat b = std::move(stack.back()); stack.pop_back();
      stack.push_back(mul(a, b));
    }
  }

  Nat e(1);
  for (Nat& part : stack) { e = mul(e, part); }

  // Fold in the 2*p that every factor's order can contain.
  e = mul(e, Nat(2));
  e = mul(e, Nat(exponent));
  return e;
}

std::vector<u64> findCompletedB1(u32 exponent, u64 below) {
  // Completed checkpoints are named pm1_<exponent>_b1_<B1>.save. Return the
  // usable ones (B1 < below) largest first, so the cheapest extension wins.
  std::vector<u64> out;
  const string prefix = "pm1_" + to_string(exponent) + "_b1_";
  std::error_code ec;
  for (const auto& e : std::filesystem::directory_iterator(".", ec)) {
    const string name = e.path().filename().string();
    if (name.rfind(prefix, 0) != 0 || name.size() <= prefix.size() + 5) { continue; }
    if (name.substr(name.size() - 5) != ".save") { continue; }
    const string mid = name.substr(prefix.size(), name.size() - prefix.size() - 5);
    char* end = nullptr;
    const u64 v = strtoull(mid.c_str(), &end, 10);
    if (end && *end == 0 && v && v < below) { out.push_back(v); }
  }
  std::sort(out.begin(), out.end(), std::greater<u64>());
  return out;
}

Nat stage1ExponentDelta(u64 b1From, u64 b1To, u32 exponent) {
  // R such that E(b1To) == E(b1From) * R.
  //
  // Computed directly rather than by dividing the two exponents. E(b1To) is
  // ~14M bits at realistic bounds and Knuth-D division of that by a 4M-bit
  // value costs tens of seconds; building R from the prime powers is exact and
  // takes well under a second.
  //
  // For each prime q <= b1To the exponent of q in E is floor(log_q B1), which
  // is non-decreasing in B1 -- that is precisely why E(b1From) divides E(b1To).
  // R therefore carries q^(new - old); for q > b1From the old exponent is 0.
  // The 2*p factor is already in the old E and cancels, so it is NOT repeated.
  assert(b1To > b1From);
  const vector<u32> primes = primesUpTo(u32(b1To));

  vector<Nat> stack;
  for (u32 q : primes) {
    u64 newPow = q;
    while (newPow <= b1To / q) { newPow *= q; }

    u64 oldPow = 1;
    if (u64(q) <= b1From) {
      oldPow = q;
      while (oldPow <= b1From / q) { oldPow *= q; }
    }
    if (newPow == oldPow) { continue; }          // this prime contributes nothing

    assert(newPow % oldPow == 0);
    stack.push_back(Nat(newPow / oldPow));
    while (stack.size() >= 2 &&
           stack[stack.size() - 1].size() >= stack[stack.size() - 2].size()) {
      Nat a = std::move(stack.back()); stack.pop_back();
      Nat b = std::move(stack.back()); stack.pop_back();
      stack.push_back(mul(a, b));
    }
  }

  Nat r(1);
  for (Nat& part : stack) { r = mul(r, part); }
  return r;
}

PM1Result runPM1Stage1(Gpu& gpu, const Config& cfg, u64 b1, bool showProgress,
                       bool doGcd) {
  PM1Result res;
  res.b1Used = b1;

  printf("  [%s exponent CPU] building E for B1 = %llu ... ", ph(1), (unsigned long long) b1);
  fflush(stdout);
  Timer buildTimer;
  const Nat E = stage1Exponent(b1, cfg.exponent);
  printf("%zu bits (%s)\n", E.bits(), fmtDuration(buildTimer.at()).c_str());

  // The ladder wants little-endian 64-bit limbs.
  const vector<u64> limbs = E.toVector();

  Timer stage1;
  Timer sinceReport;
  double lastAt = 0;

  // Squarings already done before this process started, on a resumed run. Rate
  // and ETA must be based on work done in THIS run: counting resumed squarings
  // against this run's elapsed time inflates the rate and makes the ETA look
  // far better than it is.
  u64 doneAtStart = 0;

  auto progress = [&](u64 done, u64 total) -> bool {
    if (showProgress && total) {
      const double el = stage1.at();
      const u64 thisRun = done > doneAtStart ? done - doneAtStart : 0;
      const double rate = (thisRun > 0 && el > 0) ? thisRun / el : 0;
      const double eta = rate > 0 ? (total - done) / rate : 0;
      if (stdoutIsTerminal()) {
        printf("\r  [%s stage 1 GPU] %6.2f%%  %llu/%llu squarings (%.0f/s)  elapsed %s  ETA %s   ",
               ph(2), 100.0 * double(done) / double(total),
               (unsigned long long) done, (unsigned long long) total,
               rate, fmtDuration(el).c_str(), fmtDuration(eta).c_str());
      } else if (el - lastAt >= 60.0 || done == total) {
        // Redirected to a file: one ordinary newline-terminated line a minute.
        // The in-place \r form turns a log into one unreadable mega-line.
        lastAt = el;
        printf("  [%s stage 1 GPU] %6.2f%%  %llu/%llu squarings (%.0f/s)  elapsed %s  ETA %s\n",
               ph(2), 100.0 * double(done) / double(total),
               (unsigned long long) done, (unsigned long long) total,
               rate, fmtDuration(el).c_str(), fmtDuration(eta).c_str());
      }
      fflush(stdout);
    }
    (void) sinceReport;
    return !gInterrupted.load();
  };

  // ---- checkpointing ------------------------------------------------------
  const string savePath = cfg.checkpointFile.empty()
                        ? defaultSavePath(cfg.exponent, b1) : cfg.checkpointFile;

  SaveState want;
  want.exponent = cfg.exponent;
  want.b1 = b1;
  want.eBits = E.bits();
  want.base = 3;
  want.eVersion = E_FORMAT_VERSION;

  SaveState loaded;
  const Words* resumeFrom = nullptr;
  u64 resumeBit = 0;
  bool alreadyComplete = false;

  if (cfg.checkpoint) {
    string err;
    if (loadState(savePath, want, loaded, err)) {
      if (loaded.complete) {
        alreadyComplete = true;
        printf("  [%s stage 1 GPU] already complete in %s -- skipping to the gcd\n", ph(2),
               savePath.c_str());
      } else {
        resumeFrom = &loaded.residue;
        resumeBit = loaded.nextBit;
        // So the rate/ETA reflect only what this run computes.
        doneAtStart = (E.bits() - 1) - loaded.nextBit;
        const double pct = 100.0 * double(want.eBits - 1 - resumeBit) / double(want.eBits - 1);
        printf("  [%s stage 1 GPU] resuming from %s at %.2f%% (bit %llu of %llu)\n", ph(2),
               savePath.c_str(), pct,
               (unsigned long long) resumeBit, (unsigned long long) (want.eBits - 1));
      }
    } else if (err != "no checkpoint") {
      // Anything other than "absent" means a file exists but cannot be trusted.
      printf("  [%s stage 1 GPU] ignoring checkpoint: %s\n", ph(2), err.c_str());
    }
  }

  // Checkpoint interval is configured in seconds, but the ladder counts
  // squarings, so convert once a rate is known. Until then use a rough guess
  // and let the first interval correct it.
  u32 saveEverySquarings = 0;
  double lastSaveAt = 0;
  auto saveFn = [&](const Words& residue, u64 nextBit) {
    if (!cfg.checkpoint) { return; }
    SaveState s = want;
    s.nextBit = nextBit;
    s.complete = false;
    s.residue = residue;
    string err;
    if (!saveState(savePath, s, err)) {
      printf("\n  WARNING: checkpoint failed: %s\n", err.c_str());
    }
    lastSaveAt = stage1.at();
  };

  if (cfg.checkpoint && cfg.checkpointSeconds) {
    // ~1400 us/squaring is the measured rate on this class of GPU; being off
    // only changes how often we save, never correctness.
    saveEverySquarings = u32(std::max<double>(1000, cfg.checkpointSeconds * 700.0));
  }

  // ---- B1 extension -------------------------------------------------------
  // If a SMALLER B1 already completed for this exponent, E(oldB1) divides
  // E(b1), so x_new = x_old ^ (E(b1)/E(oldB1)) -- only the new work is done.
  // The saving is exactly oldB1/b1.
  Words extBase;
  u64 extendFrom = 0;
  if (cfg.checkpoint && !alreadyComplete && !resumeFrom && cfg.extend) {
    for (u64 cand : findCompletedB1(cfg.exponent, b1)) {
      SaveState w2 = want;
      w2.b1 = cand;
      w2.eBits = stage1Exponent(cand, cfg.exponent).bits();
      SaveState got;
      string err;
      if (loadState(defaultSavePath(cfg.exponent, cand), w2, got, err) && got.complete) {
        extBase = got.residue;
        extendFrom = cand;
        break;                    // findCompletedB1 returns largest first
      }
    }
  }

  Words x;
  if (alreadyComplete) {
    x = loaded.residue;
    res.squarings = 0;
  } else if (extendFrom) {
    printf("  [%s stage 1 GPU] extending the completed B1=%llu result to B1=%llu\n", ph(2),
           (unsigned long long) extendFrom, (unsigned long long) b1);
    Timer rTimer;
    const Nat R = stage1ExponentDelta(extendFrom, b1, cfg.exponent);
    printf("  [%s exponent CPU] R = E(%llu)/E(%llu) is %zu bits, vs %zu from scratch"
           " -- %.0f%% saved (%s)\n",
           ph(1), (unsigned long long) b1, (unsigned long long) extendFrom,
           R.bits(), E.bits(),
           100.0 * (1.0 - double(R.bits()) / double(E.bits())),
           fmtDuration(rTimer.at()).c_str());

    want.baseB1 = extendFrom;
    x = gpu.powResidue(extBase, R.toVector(), cfg.reportEvery, progress,
                       nullptr, 0, saveEverySquarings, saveFn);
    res.squarings = R.bits() ? R.bits() - 1 : 0;
  } else {
    x = gpu.powBase3(limbs, cfg.reportEvery, progress,
                     resumeFrom, resumeBit, saveEverySquarings, saveFn);
    res.squarings = E.bits() ? E.bits() - 1 : 0;
  }
  res.stage1Secs = stage1.at();
  (void) lastSaveAt;

  if (showProgress && stdoutIsTerminal()) { printf("\r%*s\r", 100, ""); }

  if (gInterrupted.load()) {
    res.interrupted = true;
    return res;
  }

  // Stage 1 finished: record it as complete so a re-run (or a crash during the
  // gcd, which takes minutes) does not repeat hours of squaring.
  if (cfg.checkpoint && !alreadyComplete) {
    SaveState s = want;
    s.nextBit = 0;
    s.complete = true;
    s.residue = x;
    string err;
    if (!saveState(savePath, s, err)) {
      printf("  WARNING: could not record stage 1 as complete: %s\n", err.c_str());
    }
  }

  printf("  [%s stage 1 GPU] done in %s (%.0f us/squaring)\n", ph(2),
         fmtDuration(res.stage1Secs).c_str(),
         res.squarings ? res.stage1Secs * 1e6 / double(res.squarings) : 0.0);

  res.residue = x;                              // stage 2 continues from here
  Nat xn = fromWords(x);
  const Nat Mp = mersenne(cfg.exponent);
  if (xn.isZero()) { xn = Mp; }                 // x == 0 represents M_p
  res.xMinusOne = xn.isZero() ? Nat{} : sub(xn, Nat(1));

  // Callers checking a KNOWN factor only need (x-1) mod q and can stop here.
  if (!doGcd) { return res; }

  finishStage1Gcd(res, cfg.exponent, showProgress);
  return res;
}

// The stage-1 gcd, split out of runPM1Stage1 so the driver can run it on
// another thread while the GPU gets on with stage 2. Everything it needs is
// already in `res` by the time stage 1 returns (xMinusOne is always
// populated), and it writes only into `res`, so the two do not interact.
void finishStage1Gcd(PM1Result& res, u32 exponent, bool showProgress, bool announce) {
  const Nat Mp = mersenne(exponent);

  // factor = gcd(x - 1, M_p). This is CPU work. Historically the GPU was idle
  // throughout, which is why the phase is labelled; when the driver overlaps
  // it with stage 2 the GPU is busy instead and this prints nothing until it
  // finishes, to keep the two phases' output from interleaving.
  Nat g = gcdWithProgress(res.xMinusOne, Mp, 3, "x-1", showProgress, res.gcdSecs,
                          res.interrupted, announce);
  if (res.interrupted) { return; }

  // A gcd of 1 means no factor; a gcd of M_p means the whole number divided
  // out, which happens only for a degenerate B1 and is not a useful result.
  if (!g.isOne() && cmp(g, Mp) != 0) {
    res.foundFactor = true;
    res.gcdValue = g;
    res.factor = g;

    // The gcd carries EVERY factor whose k was smooth, multiplied together, so
    // it is routinely composite. Split it before reporting.
    printf("  splitting the gcd (%zu bits) into prime factors ... ", g.bits());
    fflush(stdout);
    Timer splitTimer;
    res.factors = splitFactors(g, exponent);
    printf("%zu found in %s\n", res.factors.size(),
           fmtDuration(splitTimer.at()).c_str());
  }
}

// ---------------------------------------------------------------------------
// P+1 stage 1
// ---------------------------------------------------------------------------
const u32 Pp1Start::SMALL_PRIMES[16] =
    {11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71};

Pp1Start Pp1Start::forRun(u32 nthRun, u32 exponent) {
  if (nthRun <= 1) { return {2, 7}; }
  if (nthRun == 2) { return {6, 5}; }
  // splitmix64 on (exponent, nthRun). Any decent mixer would do; what matters
  // is that it is FIXED, so the same assignment resumes with the same start.
  u64 h = (u64(exponent) << 32) ^ (u64(nthRun) * 0x9E3779B97F4A7C15ull);
  h ^= h >> 30; h *= 0xBF58476D1CE4E5B9ull;
  h ^= h >> 27; h *= 0x94D049BB133111EBull;
  h ^= h >> 31;
  // Prime95: numerator = 72 + (rand() & 0x7F), denominator = table[rand() & 0xF].
  return {u32(72 + (h & 0x7F)), SMALL_PRIMES[(h >> 8) & 0xF]};
}

PP1Result runPP1Stage1(Gpu& gpu, const Config& cfg, u64 b1, Pp1Start start,
                       bool showProgress) {
  PP1Result res;
  res.b1Used = b1;
  res.startUsed = start;

  printf("  [%s exponent CPU] building E for B1 = %llu ... ", ph(1), (unsigned long long) b1);
  fflush(stdout);
  Timer buildTimer;
  const Nat E = stage1Exponent(b1, cfg.exponent);
  printf("%zu bits in %s\n", E.bits(), fmtDuration(buildTimer.at()).c_str());

  Timer stage1;
  double lastAt = 0;
  u64 doneAtStart = 0;

  auto progress = [&](u64 done, u64 total) -> bool {
    if (showProgress && total) {
      const double el = stage1.at();
      const u64 thisRun = done > doneAtStart ? done - doneAtStart : 0;
      const double rate = (thisRun > 0 && el > 0) ? thisRun / el : 0;
      const double eta = rate > 0 ? (total - done) / rate : 0;
      if (stdoutIsTerminal()) {
        printf("\r  [%s P+1 start %s GPU] %6.2f%%  %llu/%llu steps (%.0f/s)"
               "  elapsed %s  ETA %s   ",
               ph(2), start.label().c_str(), 100.0 * double(done) / double(total),
               (unsigned long long) done, (unsigned long long) total,
               rate, fmtDuration(el).c_str(), fmtDuration(eta).c_str());
      } else if (el - lastAt >= 60.0 || done == total) {
        lastAt = el;
        printf("  [%s P+1 start %s GPU] %6.2f%%  %llu/%llu steps (%.0f/s)"
               "  elapsed %s  ETA %s\n",
               ph(2), start.label().c_str(), 100.0 * double(done) / double(total),
               (unsigned long long) done, (unsigned long long) total,
               rate, fmtDuration(el).c_str(), fmtDuration(eta).c_str());
      }
      fflush(stdout);
    }
    return !gInterrupted.load();
  };

  // Checkpoint under a start-specific name: different starting points are
  // different computations and their residues are not interchangeable. The
  // name changed shape in 1.9.4 (it was _s<seed>), so a P+1 run in flight
  // under an integer seed will not resume -- there is no start it could
  // resume AS, integer seeds no longer being used.
  char pathBuf[160];
  snprintf(pathBuf, sizeof(pathBuf), "pp1_%u_b1_%llu_n%ud%u.save",
           cfg.exponent, (unsigned long long) b1, start.num, start.den);
  const string savePath = pathBuf;

  SaveState want;
  want.exponent = cfg.exponent;
  want.b1 = b1;
  want.eBits = E.bits();
  want.base = 3;                 // unused by P+1, but must match on reload
  want.seed = start.id();
  want.eVersion = E_FORMAT_VERSION;

  SaveState loaded;
  const Words *rA = nullptr, *rB = nullptr;
  u64 resumeBit = 0;
  bool alreadyComplete = false;

  if (cfg.checkpoint) {
    string err;
    if (loadState(savePath, want, loaded, err)) {
      if (loaded.complete) {
        alreadyComplete = true;
        printf("  [%s P+1 start %s GPU] already complete -- skipping to the gcd\n",
               ph(2), start.label().c_str());
      } else if (!loaded.residue2.empty()) {
        rA = &loaded.residue;
        rB = &loaded.residue2;
        resumeBit = loaded.nextBit;
        doneAtStart = (E.bits() - 1) - resumeBit;
        printf("  [%s P+1 start %s GPU] resuming at %.2f%%\n", ph(2), start.label().c_str(),
               100.0 * double(doneAtStart) / double(E.bits() - 1));
      }
    } else if (err != "no checkpoint") {
      printf("  [%s P+1 start %s GPU] ignoring checkpoint: %s\n",
             ph(2), start.label().c_str(), err.c_str());
    }
  }

  auto saveFn = [&](const Words& a, const Words& b, u64 nextBit) {
    if (!cfg.checkpoint) { return; }
    SaveState s = want;
    s.nextBit = nextBit;
    s.complete = false;
    s.residue = a;
    s.residue2 = b;
    string err;
    if (!saveState(savePath, s, err)) {
      printf("\n  WARNING: checkpoint failed: %s\n", err.c_str());
    }
  };
  const u32 saveEverySteps = (cfg.checkpoint && cfg.checkpointSeconds)
      ? u32(std::max<double>(1000, cfg.checkpointSeconds * 350.0)) : 0;

  Words v;
  if (alreadyComplete) {
    v = loaded.residue;
  } else {
    // V_1 = num/den mod M_p. A rational start, per Prime95 -- see Pp1Start.
    // Cheap: ratioModMersenne avoids an extended GCD over p bits.
    const Words v1 = toWords(ratioModMersenne(start.num, start.den, cfg.exponent),
                             cfg.exponent);
    v = gpu.lucasVBase(v1, E.toVector(), cfg.reportEvery, progress, rA, rB, resumeBit,
                       saveEverySteps, saveFn);
    res.squarings = E.bits() ? E.bits() - 1 : 0;
  }
  res.stage1Secs = stage1.at();
  res.residue = v;             // stage 2 continues from here
  if (showProgress && stdoutIsTerminal()) { printf("\r%*s\r", 100, ""); }

  if (gInterrupted.load()) { res.interrupted = true; return res; }

  printf("  [%s P+1 start %s GPU] done in %s (%.0f us/step)\n", ph(2), start.label().c_str(),
         fmtDuration(res.stage1Secs).c_str(),
         res.squarings ? res.stage1Secs * 1e6 / double(res.squarings) : 0.0);

  if (cfg.checkpoint && !alreadyComplete) {
    SaveState s = want;
    s.nextBit = 0;
    s.complete = true;
    s.residue = v;
    string err;
    if (!saveState(savePath, s, err)) {
      printf("  WARNING: could not record P+1 stage 1 as complete: %s\n", err.c_str());
    }
  }

  // factor = gcd(V - 2, M_p)
  const Nat Mp = mersenne(cfg.exponent);
  Nat vn = fromWords(v);
  if (vn.isZero()) { vn = Mp; }
  Nat vm2;
  if (cmp(vn, Nat(2)) >= 0) {
    vm2 = sub(vn, Nat(2));
  } else {
    // V < 2: work modulo M_p rather than going negative.
    vm2 = sub(add(vn, Mp), Nat(2));
  }

  if (vm2.isZero()) {
    // V == 2 exactly: the sequence collapsed, so the gcd would be M_p itself.
    printf("  [%s gcd CPU] V == 2, nothing to extract for this start\n", ph(3));
    return res;
  }

  // Was its own hand-copied inline duplicate of gcdWithProgress -- which
  // meant it also had its own, separate copy of the bug where the gcd could
  // not be interrupted (see gcdWithProgress's own comment). Routed through
  // the shared helper now instead of fixing the duplicate in parallel.
  Nat g = gcdWithProgress(vm2, Mp, 3, "V-2", showProgress, res.gcdSecs, res.interrupted);
  if (res.interrupted) { return res; }

  if (!g.isOne() && cmp(g, Mp) != 0) {
    res.foundFactor = true;
    res.gcdValue = g;
    printf("  splitting the gcd (%zu bits) into prime factors ... ", g.bits());
    fflush(stdout);
    res.factors = splitFactors(g, cfg.exponent);
    printf("%zu found\n", res.factors.size());
  }
  return res;
}

// ---------------------------------------------------------------------------
// P+1 stage 2
//
// Continues from the stage-1 residue y1 = V_E(V_1,1). See Gpu.h
// (Gpu::pp1Stage2) for the pairing derivation; this is the driver around it --
// resume, checkpoint, progress, gcd, factor splitting -- structured like
// runPM1Stage2 but simpler: no B2-extension search (out of scope for this
// version, see Pp1Stage2Save.h), so the checkpoint is deleted on clean finish
// exactly like P-1's stage-2 checkpoint did before B2 extension existed --
// nothing will ever read a completed record as a seed, so keeping one around
// would just be a large unused file.
// ---------------------------------------------------------------------------
PP1Stage2Result runPP1Stage2(Gpu& gpu, const Config& cfg, const Words& y1,
                             Pp1Start start, u64 b1, u64 b2, u32 d, u32 w,
                             bool showProgress) {
  PP1Stage2Result res;
  res.b1 = b1;
  res.b2 = b2;
  res.d = d;
  res.w = w;
  res.seed = start.id();

  // Binds every accumulator this run reads or writes to the stage-1 residue
  // it belongs to -- P+1's analogue of P-1's xRes64.
  const u64 yRes64 = residue(y1);

  const string savePath = cfg.checkpointFile.empty()
      ? defaultPp1Stage2Path(cfg.exponent, b1, b2, start.id())
      : cfg.checkpointFile + ".pp1s2";

  Pp1Stage2State want;
  want.exponent = cfg.exponent;
  want.b1 = b1;
  want.b2 = b2;
  want.d = d;
  want.w = w;
  want.seed = start.id();
  want.yRes64 = yRes64;

  // ---- already complete? --------------------------------------------------
  // A prior run of this exact (b1, b2] for this seed left a COMPLETED record
  // behind rather than deleting it (see below) -- the same bargain P-1's own
  // stage 2 makes. Re-running the identical bounds, the ordinary case of
  // re-launching on the same exponent, should not re-walk it: check this
  // before building the plan at all, since sieving it is the wasted work.
  bool alreadyComplete = false;
  bool haveResume = false;
  Gpu::Stage2Pos resumePos;
  Words completedAcc;

  if (cfg.checkpoint) {
    Pp1Stage2State got;
    string err;
    if (loadPp1Stage2(savePath, want, got, err)) {
      if (got.complete) {
        alreadyComplete = true;
        completedAcc = std::move(got.acc);
        printf("  [%s stage 2 GPU] already complete in %s -- skipping to the gcd\n",
               ph(4), savePath.c_str());
      } else {
        resumePos.m = got.m;
        resumePos.jIdx = got.jIdx;
        resumePos.done = got.done;
        resumePos.acc = std::move(got.acc);
        resumePos.a = std::move(got.a);
        resumePos.s = std::move(got.s);
        haveResume = true;
      }
    } else if (err != "no stage-2 checkpoint") {
      printf("  [%s stage 2 GPU] ignoring checkpoint: %s\n", ph(4), err.c_str());
    }
  }
  res.reusedComplete = alreadyComplete;

  Words accWords;

  if (alreadyComplete) {
    accWords = std::move(completedAcc);
  } else {
    const Stage2Plan plan = buildStage2Plan(b1, b2, d, w);
    const u64 total = plan.nSlots + 2 * (plan.nBlocks() - 1);

    printf("  [%s stage 2 GPU] %s\n", ph(4), plan.describe().c_str());
    printf("  [%s stage 2 GPU] %llu primes in (%llu, %llu], %zu T-buffers\n", ph(4),
           (unsigned long long) plan.nPrimes, (unsigned long long) b1,
           (unsigned long long) b2, plan.jset.size());
    fflush(stdout);

    if (haveResume) {
      printf("  [%s stage 2 GPU] resuming from %s at %.2f%% (block %llu of %llu)\n",
             ph(4), savePath.c_str(), total ? 100.0 * double(resumePos.done) / double(total) : 0.0,
             (unsigned long long) (resumePos.m - plan.mFirst),
             (unsigned long long) plan.nBlocks());
    }

    Timer timer;
    Timer report;
    double lastAt = -1e9;
    const u64 doneAtStart = haveResume ? resumePos.done : 0;
    bool stopped = false;

    auto progress = [&](u64 done, u64 tot) {
      if (gInterrupted.load()) { stopped = true; return false; }
      if (!showProgress || !tot) { return true; }
      const double el = report.at();
      const u64 didHere = done > doneAtStart ? done - doneAtStart : 0;
      const double rate = el > 0 ? double(didHere) / el : 0;
      const double eta = rate > 0 ? double(tot - done) / rate : 0;
      if (stdoutIsTerminal()) {
        printf("\r  [%s stage 2 GPU] %6.2f%%  %llu/%llu muls (%.0f/s)  elapsed %s  ETA %s   ",
               ph(4), 100.0 * double(done) / double(tot),
               (unsigned long long) done, (unsigned long long) tot,
               rate, fmtDuration(el).c_str(), fmtDuration(eta).c_str());
        fflush(stdout);
      } else if (el - lastAt >= 60.0 || done == tot) {
        lastAt = el;
        printf("  [%s stage 2 GPU] %6.2f%%  %llu/%llu muls (%.0f/s)  elapsed %s  ETA %s\n",
               ph(4), 100.0 * double(done) / double(tot),
               (unsigned long long) done, (unsigned long long) tot,
               rate, fmtDuration(el).c_str(), fmtDuration(eta).c_str());
        fflush(stdout);
      }
      return true;
    };

    u32 saveEvery = 0;
    if (cfg.checkpoint && cfg.checkpointSeconds) {
      const double usPerMul = 2500.0;
      saveEvery = u32(std::clamp(double(cfg.checkpointSeconds) * 1e6 / usPerMul, 200.0, 200000.0));
    }

    auto save = [&](const Gpu::Stage2Pos& p) {
      Pp1Stage2State st = want;
      st.m = p.m;
      st.jIdx = p.jIdx;
      st.done = p.done;
      st.acc = p.acc;
      st.a = p.a;
      st.s = p.s;
      string err;
      if (!savePp1Stage2(savePath, st, err)) {
        printf("\n  WARNING: could not write P+1 stage-2 checkpoint: %s\n", err.c_str());
      }
    };

    Gpu::Stage2Pos stoppedAt;
    accWords =
        gpu.pp1Stage2(y1, plan, cfg.reportEvery ? cfg.reportEvery : 1000, progress,
                      nullptr, false, haveResume ? &resumePos : nullptr,
                      saveEvery, cfg.checkpoint ? save : std::function<void(const Gpu::Stage2Pos&)>{},
                      &stoppedAt);

    if (showProgress && stdoutIsTerminal()) { printf("\r%*s\r", 110, ""); }
    res.stage2Secs = timer.at();
    res.muls = total;

    if (stopped) {
      if (cfg.checkpoint && !stoppedAt.acc.empty()) { save(stoppedAt); }
      res.interrupted = true;
      printf("  [%s stage 2 GPU] interrupted; progress saved to %s\n", ph(4), savePath.c_str());
      return res;
    }

    // The walk is over: replace the resume point with a COMPLETED record --
    // the same bargain P-1's own stage 2 makes -- so a later identical run
    // skips straight to the gcd above instead of re-walking (b1, b2].
    if (cfg.checkpoint) {
      Pp1Stage2State st = want;
      st.complete = true;
      st.m = plan.mLast + 1;
      st.jIdx = 0;
      st.done = total;
      st.acc = accWords;
      string err;
      if (!savePp1Stage2(savePath, st, err)) {
        printf("  WARNING: could not record P+1 stage 2 as complete: %s\n", err.c_str());
      }
    }

    res.accRes64 = residue(accWords);
    printf("  [%s stage 2 GPU] done in %s (%.0f us/mul), acc res64 %016llx\n", ph(4),
           fmtDuration(res.stage2Secs).c_str(),
           total ? res.stage2Secs * 1e6 / double(total) : 0.0,
           (unsigned long long) res.accRes64);
  }

  res.accRes64 = residue(accWords);

  const Nat Mp = mersenne(cfg.exponent);
  Nat acc = fromWords(accWords);
  if (acc.isZero()) { acc = Mp; }

  Nat g = gcdWithProgress(acc, Mp, 5, "acc", showProgress, res.gcdSecs, res.interrupted);
  if (res.interrupted) { return res; }

  if (!g.isOne() && cmp(g, Mp) != 0) {
    res.foundFactor = true;
    res.gcdValue = g;
    printf("  splitting the gcd (%zu bits) into prime factors ... ", g.bits());
    fflush(stdout);
    Timer splitTimer;
    res.factors = splitFactors(g, cfg.exponent);
    printf("%zu found in %s\n", res.factors.size(), fmtDuration(splitTimer.at()).c_str());
  }

  // The checkpoint is deliberately NOT deleted here -- it now holds a
  // completed record, which is what lets a later identical run skip the walk
  // entirely instead of always starting stage 2 over from scratch.
  return res;
}

// ---------------------------------------------------------------------------
// P-1 stage 2
//
// Continues from the stage-1 residue x. See Stage2Plan.h for why one multiply
// covers two primes; this is the driver around it -- resume, checkpoint,
// progress, gcd, factor splitting.
// ---------------------------------------------------------------------------
// `stage1Residue`, not `residue`: there is a free function residue(Words) that a
// parameter of that name silently shadows, and the accumulator's res64 is
// computed with it below.
PM1Stage2Result runPM1Stage2(Gpu& gpu, const Config& cfg, const Words& stage1Residue,
                             u64 b1, u64 b2, u32 d, u32 w, bool showProgress,
                             const std::atomic<bool>* abortIf) {
  PM1Stage2Result res;
  res.b1 = b1;
  res.b2 = b2;
  res.d = d;
  res.w = w;

  // Binds every accumulator this run reads or writes to the stage-1 residue it
  // belongs to. Cheap, and it makes an extension across a changed E impossible.
  const u64 xRes64 = residue(stage1Residue);

  const string savePath = cfg.checkpointFile.empty()
      ? defaultStage2Path(cfg.exponent, b1, b2)
      : cfg.checkpointFile + ".s2";

  // ---- B2 extension -------------------------------------------------------
  // A COMPLETED stage 2 for this (exponent, b1) at a smaller B2 is exactly the
  // seed this run wants: the accumulator is a product over slots, so walking
  // (thatB2, b2] on top of it yields the product for the whole of (b1, b2].
  // The saving is exactly thatB2/b2 of the primes.
  //
  // The earlier run's pairing shape does NOT have to match this one's -- each
  // range gets its own plan and the product does not care -- so (D, w) are read
  // out for reporting rather than required to agree.
  //
  // Skipped when checkpoint_file is set: the search is over the default
  // pm1_<e>_s2_<b1>_<b2>.save naming, and an override means the user is placing
  // files deliberately.
  u64 lo = b1;
  Words seedAcc;
  string supersededPath;
  bool alreadyComplete = false;

  if (cfg.checkpoint && cfg.extend && cfg.checkpointFile.empty()) {
    for (u64 cand : findStage2Saves(cfg.exponent, b1)) {
      const string path = defaultStage2Path(cfg.exponent, b1, cand);
      Stage2State got;
      string err;
      const bool usable = loadCompletedStage2(path, cfg.exponent, b1, xRes64, got, err);

      if (cand > b2) {
        // Covers more than was asked for. Using it would mean reporting a B2
        // this run did not request, so it is left alone and merely mentioned.
        if (usable) {
          printf("  [%s stage 2 GPU] note: %s already covers B2=%llu -- raise b2 to reuse it\n",
                 ph(4), path.c_str(), (unsigned long long) cand);
        }
        continue;
      }
      if (!usable) {
        // A partial checkpoint for this run's own B2 is the normal
        // interrupted-run case, handled by the resume path below; saying
        // "ignoring" about it would be misleading.
        if (cand != b2) {
          printf("  [%s stage 2 GPU] not usable as a seed: %s (%s)\n", ph(4), path.c_str(), err.c_str());
        }
        continue;
      }

      if (cand == b2) {
        alreadyComplete = true;
        seedAcc = std::move(got.acc);
        printf("  [%s stage 2 GPU] already complete in %s -- skipping to the gcd\n",
               ph(4), path.c_str());
      } else {
        lo = cand;
        seedAcc = std::move(got.acc);
        supersededPath = path;
        res.fromB2 = cand;
      }
      break;                              // findStage2Saves returns largest first
    }
  }
  res.reusedComplete = alreadyComplete;

  Words accWords;

  if (alreadyComplete) {
    accWords = std::move(seedAcc);
  } else {
    const Stage2Plan plan = buildStage2Plan(lo, b2, d, w);
    const u64 total = plan.nSlots + 2 * (plan.nBlocks() - 1);

    if (res.fromB2) {
      // No "vs from scratch" prime count here: getting one means building the
      // full (b1, b2] plan, which sieves and greedily matches the entire range.
      // That is seconds of work at realistic bounds, for a log line.
      printf("  [%s stage 2 GPU] extending the completed B2=%llu result to B2=%llu"
             " -- only (%llu, %llu] is left to walk\n", ph(4),
             (unsigned long long) res.fromB2, (unsigned long long) b2,
             (unsigned long long) res.fromB2, (unsigned long long) b2);
    }
    printf("  [%s stage 2 GPU] %s\n", ph(4), plan.describe().c_str());
    printf("  [%s stage 2 GPU] %llu primes in (%llu, %llu], %zu T-buffers\n", ph(4),
           (unsigned long long) plan.nPrimes, (unsigned long long) plan.b1,
           (unsigned long long) plan.b2, plan.jset.size());
    fflush(stdout);

    // ---- resume ------------------------------------------------------------
    Stage2State want;
    want.exponent = cfg.exponent;
    want.b1 = b1;
    want.b2 = b2;
    want.d = d;
    want.w = w;
    want.fromB2 = res.fromB2;
    want.xRes64 = xRes64;

    Gpu::Stage2Pos resumePos;
    bool haveResume = false;
    if (cfg.checkpoint) {
      Stage2State got;
      string err;
      if (loadStage2(savePath, want, got, err)) {
        resumePos.m = got.m;
        resumePos.jIdx = got.jIdx;
        resumePos.done = got.done;
        resumePos.acc = std::move(got.acc);
        resumePos.a = std::move(got.a);
        resumePos.s = std::move(got.s);
        haveResume = true;
        printf("  [%s stage 2 GPU] resuming from %s at %.2f%% (block %llu of %llu)\n",
               ph(4), savePath.c_str(), total ? 100.0 * double(resumePos.done) / double(total) : 0.0,
               (unsigned long long) (resumePos.m - plan.mFirst),
               (unsigned long long) plan.nBlocks());
      } else if (err != "no stage-2 checkpoint") {
        printf("  [%s stage 2 GPU] ignoring checkpoint: %s\n", ph(4), err.c_str());
      }
    }

    // ---- the walk ----------------------------------------------------------
    Timer timer;
    Timer report;
    double lastAt = -1e9;
    // Squarings already done before this run started must not count against this
    // run's elapsed time, or the rate and ETA come out badly wrong on a resume.
    const u64 doneAtStart = haveResume ? resumePos.done : 0;
    bool stopped = false;

    auto progress = [&](u64 done, u64 tot) {
      if (gInterrupted.load()) { stopped = true; return false; }
      if (!showProgress || !tot) { return true; }
      const double el = report.at();
      const u64 didHere = done > doneAtStart ? done - doneAtStart : 0;
      const double rate = el > 0 ? double(didHere) / el : 0;
      const double eta = rate > 0 ? double(tot - done) / rate : 0;
      if (stdoutIsTerminal()) {
        printf("\r  [%s stage 2 GPU] %6.2f%%  %llu/%llu muls (%.0f/s)  elapsed %s  ETA %s   ",
               ph(4), 100.0 * double(done) / double(tot),
               (unsigned long long) done, (unsigned long long) tot,
               rate, fmtDuration(el).c_str(), fmtDuration(eta).c_str());
        fflush(stdout);
      } else if (el - lastAt >= 60.0 || done == tot) {
        lastAt = el;
        printf("  [%s stage 2 GPU] %6.2f%%  %llu/%llu muls (%.0f/s)  elapsed %s  ETA %s\n",
               ph(4), 100.0 * double(done) / double(tot),
               (unsigned long long) done, (unsigned long long) tot,
               rate, fmtDuration(el).c_str(), fmtDuration(eta).c_str());
        fflush(stdout);
      }
      return true;
    };

    // Checkpoint on a multiply count rather than a timer, because the saver has
    // to run at a slot boundary. Convert the configured interval using the rate
    // measured in M6b, then clamp so it is neither constant I/O nor never.
    u32 saveEvery = 0;
    if (cfg.checkpoint && cfg.checkpointSeconds) {
      const double usPerMul = 2500.0;
      saveEvery = u32(std::clamp(double(cfg.checkpointSeconds) * 1e6 / usPerMul, 200.0, 200000.0));
    }

    auto save = [&](const Gpu::Stage2Pos& p) {
      Stage2State st = want;
      st.m = p.m;
      st.jIdx = p.jIdx;
      st.done = p.done;
      st.acc = p.acc;
      st.a = p.a;
      st.s = p.s;
      string err;
      if (!saveStage2(savePath, st, err)) {
        printf("\n  WARNING: could not write stage-2 checkpoint: %s\n", err.c_str());
      }
    };

    Gpu::Stage2Pos stoppedAt;
    accWords =
        gpu.stage2(stage1Residue, plan, cfg.reportEvery ? cfg.reportEvery : 1000, progress,
                   nullptr, false, haveResume ? &resumePos : nullptr,
                   saveEvery, cfg.checkpoint ? save : std::function<void(const Gpu::Stage2Pos&)>{},
                   &stoppedAt,
                   // Already folded into the saved accumulator on a resume.
                   (!haveResume && !seedAcc.empty()) ? &seedAcc : nullptr,
                   abortIf);

    if (showProgress && stdoutIsTerminal()) { printf("\r%*s\r", 110, ""); }
    res.stage2Secs = timer.at();
    res.muls = total;

    if (stopped) {
      // Save exactly where the walk stopped, so a resume neither repeats nor
      // skips a slot.
      if (cfg.checkpoint && !stoppedAt.acc.empty()) { save(stoppedAt); }
      res.interrupted = true;
      printf("  [%s stage 2 GPU] interrupted; progress saved to %s\n", ph(4), savePath.c_str());
      return res;
    }

    // Checked AFTER `stopped`: if the user interrupted as well, the interrupt
    // is the answer the caller has to act on. Reaching here means the
    // overlapped stage-1 gcd raised the flag -- it already has a factor, so
    // this whole walk is moot and the job ends on stage 1's result.
    //
    // stoppedAt is empty when the abort landed inside the T_j table build,
    // before any slot was walked: nothing to save, nothing to resume, and the
    // best possible case for this flag.
    if (abortIf && abortIf->load(std::memory_order_relaxed)) {
      res.abandoned = true;
      if (stoppedAt.acc.empty()) {
        printf("  [%s stage 2 GPU] stage 1's gcd found a factor -- abandoned during\n"
               "      table setup, before any of the walk ran\n", ph(4));
      } else {
        const double pct = total ? 100.0 * double(stoppedAt.done) / double(total) : 0.0;
        printf("  [%s stage 2 GPU] stage 1's gcd found a factor -- abandoned at %.1f%%\n"
               "      (%llu of %llu muls); the rest of the walk is skipped\n",
               ph(4), pct, (unsigned long long) stoppedAt.done,
               (unsigned long long) total);
        // Real work, and a legitimate resume point if this exponent is ever
        // re-run with that factor declared known -- the accumulator is bound
        // to the stage-1 residue, so it cannot be misapplied elsewhere.
        if (cfg.checkpoint) { save(stoppedAt); }
      }
      return res;
    }

    // The walk is over: replace the resume point with a COMPLETED record. That
    // is what a later run extends to a larger B2, and writing it here rather
    // than after the gcd means a crash during the gcd -- which takes minutes --
    // costs nothing.
    if (cfg.checkpoint) {
      Stage2State st = want;
      st.complete = true;
      st.m = plan.mLast + 1;
      st.jIdx = 0;
      st.done = total;
      st.acc = accWords;
      // A and S are per-block state of a walk that is over: ~20 MB of nothing.
      string err;
      if (!saveStage2(savePath, st, err)) {
        printf("  WARNING: could not record stage 2 as complete: %s\n", err.c_str());
      } else if (!supersededPath.empty()) {
        // Only now, with the wider result safely on disk: a completed record for
        // a smaller B2 is strictly dominated by this one.
        std::error_code ec;
        if (filesystem::remove(supersededPath, ec)) {
          printf("  [%s stage 2 GPU] removed %s, superseded by B2=%llu\n", ph(4),
                 supersededPath.c_str(), (unsigned long long) b2);
        }
      }
    }

    // The accumulator's low bits are the cheapest way to tell two runs apart --
    // an interrupted-and-resumed walk must reproduce this exactly.
    res.accRes64 = residue(accWords);
    printf("  [%s stage 2 GPU] done in %s (%.0f us/mul), acc res64 %016llx\n", ph(4),
           fmtDuration(res.stage2Secs).c_str(),
           total ? res.stage2Secs * 1e6 / double(total) : 0.0,
           (unsigned long long) res.accRes64);
  }

  // ---- the gcd -----------------------------------------------------------
  res.accRes64 = residue(accWords);

  const Nat Mp = mersenne(cfg.exponent);
  Nat acc = fromWords(accWords);
  if (acc.isZero()) { acc = Mp; }

  Nat g = gcdWithProgress(acc, Mp, 5, "acc", showProgress, res.gcdSecs, res.interrupted);
  if (res.interrupted) { return res; }

  if (!g.isOne() && cmp(g, Mp) != 0) {
    res.foundFactor = true;
    res.gcdValue = g;
    printf("  splitting the gcd (%zu bits) into prime factors ... ", g.bits());
    fflush(stdout);
    Timer splitTimer;
    res.factors = splitFactors(g, cfg.exponent);
    printf("%zu found in %s\n", res.factors.size(), fmtDuration(splitTimer.at()).c_str());
  }

  // The checkpoint is deliberately NOT deleted here. It now holds a completed
  // record, and that accumulator is what lets a later run raise B2 without
  // re-walking (b1, b2] -- the same bargain stage 1 makes by keeping its own
  // completed residue.
  return res;
}

