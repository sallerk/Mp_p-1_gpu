// Copyright (C) Mp_p-1_gpu
//
// See Gcd.h.
//
// Lehmer's method in one paragraph: a full division of two n-limb numbers costs
// O(n) but only removes ~64 bits. The quotient sequence, though, is almost
// always determined by the *leading* limbs alone. So we run Euclid on the two
// top limbs, accumulating the 2x2 cofactor matrix [[A,B],[C,D]] as we go, and
// stop as soon as the next quotient becomes ambiguous. Applying that one matrix
// to the full-length operands then performs many Euclid steps for a single O(n)
// pass. When the leading limbs cannot decide even one step (the matrix comes
// back as the identity), we fall back to a single real division so that
// progress is always made and the loop cannot stall.

#include "Gcd.h"
#include "Parallel.h"
#include "timeutil.h"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdlib>
#include <functional>
#include <utility>
#include <vector>

using namespace std;

std::atomic<u64> gMulNanos{0};
std::atomic<u64> gDivNanos{0};

namespace {
u64 nanosSince(const Timer& t) {
  return u64(t.at() * 1e9);
}
} // namespace

// Instrumentation: if the leading-limb loop never resolves a step, gcdLehmer
// silently degenerates into plain Euclid and every test still passes while the
// optimisation does nothing. These counters make that visible.
std::function<bool(size_t, size_t, u64)> gGcdProgress;

u64 gLehmerFastSteps = 0;
u64 gLehmerFallbackDivisions = 0;

// Thin alias: the threading budget itself is shared with BigInt.cpp's Toom
// tiers (see Parallel.h) since mulToom3/mulToom4 calls can nest inside an
// already-parallel matCompose batch.
void setGcdThreads(unsigned n) { setParallelThreads(n); }

Nat gcdEuclid(Nat a, Nat b) {
  while (!b.isZero()) {
    Nat r = mod(a, b);
    a = std::move(b);
    b = std::move(r);
  }
  return a;
}

namespace {

// The 64 bits of x starting at bit position s (zero above the top of x).
u64 bitsAt(const Nat& x, size_t s) {
  const size_t limb = s / 64, off = s % 64;
  const u64 lo = x[limb] >> off;
  const u64 hi = off ? (x[limb + 1] << (64 - off)) : 0;
  return lo | hi;
}

// A*a + B*b, where the caller knows the result is non-negative. Lehmer's
// cofactors alternate in sign, so at most one of A,B is negative and the
// combination reduces to a single add or subtract of two products.
Nat combine(i64 A, const Nat& a, i64 B, const Nat& b) {
  auto mag = [](i64 v) { return Nat(u64(v < 0 ? -v : v)); };
  if (A >= 0 && B >= 0) { return add(mul(a, mag(A)), mul(b, mag(B))); }
  if (A >= 0)           { return sub(mul(a, mag(A)), mul(b, mag(B))); }
  return sub(mul(b, mag(B)), mul(a, mag(A)));
}

// Cofactors are kept under 2^30 and the leading words under 2^32, so every
// intermediate product stays inside i64. That caps progress at ~30 bits per
// pass instead of ~60, which is a fair trade for not having to reason about
// overflow in the inner loop.
const i64 COF_LIMIT = i64(1) << 30;

} // namespace

Nat gcdLehmer(Nat a, Nat b) {
  if (cmp(a, b) < 0) { std::swap(a, b); }

  while (!b.isZero()) {
    // Once the operands fit in a single limb, finish with plain Euclid on u64.
    if (a.size() <= 1) {
      u64 x = a[0], y = b[0];
      while (y) { u64 t = x % y; x = y; y = t; }
      return Nat(x);
    }

    // Leading 32 bits of a, and the bits of b at the SAME absolute positions.
    // Aligning on bit position rather than limb index is essential: reading
    // b's limb at a's top index yields zero whenever b is shorter, which
    // silently disables the whole method.
    const size_t abits = a.bits();
    const size_t s = abits > 32 ? abits - 32 : 0;
    i64 ah = i64(bitsAt(a, s));
    i64 bh = i64(bitsAt(b, s));

    // Track a = A*a0 + B*b0, b = C*a0 + D*b0 across the simulated steps.
    i64 A = 1, B = 0, C = 0, D = 1;
    bool progressed = false;
    while (true) {
      // The next quotient is only trustworthy if it comes out the same at both
      // ends of the interval the true values are known to lie in.
      if (bh + C == 0 || bh + D == 0) { break; }
      const i64 q = (ah + A) / (bh + C);
      if (q <= 0 || q != (ah + B) / (bh + D)) { break; }

      // One Euclid step: (a,b) -> (b, a - q*b), carried on both the leading
      // words and the cofactors.
      i64 nC = A - q * C, nD = B - q * D;
      A = C; B = D; C = nC; D = nD;
      const i64 nbh = ah - q * bh;
      ah = bh; bh = nbh;
      progressed = true;

      if (bh <= 0) { break; }
      if (std::abs(A) > COF_LIMIT || std::abs(B) > COF_LIMIT ||
          std::abs(C) > COF_LIMIT || std::abs(D) > COF_LIMIT) { break; }
    }

    if (!progressed) {
      // Leading limbs could not resolve a single step: do one honest division.
      ++gLehmerFallbackDivisions;
      Nat r = mod(a, b);
      a = std::move(b);
      b = std::move(r);
      continue;
    }
    ++gLehmerFastSteps;

    // One O(n) pass performs every Euclid step the inner loop just validated.
    // Both results must be formed from the ORIGINAL a and b.
    Nat na = combine(A, a, B, b);
    Nat nb = combine(C, a, D, b);
    a = std::move(na);
    b = std::move(nb);

    if (cmp(a, b) < 0) { std::swap(a, b); }
  }
  return a;
}

// ===========================================================================
// Half-GCD
// ===========================================================================
//
// Representation. After an EVEN number of Euclid steps the cofactor matrix has
// a fixed sign pattern, so it can be stored with four non-negative entries:
//
//     a' = A*a - B*b        b' = D*b - C*a        with  A*D - B*C = 1
//
// Every routine below produces and consumes matrices in this form, and each
// composes two of them without any sign bookkeeping. Keeping the step count
// even is what buys that simplicity, so reduceDirect() is careful to stop on
// an even step.
//
// Correctness is not staked on the recursion constants. A matrix derived from
// truncated operands is always checked against the full-length values before
// being accepted (matApply returns false if either result would go negative,
// and the caller additionally requires a' > b'). If the check fails the matrix
// is discarded and the level falls back to direct Lehmer steps, which are
// already trusted. So GUARD_BITS and the thresholds affect only speed.

namespace {

// ---------------------------------------------------------------------------
// Parallelism.
//
// The half-GCD recursion is sequential -- each half needs the previous half's
// result -- but the time goes into large multiplications, and those come in
// independent batches: 8 in matCompose, 4 in matApply. Running those across
// cores is most of the available win. The dispatcher itself (runParallel) and
// its global thread budget live in Parallel.h/.cpp, shared with BigInt.cpp's
// Toom-tier sub-products -- see that header for why the budget must be shared
// rather than per-file.
// ---------------------------------------------------------------------------

// Progress bookkeeping. Updated on the recursion spine (only the individual
// multiplications are offloaded to workers), so no locking is needed.
u64 gMulCount = 0;
size_t gTopBits = 0, gStartBits = 0;

void tick(u64 muls) {
  gMulCount += muls;
  // On the recursion spine only (see the comment above gMulCount) -- this
  // never runs inside a runParallel worker task, so throwing here unwinds
  // through an ordinary single-thread call chain, not across threads.
  if (gGcdProgress && !gGcdProgress(gTopBits, gStartBits, gMulCount)) {
    throw GcdAborted{};
  }
}

struct Mat {
  Nat A, B, C, D;
};

Mat matIdentity() {
  Mat m;
  m.A = Nat(1); m.B = Nat(0); m.C = Nat(0); m.D = Nat(1);
  return m;
}

bool matIsIdentity(const Mat& m) {
  return m.A.isOne() && m.B.isZero() && m.C.isZero() && m.D.isOne();
}

// Timed wrappers, used only where operands are large enough that the split
// between multiplication and division time is worth measuring (see Gcd.h's
// gMulNanos/gDivNanos). Timer overhead is negligible next to a large-operand
// mul/divrem, so these are safe to leave permanently instrumented.
Nat timedMul(const Nat& a, const Nat& b) {
  Timer t;
  Nat r = mul(a, b);
  gMulNanos.fetch_add(nanosSince(t));
  return r;
}

void timedDivrem(const Nat& u, const Nat& v, Nat& q, Nat& r) {
  Timer t;
  divrem(u, v, q, r);
  gDivNanos.fetch_add(nanosSince(t));
}

Nat timedMod(const Nat& a, const Nat& b) {
  Timer t;
  Nat r = mod(a, b);
  gDivNanos.fetch_add(nanosSince(t));
  return r;
}

// m1 applied first, then m2. Derived by substituting
//   a1 = A1*a - B1*b, b1 = D1*b - C1*a
// into m2's formulas; every coefficient comes out non-negative.
Mat matCompose(const Mat& m1, const Mat& m2) {
  Nat t[8];
  const bool big = m1.A.size() + m2.A.size() >= PARALLEL_MIN_LIMBS;
  runParallel({
    [&]{ t[0] = timedMul(m2.A, m1.A); }, [&]{ t[1] = timedMul(m2.B, m1.C); },
    [&]{ t[2] = timedMul(m2.A, m1.B); }, [&]{ t[3] = timedMul(m2.B, m1.D); },
    [&]{ t[4] = timedMul(m2.D, m1.C); }, [&]{ t[5] = timedMul(m2.C, m1.A); },
    [&]{ t[6] = timedMul(m2.D, m1.D); }, [&]{ t[7] = timedMul(m2.C, m1.B); },
  }, big);
  Mat r;
  r.A = add(t[0], t[1]);
  r.B = add(t[2], t[3]);
  r.C = add(t[4], t[5]);
  r.D = add(t[6], t[7]);
  if (big) { tick(8); }
  return r;
}

// na = A*a - B*b, nb = D*b - C*a. False if either would be negative, which is
// how a matrix computed from truncated inputs is rejected.
bool matApply(const Mat& m, const Nat& a, const Nat& b, Nat& na, Nat& nb) {
  Nat p1, p2, q1, q2;
  const bool big = a.size() >= PARALLEL_MIN_LIMBS;
  runParallel({
    [&]{ p1 = timedMul(m.A, a); }, [&]{ p2 = timedMul(m.B, b); },
    [&]{ q1 = timedMul(m.D, b); }, [&]{ q2 = timedMul(m.C, a); },
  }, big);
  if (big) { tick(4); }
  if (cmp(p1, p2) < 0) { return false; }
  if (cmp(q1, q2) < 0) { return false; }
  na = sub(p1, p2);
  nb = sub(q1, q2);
  return true;
}

// Two Euclid steps with quotients q1 then q2, in even-parity form. A single
// step is odd parity and has no representation here, which is deliberate --
// keeping every matrix even is what removes the sign bookkeeping.
//   a -> b -> a - q1*b               (after step 1: a1=b,       b1=a-q1*b)
//   then                              (after step 2: a2=b1,      b2=a1-q2*b1)
//   a2 = a - q1*b                     = 1*a - q1*b
//   b2 = b - q2*(a - q1*b)            = (1 + q1*q2)*b - q2*a
Mat matStep2(const Nat& q1, const Nat& q2) {
  Mat m;
  m.A = Nat(1);
  m.B = q1;
  m.C = q2;
  m.D = add(Nat(1), mul(q1, q2));
  return m;
}

// Extra leading bits carried into a sub-problem so the quotient prefix it
// computes is still valid for the full-length operands.
const size_t GUARD_BITS = 64;

// Below this, recursion costs more than it saves and direct steps win.
const size_t DIRECT_BITS = 2048;

// Perform Euclid/Lehmer steps on (a,b) until bits(a) <= target, accumulating
// the (even-parity) cofactor matrix. Always correct; this is the fallback that
// every other path can retreat to.
Mat reduceDirect(Nat& a, Nat& b, size_t target) {
  Mat acc = matIdentity();
  while (!b.isZero() && a.bits() > target) {
    // Two divisions give an even-parity matrix directly. Using divrem rather
    // than the leading-word loop keeps this simple; it is only ever called on
    // operands that are already small.
    Nat q1, r1;
    timedDivrem(a, b, q1, r1);
    if (r1.isZero()) {
      // gcd found exactly; one more step would divide by zero.
      // Represent the single step and stop.
      Nat na = b, nb = r1;
      // odd number of steps: fold in a trivial second step with q2 = 0
      Mat m = matStep2(q1, Nat(0));
      Nat ta, tb;
      if (matApply(m, a, b, ta, tb)) { a = ta; b = tb; acc = matCompose(acc, m); }
      else { a = na; b = nb; }
      break;
    }
    Nat q2, r2;
    timedDivrem(b, r1, q2, r2);
    Mat m = matStep2(q1, q2);
    Nat na, nb;
    if (!matApply(m, a, b, na, nb)) { break; }   // should not happen
    a = std::move(na);
    b = std::move(nb);
    acc = matCompose(acc, m);
  }
  return acc;
}

// Reduce until bits(a) <= target, recursing on truncated operands.
Mat reduceRec(Nat& a, Nat& b, size_t target) {
  if (b.isZero() || a.bits() <= target) { return matIdentity(); }

  const size_t n = a.bits();
  if (n - target <= DIRECT_BITS) { return reduceDirect(a, b, target); }

  Mat acc = matIdentity();

  // Halve the remaining distance, twice.
  for (int half = 0; half < 2; ++half) {
    if (b.isZero() || a.bits() <= target) { break; }
    const size_t cur = a.bits();
    const size_t mid = (half == 0) ? target + (cur - target) / 2 : target;
    if (cur <= mid) { continue; }

    // To determine the quotients that shave (cur - mid) bits, the leading
    // 2*(cur-mid) + guard bits of the operands are enough.
    const size_t d = cur - mid;
    const size_t keep = 2 * d + GUARD_BITS;
    const size_t k = (cur > keep) ? cur - keep : 0;

    bool applied = false;
    if (k > 0) {
      Nat ah = shr(a, k), bh = shr(b, k);
      if (!bh.isZero() && ah.bits() > (mid > k ? mid - k : 0)) {
        Mat m = reduceRec(ah, bh, mid > k ? mid - k : 0);
        if (!matIsIdentity(m)) {
          Nat na, nb;
          // Accept only if it is valid for the FULL-length operands.
          if (matApply(m, a, b, na, nb) && cmp(na, nb) > 0) {
            a = std::move(na);
            b = std::move(nb);
            acc = matCompose(acc, m);
            applied = true;
          }
        }
      }
    }

    // Whether or not the truncated recursion helped, finish this half with
    // steps that are unconditionally correct.
    if (a.bits() > mid && !b.isZero()) {
      Mat m = (a.bits() - mid <= DIRECT_BITS) ? reduceDirect(a, b, mid)
                                              : reduceRec(a, b, mid);
      acc = matCompose(acc, m);
    }
    (void) applied;
  }

  return acc;
}

} // namespace

Nat gcdHalf(Nat a, Nat b) {
  if (cmp(a, b) < 0) { std::swap(a, b); }

  const size_t startBits = a.bits();
  gStartBits = startBits;
  gMulCount = 0;
  while (!b.isZero() && a.bits() > DIRECT_BITS) {
    const size_t n = a.bits();
    gTopBits = n;
    if (gGcdProgress && !gGcdProgress(n, startBits, gMulCount)) { throw GcdAborted{}; }
    // Drive a down to half its length, then let the loop repeat.
    reduceRec(a, b, n / 2);
    if (cmp(a, b) < 0) { std::swap(a, b); }
    if (a.bits() >= n) {
      // No progress at all: take one honest division so the loop cannot stall.
      Nat r = timedMod(a, b);
      a = std::move(b);
      b = std::move(r);
    }
  }
  return gcdLehmer(a, b);
}
