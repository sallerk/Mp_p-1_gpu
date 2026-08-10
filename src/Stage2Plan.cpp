// Copyright (C) Mp_p-1_gpu

#include "Stage2Plan.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

// D must be EVEN. With D even every j coprime to D is odd, so m*D +- j is odd
// and no slot is wasted on an even (hence composite) candidate. An odd D would
// spend half its j values on numbers that cannot be prime.
//
// Within one radical -- 210, 420, 630 ... all being 2*3*5*7 times something --
// pairing quality is identical, because it depends on phi(D)/D. What a larger D
// buys there is fewer blocks, and so less of the 2-multiplies-per-block A/S
// recurrence overhead; what it costs is more T_j buffers. Bringing in 11 (2310)
// is the only step in this ladder that genuinely improves pairing:
// phi/D falls from 0.229 to 0.208.
const u32 STAGE2_D_CANDIDATES[] = {210, 420, 630, 840, 1050, 1260, 1680, 2310, 4620};
const size_t STAGE2_D_COUNT = sizeof(STAGE2_D_CANDIDATES) / sizeof(STAGE2_D_CANDIDATES[0]);

u32 stage2NumJ(u32 d, u32 w) {
  // phi(d) / 2 per half-D of window. Counting directly rather than via the
  // formula keeps this honest for any (d, w), including even w.
  u32 n = 0;
  const u64 jmax = u64(w) * d / 2;
  for (u64 j = 1; j < jmax; ++j) {
    if (std::gcd(j, u64(d)) == 1) { ++n; }
  }
  return n;
}

namespace {

std::vector<u32> basePrimes(u32 limit) {
  std::vector<bool> comp(size_t(limit) + 1, false);
  std::vector<u32> out;
  for (u32 i = 2; i <= limit; ++i) {
    if (comp[i]) { continue; }
    out.push_back(i);
    for (u64 k = u64(i) * i; k <= limit; k += i) { comp[size_t(k)] = true; }
  }
  return out;
}

// f(q) for every prime q in (lo, hi], ascending. Segmented over odds only, so
// the working set stays in cache even for a B2 in the hundreds of millions.
//
// A local sieve rather than the existing `Primes` class: that one answers
// "is n prime" by trial division and caps out at 10^10, whereas what is needed
// here is bulk enumeration of a range.
template <class F>
void forEachPrime(u64 lo, u64 hi, F f) {
  if (hi < 2 || lo >= hi) { return; }
  if (lo < 2) { f(2); }

  u32 root = u32(std::sqrt(double(hi))) + 2;
  while (u64(root) * root > hi) { --root; }
  const std::vector<u32> base = basePrimes(root + 1);

  u64 segLo = lo + 1;
  if ((segLo & 1) == 0) { ++segLo; }
  if (segLo < 3) { segLo = 3; }

  const u64 SPAN = 1u << 21;   // numbers spanned per segment; half that many odds
  std::vector<bool> comp;

  for (; segLo <= hi; segLo += SPAN) {
    const u64 segHi = std::min(segLo + SPAN - 2, hi);
    if (segHi < segLo) { break; }
    const size_t n = size_t((segHi - segLo) / 2) + 1;
    comp.assign(n, false);

    for (u32 p : base) {
      if (p == 2) { continue; }
      u64 start = u64(p) * p;
      if (start > segHi) { break; }
      if (start < segLo) {
        start = ((segLo + p - 1) / p) * u64(p);
        if ((start & 1) == 0) { start += p; }   // keep it odd
      }
      for (u64 k = start; k <= segHi; k += 2ull * p) { comp[size_t((k - segLo) / 2)] = true; }
    }

    for (size_t i = 0; i < n; ++i) {
      if (!comp[i]) { f(segLo + 2ull * i); }
    }
  }
}

// A bitmap of the primes in (lo, hi], over odds only. Needed because the
// greedy matching below must ask "is 2mD-q also prime" millions of times, in no
// particular order.
struct PrimeBits {
  u64 lo = 0;                 // odd, <= every member
  std::vector<u64> bits;

  void build(u64 b1, u64 b2) {
    lo = (b1 + 1) | 1;
    if (b2 < lo) { return; }
    bits.assign(size_t(((b2 - lo) / 2 + 1 + 63) / 64), 0);
    forEachPrime(b1, b2, [&](u64 q) {
      if (q < lo) { return; }   // only q == 2, which cannot occur here
      const u64 i = (q - lo) >> 1;
      bits[size_t(i >> 6)] |= 1ull << (i & 63);
    });
  }

  bool has(u64 n) const {
    if (n < lo || (n & 1) == 0) { return false; }
    const u64 i = (n - lo) >> 1;
    if ((i >> 6) >= bits.size()) { return false; }
    return (bits[size_t(i >> 6)] >> (i & 63)) & 1;
  }

  template <class F> void forEach(F f) const {
    for (size_t k = 0; k < bits.size(); ++k) {
      u64 v = bits[k];
      while (v) {
        const u32 b = u32(std::countr_zero(v));
        v &= v - 1;
        f(lo + 2ull * (u64(k) * 64 + b));
      }
    }
  }
};

} // namespace

bool Stage2Plan::slot(u64 m, size_t jIdx) const {
  if (m < mFirst || m > mLast || jIdx >= jset.size()) { return false; }
  const u64 idx = (m - mFirst) * jset.size() + jIdx;
  return (bits[size_t(idx >> 6)] >> (idx & 63)) & 1;
}

double Stage2Plan::mulsPerPrime() const {
  return nPrimes ? double(muls()) / double(nPrimes) : 0.0;
}

double Stage2Plan::pairing() const {
  return nSlots ? double(nPrimes) / double(nSlots) : 0.0;
}

std::string Stage2Plan::describe() const {
  char buf[256];
  snprintf(buf, sizeof(buf),
           "D=%-5u w=%u  %4zu bufs  %7llu blocks  %9llu primes -> %9llu muls "
           "(%.3f/prime, %.2f primes/mul)",
           d, w, jset.size(), (unsigned long long) nBlocks(),
           (unsigned long long) nPrimes, (unsigned long long) muls(),
           mulsPerPrime(), pairing());
  return buf;
}

Stage2Plan buildStage2Plan(u64 b1, u64 b2, u32 d, u32 w) {
  if (b2 <= b1) { throw std::runtime_error("stage 2: B2 must exceed B1"); }
  if (d < 2 || (d & 1)) { throw std::runtime_error("stage 2: D must be even"); }
  if (w < 1) { throw std::runtime_error("stage 2: window must be at least 1"); }

  const u64 jmax = u64(w) * d / 2;
  // Every prime in range must exceed jmax, or j could reach 0 / wrap past the
  // table and the m*D +- j decomposition stops being well defined.
  if (jmax > b1) { throw std::runtime_error("stage 2: w*D/2 must not exceed B1"); }

  Stage2Plan p;
  p.b1 = b1;
  p.b2 = b2;
  p.d = d;
  p.w = w;

  std::vector<int> idxOf(size_t(jmax), -1);
  for (u64 j = 1; j < jmax; ++j) {
    if (std::gcd(j, u64(d)) == 1) {
      idxOf[size_t(j)] = int(p.jset.size());
      p.jset.push_back(u32(j));
    }
  }
  const u64 J = p.jset.size();

  p.mFirst = (b1 - jmax) / d;
  if (p.mFirst == 0) { p.mFirst = 1; }
  p.mLast = (b2 + jmax) / d;
  if (p.mLast < p.mFirst) { p.mLast = p.mFirst; }

  p.bits.assign(size_t((p.nBlocks() * J + 63) / 64), 0);

  PrimeBits primes;
  primes.build(b1, b2);

  auto isSet = [&](u64 idx) { return (p.bits[size_t(idx >> 6)] >> (idx & 63)) & 1; };

  // Greedy matching, primes ascending.
  //
  //   1. If any candidate slot is already open, take it -- that slot was opened
  //      by the only other prime it can cover, so this is a completed pair and
  //      costs no multiply at all.
  //   2. Otherwise open the slot whose partner 2mD-q is itself a prime in
  //      range, reserving a pairing for later.
  //   3. Otherwise open any candidate slot.
  //
  // Processing in ascending order is what makes step 1 correct: a slot covering
  // q is only ever open because its other member, which is smaller, opened it.
  primes.forEach([&](u64 q) {
    ++p.nPrimes;

    const u64 mLo = std::max<u64>(p.mFirst, (q - jmax) / d + 1);
    const u64 mHi = std::min<u64>(p.mLast, (q + jmax - 1) / d);

    u64 reuse = 0, reserve = 0, any = 0;
    bool haveReuse = false, haveReserve = false, haveAny = false;

    for (u64 m = mLo; m <= mHi; ++m) {
      const u64 md = m * d;
      if (md == q) { continue; }
      const u64 j = (q > md) ? (q - md) : (md - q);
      if (j >= jmax) { continue; }
      const int ji = idxOf[size_t(j)];
      if (ji < 0) { continue; }   // cannot happen: j == +-q (mod d), so gcd(j,d)=1

      const u64 idx = (m - p.mFirst) * J + u64(ji);
      if (isSet(idx)) { reuse = idx; haveReuse = true; break; }
      // Only a partner ABOVE q is worth reserving for. A prime below q has
      // already been placed and will never come back to claim this slot, so
      // reserving against it wastes the multiply -- and with w > 1 most
      // candidates are below q, which is enough to make a wider window score
      // worse than a narrow one.
      if (!haveReserve && md > q && primes.has(2 * md - q)) { reserve = idx; haveReserve = true; }
      if (!haveAny) { any = idx; haveAny = true; }
    }

    if (haveReuse) { ++p.nPaired; return; }

    const u64 idx = haveReserve ? reserve : any;
    if (!haveReserve && !haveAny) {
      throw std::runtime_error("stage 2: prime " + std::to_string(q) + " has no slot");
    }
    p.bits[size_t(idx >> 6)] |= 1ull << (idx & 63);
    ++p.nSlots;
  });

  return p;
}

Stage2Shape chooseStage2Shape(u64 budgetBytes, u64 residueBytes, u64 b1) {
  // Efficiency is measured on a fixed sample range rather than the job's real
  // (B1,B2): the ranking of shapes is stable, and sieving the real range once
  // per candidate would cost more than the choice is worth.
  const u64 sampleB1 = 1000000, sampleB2 = 4000000;

  Stage2Shape best;
  double bestCost = 1e300;

  for (size_t i = 0; i < STAGE2_D_COUNT; ++i) {
    const u32 d = STAGE2_D_CANDIDATES[i];
    for (u32 w : {1u, 3u, 5u, 7u, 9u}) {
      if (u64(w) * d / 2 > std::min(b1, sampleB1)) { break; }
      const u64 bufs = stage2NumJ(d, w);
      if (bufs * residueBytes > budgetBytes) { break; }
      const double cost = buildStage2Plan(sampleB1, sampleB2, d, w).mulsPerPrime();
      if (cost < bestCost) { bestCost = cost; best = {d, w}; }
    }
  }
  return best;
}

// ---------------------------------------------------------------------------

namespace {

// Full independent audit of a plan. With a window every prime has SEVERAL
// candidate slots, so this cannot assume where a prime went; it checks the
// properties the GPU engine actually depends on:
//
//   A. the j table is exactly the coprime residues below w*D/2
//   B. nPrimes matches an independent count
//   C. nPrimes == nSlots + nPaired (every prime opens a slot or reuses one)
//   D. every prime is covered by at least one OPEN slot -- nothing is missed
//   E. every OPEN slot covers at least one prime in range -- nothing is wasted
//
// D and E are the load-bearing ones, and E in particular exercises the
// index -> (m, j) decode the GPU engine will use: get that wrong and md +- j
// lands on arbitrary numbers, which are almost never prime.
//
// There is deliberately no "check (mD)^2 - j^2 is divisible by q" here. That
// is the algebraic identity (mD-j)(mD+j) and holds for any integers at all, so
// testing it would prove nothing. The claim that actually needs arithmetic --
// that the accumulated product is divisible by a real factor of M_p -- can only
// be tested against a real residue, and is covered by the GPU stage-2 test.
int checkPlan(u64 b1, u64 b2, u32 d, u32 w, bool checkDecode) {
  int fails = 0;
  const Stage2Plan p = buildStage2Plan(b1, b2, d, w);
  const u64 jmax = u64(w) * d / 2;
  const u64 J = p.jset.size();

  // A
  if (J != stage2NumJ(d, w)) {
    printf("  FAIL D=%u w=%u: %llu j-values, expected %u\n", d, w,
           (unsigned long long) J, stage2NumJ(d, w));
    ++fails;
  }
  for (size_t i = 0; i < p.jset.size(); ++i) {
    const u32 j = p.jset[i];
    if (j == 0 || j >= jmax || std::gcd(u64(j), u64(d)) != 1 ||
        (i && j <= p.jset[i - 1])) {
      printf("  FAIL D=%u w=%u: bad j value %u at %zu\n", d, w, j, i);
      ++fails;
      break;
    }
  }

  std::vector<int> idxOf(size_t(jmax), -1);
  for (size_t i = 0; i < p.jset.size(); ++i) { idxOf[p.jset[i]] = int(i); }

  PrimeBits primes;
  primes.build(b1, b2);

  // B and D
  u64 count = 0;
  int missFails = 0;
  primes.forEach([&](u64 q) {
    ++count;
    const u64 mLo = std::max<u64>(p.mFirst, (q - jmax) / d + 1);
    const u64 mHi = std::min<u64>(p.mLast, (q + jmax - 1) / d);
    for (u64 m = mLo; m <= mHi; ++m) {
      const u64 md = m * d;
      if (md == q) { continue; }
      const u64 j = (q > md) ? (q - md) : (md - q);
      if (j >= jmax) { continue; }
      const int ji = idxOf[size_t(j)];
      if (ji >= 0 && p.slot(m, size_t(ji))) { return; }
    }
    if (missFails++ < 3) {
      printf("  FAIL D=%u w=%u: prime %llu is in no open slot\n", d, w, (unsigned long long) q);
    }
  });
  fails += missFails;

  if (count != p.nPrimes) {
    printf("  FAIL D=%u w=%u: counted %llu primes, plan says %llu\n", d, w,
           (unsigned long long) count, (unsigned long long) p.nPrimes);
    ++fails;
  }

  // C
  if (p.nSlots + p.nPaired != p.nPrimes) {
    printf("  FAIL D=%u w=%u: %llu slots + %llu paired != %llu primes\n", d, w,
           (unsigned long long) p.nSlots, (unsigned long long) p.nPaired,
           (unsigned long long) p.nPrimes);
    ++fails;
  }

  // E and F
  u64 set = 0, wasted = 0;
  int divFails = 0;
  for (u64 i = 0; i < p.nBlocks() * J; ++i) {
    if (!((p.bits[size_t(i >> 6)] >> (i & 63)) & 1)) { continue; }
    ++set;
    const u64 m = p.mFirst + i / J;
    const u64 j = p.jset[size_t(i % J)];
    const u64 md = m * d;
    const bool loPrime = md > j && primes.has(md - j);
    const bool hiPrime = primes.has(md + j);
    if (!loPrime && !hiPrime) { ++wasted; }
    if (checkDecode) {
      // The decoded j must be a real table entry and coprime to D, and the
      // decoded m must lie in range. Cheap, and it catches an off-by-one in the
      // index packing that E might otherwise only show statistically.
      if (m < p.mFirst || m > p.mLast || j == 0 || j >= jmax ||
          std::gcd(j, u64(d)) != 1) {
        if (divFails++ < 3) {
          printf("  FAIL D=%u w=%u: slot index %llu decodes to invalid m=%llu j=%llu\n",
                 d, w, (unsigned long long) i, (unsigned long long) m, (unsigned long long) j);
        }
      }
    }
  }
  fails += divFails;

  if (wasted) {
    printf("  FAIL D=%u w=%u: %llu open slots cover no prime\n", d, w, (unsigned long long) wasted);
    ++fails;
  }
  if (set != p.nSlots) {
    printf("  FAIL D=%u w=%u: %llu bits set, nSlots says %llu\n", d, w,
           (unsigned long long) set, (unsigned long long) p.nSlots);
    ++fails;
  }
  if (p.pairing() < 1.0 || p.pairing() > 2.0) {
    printf("  FAIL D=%u w=%u: pairing %.3f outside [1,2]\n", d, w, p.pairing());
    ++fails;
  }

  printf("  %-4s %s\n", fails ? "FAIL" : "ok", p.describe().c_str());
  return fails;
}

} // namespace

int runStage2PlanTests() {
  printf("stage-2 plan self-test\n");
  int fails = 0;

  printf("\n  correctness audit, B1=100K B2=3M\n");
  for (u32 d : {210u, 420u, 2310u}) {
    for (u32 w : {1u, 3u, 5u}) { fails += checkPlan(100000, 3000000, d, w, true); }
  }

  // What the pairing window actually buys, at bounds a real job would use.
  // This table is where the bounds model's cost-per-prime comes from.
  printf("\n  pairing window, B1=1M B2=30M\n");
  for (u32 d : {210u, 420u, 1050u, 2310u}) {
    for (u32 w : {1u, 3u, 5u, 7u, 9u}) {
      if (stage2NumJ(d, w) > 600) { continue; }
      fails += checkPlan(1000000, 30000000, d, w, false);
    }
  }

  printf("\n  shape chosen by memory budget (18 MB per residue, B1=1M)\n");
  const u64 res = 18ull << 20;
  for (u64 gb : {1ull, 2ull, 4ull, 6ull, 8ull}) {
    const Stage2Shape s = chooseStage2Shape(gb << 30, res, 1000000);
    const u32 bufs = s.d ? stage2NumJ(s.d, s.w) : 0;
    printf("    %llu GB -> D=%u w=%u (%u buffers, %.1f GB)\n", (unsigned long long) gb,
           s.d, s.w, bufs, double(bufs) * double(res) / double(1u << 30));
  }

  printf("\n%s\n", fails ? "stage-2 plan: FAILED" : "stage-2 plan: all tests passed");
  return fails ? 1 : 0;
}
