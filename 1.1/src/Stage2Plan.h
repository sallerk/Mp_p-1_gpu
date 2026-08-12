// Copyright (C) Mp_p-1_gpu
//
// Stage-2 plan: which multiplies to do, in what order.
//
// Stage 1 finds q = 2kp+1 only when k is entirely B1-powersmooth. Stage 2 also
// allows ONE prime factor of k in (B1, B2], which is by far the most common
// near-miss -- k's largest prime factor is typically much larger than its
// second-largest.
//
// With x the stage-1 residue, a q whose k is smooth apart from one prime
// q' in (B1,B2] satisfies ord_q(x) | q'. So we want
//
//   gcd( product over primes q' in (B1,B2] of (x^q' - 1),  M_p )
//
// MONTGOMERY PAIRING gets two primes per multiply. Pick a highly composite D.
// Every prime q' > D is uniquely q' = m*D +- j with 0 < j < D/2 and
// gcd(j,D) = 1, taking m = round(q'/D). Then with
//
//   A_m = x^((m*D)^2)      T_j = x^(j^2)
//   A_m - T_j = x^(j^2) * ( x^((m*D)^2 - j^2) - 1 )
//   (m*D)^2 - j^2 = (m*D - j) * (m*D + j)
//
// the single term (A_m - T_j) is divisible by BOTH m*D-j and m*D+j whenever
// those are prime. x^(j^2) is a unit, so it contributes nothing spurious.
//
// This file only decides the (m, j) slots -- it touches no GPU state, so it is
// exactly testable on its own.

#pragma once

#include "common.h"

#include <string>
#include <vector>

// THE PAIRING WINDOW.
//
// With j restricted to (0, D/2) every prime has exactly ONE possible slot: the
// partner 2mD-q is the unique member of q's class mod 2D within distance D. So
// there is no choice to make, and the pairing rate is whatever luck gives --
// measured 1.17 primes per multiply at D=2310, which matches the closed form
// 2p/(1-(1-p)^2) with p = P(partner is prime) exactly.
//
// Letting j run up to w*D/2 for odd w > 1 gives each prime about w candidate
// slots, which turns slot assignment into a matching problem with real freedom
// and raises the pairing rate. The cost is w times as many T_j buffers.
// w = 1 reproduces the no-choice scheme above.
struct Stage2Plan {
  u64 b1 = 0;
  u64 b2 = 0;
  u32 d = 0;
  u32 w = 1;                  // pairing window, odd; j < w*d/2

  // The j values, ascending: 0 < j < w*d/2 with gcd(j,d) = 1.
  // Size is about phi(d)*w/2, and equals the number of GPU residue buffers.
  std::vector<u32> jset;

  u64 mFirst = 0;             // inclusive block range, m = round(q/d)
  u64 mLast = 0;

  // One bit per (m, j) slot, index (m - mFirst) * jset.size() + jIdx.
  std::vector<u64> bits;

  u64 nPrimes = 0;            // primes in (b1, b2]
  u64 nSlots = 0;             // set bits == accumulator multiplies
  u64 nPaired = 0;            // primes that reused a slot another prime opened

  u64 nBlocks() const { return mLast - mFirst + 1; }
  size_t nJ() const { return jset.size(); }

  bool slot(u64 m, size_t jIdx) const;

  // Accumulator multiplies plus the A/S recurrence, which costs 2 per block.
  u64 muls() const { return nSlots + 2 * nBlocks(); }

  // Multiplies per prime -- this is the `factorP2` the bounds model needs, and
  // it must come from measurement here rather than from another program's
  // numbers.
  double mulsPerPrime() const;

  // Primes covered per accumulator multiply, in [1, 2]. 2 would be perfect
  // pairing.
  double pairing() const;

  std::string describe() const;
};

// Highly composite D values, ascending. Larger D pairs better and needs fewer
// blocks, but the T_j table costs phi(d)/2 full residues of GPU memory.
extern const u32 STAGE2_D_CANDIDATES[];
extern const size_t STAGE2_D_COUNT;

// The number of T_j buffers (d, w) needs: the count of j in (0, w*d/2) coprime
// to d. Exactly phi(d)*w/2 for odd w.
u32 stage2NumJ(u32 d, u32 w = 1);

// The (D, w) combination that minimises predicted multiplies while its T_j
// table fits in `budgetBytes`, given one residue costs `residueBytes`. Returns
// d = 0 if even the smallest combination does not fit.
struct Stage2Shape { u32 d = 0; u32 w = 1; };
Stage2Shape chooseStage2Shape(u64 budgetBytes, u64 residueBytes, u64 b1);

// Throws std::runtime_error if b2 <= b1, d is odd, or w*d/2 > b1.
Stage2Plan buildStage2Plan(u64 b1, u64 b2, u32 d, u32 w = 1);

// --selftest stage2plan. Returns 0 on success.
int runStage2PlanTests();
