// Copyright (C) Mp_p-1_gpu
//
// Choosing B1 and B2 automatically.
//
// A factor q of M_p has the form 2kp+1. Stage 1 finds it when k is B1-smooth;
// stage 2 also finds it when k is B1-smooth apart from a single prime in
// (B1, B2]. The probability of each is given by Dickman's rho function,
// integrated over the plausible sizes of q, and the bounds worth using are the
// pair minimising expected total work.
//
// The PROBABILITY half of the model is ported from gpuowl's pm1/pm1.cpp (Preda)
// and is checked against that tool's own output in --selftest=bounds.
//
// The COST half is deliberately not ported. gpuowl's constants describe
// gpuowl's stage 2 on its hardware; this is a different implementation with a
// different inner loop, so its costs are measured here and the bounds follow
// from those measurements.

#pragma once

#include "common.h"

#include <string>

// Everything the work model needs about how expensive this build actually is,
// in units of one PRP iteration (one full-size squaring).
struct CostModel {
  // A stage-1 iteration is a squaring with the multiply by 3 folded into the
  // carry: a PRP iteration plus a little.
  double p1IterCost = 1.12;

  // A stage-2 accumulator multiply, relative to a squaring. MEASURED at 1.56 by
  // --selftest=stage2: a multiply needs two forward transforms plus an inverse,
  // where the fused-carry squaring needs one of each.
  double s2MulCost = 1.56;

  // Accumulator multiplies per prime in (B1,B2], including the per-block A/S
  // recurrence. Comes from the (D, w) shape actually chosen for the job --
  // Stage2Plan::mulsPerPrime() -- because it depends on how much GPU memory the
  // T_j table gets.
  double mulsPerPrime = 0.76;

  // One gcd, in PRP iterations. This is fixed per exponent rather than
  // proportional to B1, so it pushes the optimum towards larger bounds. With a
  // stage 2 there are two of them.
  double gcdIters = 0;

  // How much worse than the minimum expected cost is still acceptable, as a
  // fraction. Among everything inside that band, the bounds with the highest
  // P(factor) win.
  //
  // This matters because the minimum is extremely FLAT. Measured for M81679223
  // at TF76: going from the argmin (B1=0.2M, B2=6M) to 4x/8x those bounds raises
  // expected cost by 2.6% while nearly DOUBLING P(factor), 2.36% -> 4.39%. Total
  // work is 2.2% of a PRP test either way. Taking the exact argmin of a curve
  // that flat is false precision -- the model's own inputs are not good to 2.6%
  // (mulsPerPrime alone was out by 30% on the first real run), so a difference
  // that small is below the noise floor, while the probability difference is
  // real and large.
  //
  // 0 restores the strict argmin.
  double tolerance = 0.05;
};

// One gcd in PRP-iteration units. Calibrated from this program's own timings:
// ~14 min for the half-GCD at p = 86.6e6 against 1.63 ms per squaring there.
// The gcd is O(n^1.39) and a squaring roughly O(n log n) ~ O(n^1.05), so the
// ratio grows as about n^0.34.
double gcdIterCost(double exponent);

struct Bounds {
  u64 b1 = 0;
  u64 b2 = 0;              // == b1 means no stage 2 is worth running
  double probStage1 = 0;   // P(found in stage 1)
  double probStage2 = 0;   // P(found in stage 2), disjoint from probStage1
  double workStage1 = 0;   // PRP iterations, including its gcd
  double workStage2 = 0;   // PRP iterations, including its gcd; 0 if b2 == b1
  double expectedCost = 0; // what the scan minimised

  double prob() const { return probStage1 + probStage2; }
  std::string describe(double exponent) const;
};

// P(stage 1 finds it) and P(stage 2 finds it), disjoint. b2 <= b1 gives a zero
// second component.
void pm1Prob(double exponent, u32 factoredTo, double b1, double b2,
             double& p1, double& p2);

// Probability that stage 1 alone finds a factor. Kept because it is the natural
// reduction of pm1Prob at b2 == b1, and is checked to agree with it.
double pm1Stage1Prob(double exponent, u32 factoredTo, double B1);

// Same shape as pm1Prob, but for P+1's different smoothness target (q+1
// rather than q-1) -- per SEED, conditioned on that seed already being in
// the "genuine P+1" residue class. See the comment above this function's
// definition in Bounds.cpp for why the two targets need different math, not
// just different constants.
void pp1Prob(double exponent, u32 factoredTo, double b1, double b2,
             double& p1, double& p2);

// The (B1, B2) minimising expected total work. `bias` is how much a factor is
// worth relative to a PRP "composite" result: 1.0 means equally valuable, 2.0
// twice. A non-zero fixedB1 or fixedB2 pins that bound. allowStage2 = false
// restricts the same scan, same cost function, to B2 == B1 -- which is how the
// value of stage 2 is measured, rather than by comparing against a model that
// costs work differently.
Bounds chooseBounds(double exponent, u32 factoredTo, double bias,
                    const CostModel& cost, u64 fixedB1 = 0, u64 fixedB2 = 0,
                    bool allowStage2 = true);

// The B1 minimising expected total work across nSeeds independent P+1 seed
// attempts, each paying its own full stage-1 ladder (and stage-2 walk, if
// sharedB2 > the candidate b1) and gcd. B2/pairing shape are NOT chosen
// here -- P+1 has no cost model of its own for stage 2 yet, so sharedB2
// (already chosen by chooseBounds for P-1) is fixed input, mirroring how
// chooseBounds itself accepts a fixedB2. Returns the existing Bounds struct
// for symmetry with chooseBounds's result; its .b2 is an ECHO of sharedB2,
// not something this function chose.
Bounds choosePP1B1(double exponent, u32 factoredTo, double bias,
                   const CostModel& cost, u32 nSeeds, u64 sharedB2,
                   u64 fixedB1 = 0);

// There was a chooseB1() here carrying gpuowl's MERGED cost model, in which
// stage-1 squarings are credited against the PRP test that follows because
// gpuowl runs P-1 inside the PRP chain. This program does not merge, so that
// credit does not apply and the two models minimise different things -- which is
// exactly what made an early version of the bounds tests fail. It is deleted
// rather than kept "for reference": for a stage-1-only choice use
// chooseBounds(..., allowStage2 = false), which shares this program's own cost
// model and is therefore comparable.

// Approximate count of primes in (b1, b2].
double nPrimesBetween(double b1, double b2);

// --selftest=bounds. Returns 0 on success.
int runBoundsTests();

// --bounds <exponent> [tf] [bias]: print the cost/probability surface around the
// automatic choice, so the choice can be judged rather than taken on trust.
void printBoundsSurface(double exponent, u32 factoredTo, double bias,
                        const CostModel& cost);
