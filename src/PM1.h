// Copyright (C) Mp_p-1_gpu
//
// P-1 stage 1.
//
//   E = product of q^floor(log_q B1) over primes q <= B1
//   x = 3^E mod M_p
//   factor = gcd(x - 1, M_p)
//
// Any prime factor q of M_p whose k (in q = 2kp+1) is B1-smooth divides x-1,
// so it turns up in the gcd.

#pragma once

#include "BigInt.h"
#include "Config.h"
#include "common.h"

#include <atomic>

#include <string>
#include <vector>

class Gpu;

// One prime factor of M_p, as reported. `prime` and `dividesMp` are verified
// rather than assumed: a gcd can return several factors multiplied together,
// and reporting that product as "a factor" would be true but useless.
struct FoundFactor {
  Nat value;
  u64 k = 0;              // (value-1)/(2p), 0 if it does not fit in u64
  bool prime = false;     // Miller-Rabin
  bool dividesMp = false; // verified 2^p == 1 (mod value)
};

struct PM1Result {
  bool foundFactor = false;
  Nat factor;               // meaningful only if foundFactor
  bool interrupted = false;
  u64 b1Used = 0;
  double stage1Secs = 0;
  double gcdSecs = 0;
  u64 squarings = 0;
  Nat xMinusOne;            // stage-1 residue minus 1, always populated
  Nat gcdValue;             // the raw gcd, which may be a product of factors
  std::vector<FoundFactor> factors;   // gcdValue split into primes where possible
  Words residue;            // x itself, as stage 2 needs it on the GPU
};

// P-1 stage 2, run on the residue stage 1 left behind. See Stage2Plan.h for the
// pairing scheme; this is the driver around it -- checkpointing, progress, gcd
// and factor splitting.
struct PM1Stage2Result {
  bool foundFactor = false;
  bool interrupted = false;
  // Stopped early because the overlapped stage-1 gcd already had a factor,
  // which ends the job -- NOT the user interrupting, which does not. Kept
  // apart from `interrupted` because the two demand opposite things of the
  // caller: one says "resume this later", the other says "this work is moot".
  bool abandoned = false;
  u64 b1 = 0, b2 = 0;
  u32 d = 0, w = 0;
  double stage2Secs = 0;
  double gcdSecs = 0;
  u64 muls = 0;
  u64 accRes64 = 0;         // low 64 bits of the accumulated product
  u64 fromB2 = 0;           // seeded from a completed run to this B2; 0 = from B1
  bool reusedComplete = false;  // (b1, b2] was already finished; no walk at all
  Nat gcdValue;
  std::vector<FoundFactor> factors;
};

struct Stage2Plan;

// Builds its own plan rather than taking one, because the range it has to walk
// is not known until the checkpoint directory has been consulted: a completed
// run to a smaller B2 turns this into a walk over (thatB2, b2] seeded with its
// accumulator. Keeping that decision here puts all of stage 2's checkpoint
// policy in one place, next to the B1 extension it mirrors.
// abortIf, when set and raised, stops the walk (and the table build) at the
// next check and returns with `abandoned` set. Prime95 does the same thing by
// never starting stage 2: its stage-1 GCD does `goto bingo`, jumping past the
// whole stage-2 body (ecm.cpp), and running stage 2 anyway is an opt-in INI
// setting there. This program cannot skip stage 2 up front because its gcd is
// overlapped with it, so it stops it instead.
PM1Stage2Result runPM1Stage2(Gpu& gpu, const Config& cfg, const Words& stage1Residue,
                             u64 b1, u64 b2, u32 d, u32 w, bool showProgress,
                             const std::atomic<bool>* abortIf = nullptr);

// 3 when stage 1 runs alone, 5 when a stage 2 follows. Set by the driver before
// anything prints, so the "[2/5 stage 1 GPU]" tags count the real total.
extern u32 gPhaseTotal;

// Split a gcd result into the prime factors of M_p it contains. Every factor
// has the form 2kp+1, so candidates are enumerated by k rather than by trial
// division over all integers. Anything that cannot be split within maxK is
// returned whole, flagged with prime=false.
std::vector<FoundFactor> splitFactors(const Nat& g, u32 exponent, u64 maxK = 20000000);

// P+1 stage 1.
//
//   V = V_E(seed, 1) mod M_p        (Lucas sequence, Q = 1)
//   factor = gcd(V - 2, M_p)
//
// P+1 finds q when (seed^2 - 4) is a quadratic NON-residue mod q and q+1 is
// B1-smooth; when it IS a residue the computation degenerates into a P-1 run
// with base related to the seed. Which case applies is ~50/50 and cannot be
// known in advance, so several seeds are tried.
// P+1's starting point, as the rational num/den Prime95 uses rather than a
// bare integer. The Lucas sequence starts at V_1 = num/den mod M_p, which is
// just num * den^-1 -- see ratioModMersenne in BigInt.h.
//
// Peter Montgomery's choices, and Prime95's (ecm.cpp, "Convert run number
// into starting point"):
//
//   run 1   2/7     finds every factor q with q+1 divisible by 6
//   run 2   6/5     finds every factor q with q+1 divisible by 4
//   run 3   random  no special structure; q+1 divisible by 2 only
//
// which is why 2/7 is worth running before 6/5, and both before any random
// pair. An integer seed s -- what this program used before 1.9.4 -- is the
// degenerate case s/1, with none of that structure.
struct Pp1Start {
  u32 num = 2;
  u32 den = 7;

  // A single u32 identifying this start, for the save records, which have
  // always carried one u32 of "which computation is this". num <= 199 and
  // den <= 71, so the pair packs without collision and the save FORMAT is
  // unchanged -- only the value it holds.
  u32 id() const { return num * 256 + den; }

  // "2/7", exactly as it appears in results.txt.
  std::string label() const { return std::to_string(num) + "/" + std::to_string(den); }

  // The denominators Prime95 draws from for run 3 and above.
  static const u32 SMALL_PRIMES[16];

  // Prime95's table. For run 3 (and, in Prime95, anything above it) the
  // numerator and denominator are drawn with
  // rand(), having called srand(time(NULL)) -- so its start differs run to
  // run, and it has to write the pair into its save file and read it back to
  // resume (ecm.cpp's "Replace random start with random start from save
  // file"). Deriving the pair from (exponent, nthRun) instead gives the same
  // ranges and the same absence of structure, while making a run reproducible
  // and needing no save-format change at all: a resumed run recomputes the
  // start it was already using.
  static Pp1Start forRun(u32 nthRun, u32 exponent);
};

struct PP1Result {
  bool foundFactor = false;
  Nat gcdValue;
  std::vector<FoundFactor> factors;
  Pp1Start startUsed;
  bool interrupted = false;
  u64 b1Used = 0;
  double stage1Secs = 0;
  double gcdSecs = 0;
  u64 squarings = 0;
  Words residue;            // V_E(seed,1) itself, as stage 2 needs it on the GPU
};

PP1Result runPP1Stage1(Gpu& gpu, const Config& cfg, u64 b1, Pp1Start start,
                       bool showProgress);

// P+1 stage 2, run on the V_E residue stage 1 left behind. See Gpu.h
// (Gpu::pp1Stage2) for the pairing derivation; this is the driver around it --
// checkpointing, progress, gcd and factor splitting, structured the same way
// runPM1Stage2 is. An identical re-run (same exponent, b1, b2, seed) reuses a
// completed checkpoint and skips straight to the gcd, the same as P-1's own
// stage 2. Unlike P-1's driver, this does not search for a completed
// SMALLER-B2 checkpoint to extend into a wider walk: that B2-extension case is
// out of scope for this version (see Pp1Stage2Save.h) -- only the exact-match
// and interrupted-walk-resume cases are handled here.
struct PP1Stage2Result {
  bool foundFactor = false;
  bool interrupted = false;
  u64 b1 = 0, b2 = 0;
  u32 d = 0, w = 0;
  u32 seed = 0;
  double stage2Secs = 0;
  double gcdSecs = 0;
  u64 muls = 0;
  u64 accRes64 = 0;         // low 64 bits of the accumulated product
  bool reusedComplete = false;  // (b1, b2] was already finished; no walk at all
  Nat gcdValue;
  std::vector<FoundFactor> factors;
};

PP1Stage2Result runPP1Stage2(Gpu& gpu, const Config& cfg, const Words& y1,
                             Pp1Start start, u64 b1, u64 b2, u32 d, u32 w,
                             bool showProgress);

// R such that E(b1To) == E(b1From) * R -- the extra work an extension needs.
Nat stage1ExponentDelta(u64 b1From, u64 b1To, u32 exponent);

// Completed stage-1 checkpoints for this exponent with B1 < below, largest
// first.
std::vector<u64> findCompletedB1(u32 exponent, u64 below);

// Builds E for the given B1 and exponent. Exposed for testing.
Nat stage1Exponent(u64 b1, u32 exponent);

// doGcd == false stops after stage 1 and leaves the residue in xMinusOne. The
// self-tests use that: checking a KNOWN factor q only needs (x-1) mod q, which
// is a 2-limb division instead of a full multi-minute GCD.
PM1Result runPM1Stage1(Gpu& gpu, const Config& cfg, u64 b1, bool showProgress,
                       bool doGcd = true);

// The stage-1 gcd on its own, for callers that took doGcd == false in order to
// do something else first. Reads res.xMinusOne and writes res.foundFactor /
// factor / gcdValue / factors / gcdSecs; touches no GPU state, so it is safe
// to run on another thread while the GPU works (which is what the driver does
// to overlap it with stage 2).
void finishStage1Gcd(PM1Result& res, u32 exponent, bool showProgress,
                     bool announce = true);

// Human-readable duration, matching mersenne_tf: 45s, 12m34s, 1h23m, 3d04h.
std::string fmtDuration(double secs);
