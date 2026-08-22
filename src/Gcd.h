// Copyright (C) Mp_p-1_gpu
//
// gcd for the P-1/P+1 factor extraction step: gcd(x-1, 2^p-1) with p > 1e8.
//
// Four implementations, each the oracle for the one below it:
//
//   gcdEuclid  -- textbook, one division per step. Slow, but small enough to
//                 be read and believed. The oracle every other tier is
//                 checked against, directly or transitively.
//   gcdLehmer  -- does ~64 bits of reduction per O(n) pass using only the top
//                 limbs. Same asymptotics as Euclid but a much better
//                 constant; gcdHalf's own base case and fallback.
//   gcdHalf    -- subquadratic recursive half-GCD (Schoenhage). O(M(n) log n)
//                 instead of O(n^2); this file's own hand-rolled ceiling
//                 (BigInt.cpp's multiply tops out at Toom-4).
//   gcdGmp     -- delegates to GMP's mpz_gcd. Same algorithm family as
//                 gcdHalf underneath (Lehmer + subquadratic HGCD), but GMP's
//                 hand-tuned assembly multiply (including an FFT tier this
//                 file has none of) measured ~11x faster at production scale
//                 -- see gcdGmp's own comment in Gcd.cpp for the numbers and
//                 what is given up (mid-computation interruptibility) to get
//                 them. This is gcd()'s actual implementation below.
//
// gcdHalf is not dead code for having been superseded as the production
// path: it stays fully self-tested as the reversion path if gcdGmp is ever
// pulled, and gcdLehmer/gcdEuclid remain its own correctness oracle.

#pragma once

#include "BigInt.h"

#include <atomic>
#include <functional>

Nat gcdEuclid(Nat a, Nat b);
Nat gcdLehmer(Nat a, Nat b);

// Diagnostic only: nanoseconds spent inside mul() (matCompose/matApply) vs.
// divrem()/mod() (reduceDirect, gcdHalf's no-progress fallback) during the
// large-operand part of gcdHalf. Read with .load(), reset with .store(0).
extern std::atomic<u64> gMulNanos, gDivNanos;

// Subquadratic: recursive half-GCD (Schoenhage). Reduces the operands by half
// their length per level using only their leading bits, so the cost is
// O(M(n) log n) instead of O(n^2). No longer gcd()'s implementation (see the
// file comment above) but kept fully live as the reversion path.
Nat gcdHalf(Nat a, Nat b);

// GMP-backed. gcd()'s actual implementation -- see the file comment above
// and this function's own comment in Gcd.cpp.
Nat gcdGmp(const Nat& a, const Nat& b);

// Optional progress hook. Called frequently -- after each batch of large
// multiplications, not merely once per halving -- with:
//   bitsNow    top-level operand size right now
//   bitsStart  operand size when the gcd began
//   muls       large multiplications completed so far
// Reporting purely by bitsNow is misleading: the top-level size halves, so a
// bits-linear percentage jumps 0 -> 50 -> 75 with each step taking a quarter
// the time of the last. Callers should weight by work; see GCD_WORK_EXPONENT.
// thread_local, and that is load-bearing rather than tidiness. Since 1.6 the
// driver can have TWO gcds in flight at once (stage 1's, overlapped with stage
// 2, plus stage 2's own). A single global hook would mean one gcd invoking the
// other's callback -- a lambda capturing the other call's locals by reference,
// on the wrong thread. Per-thread, each gcd sees only its own hook, and a gcd
// started on a thread that never installed one simply reports nothing.
// gcdHalf's own worker threads are unaffected: progress is ticked on the
// recursion spine, not in the offloaded multiplications.
//
// Return false to abort: gcdHalf throws GcdAborted from the next call site it
// reaches, on the SAME thread (each thread's hook only ever sees that
// thread's own gcd, so the throw unwinds a normal single-thread call chain
// even with two gcds concurrently in flight). Before this existed, the gcd
// phases were the only part of a run Ctrl-C could not stop: a full gcd is
// minutes, the GPU phases already checked an equivalent hook every
// reportEvery iterations, and a request to stop mid-gcd was silently ignored
// until it finished on its own -- at which point the job looked like it had
// completed normally and (with a stage 2 configured) got removed from
// worktodo.txt without ever writing a result.
//
// gcdGmp (gcd()'s actual implementation, see the file comment above) has no
// equivalent: mpz_gcd is one opaque call with no progress callback, so
// nothing can throw GcdAborted out of it. PM1.cpp's gcdWithProgress checks
// gInterrupted itself immediately before calling gcd() instead, so a Ctrl-C
// that lands before the call still takes effect; one already in flight now
// runs to completion, capped at gcdGmp's own measured worst case (~12s at
// production scale) rather than gcdHalf's (~140s) -- see gcdGmp's comment in
// Gcd.cpp. This hook and GcdAborted remain fully live for anything that
// calls gcdHalf directly (its own self-tests do). See main.cpp's
// finishStage1Gcd for the overlapped background gcd's own join.
extern thread_local std::function<bool(size_t bitsNow, size_t bitsStart, u64 muls)> gGcdProgress;

// Thrown by gcdHalf (via the tick() and top-of-loop call sites) when
// gGcdProgress returns false. Caught by whoever installed the hook -- gcdHalf
// itself has no notion of "interrupted," only "the hook said stop." Callers
// that do not install a hook (differential tests, --selftest, etc.) never see
// this: gGcdProgress defaults to empty, and an empty std::function is never
// called at all, let alone returns false.
struct GcdAborted {};

// Measured cost exponent of gcdHalf (4x size -> ~6.9x time). Work done by the
// time the operands have shrunk from N to n is about 1 - (n/N)^GCD_WORK_EXPONENT,
// which is a far better progress estimate than the bit ratio: at the halfway
// point in bits, roughly 62% of the work is actually done, not 50%.
inline constexpr double GCD_WORK_EXPONENT = 1.39;

// Cap the worker threads used by the half-GCD's independent multiplications.
// 0 restores the default (all hardware threads).
void setGcdThreads(unsigned n);

inline Nat gcd(const Nat& a, const Nat& b) { return gcdGmp(a, b); }
