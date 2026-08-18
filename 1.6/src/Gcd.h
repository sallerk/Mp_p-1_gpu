// Copyright (C) Mp_p-1_gpu
//
// gcd for the P-1/P+1 factor extraction step: gcd(x-1, 2^p-1) with p > 1e8.
//
// Two implementations on purpose:
//
//   gcdEuclid  -- textbook, one division per step. Slow, but small enough to be
//                 read and believed. Used as the oracle in the tests.
//   gcdLehmer  -- does ~64 bits of reduction per O(n) pass using only the top
//                 limbs. Same asymptotics as Euclid but a much better constant;
//                 this is what runs today.
//
// A subquadratic half-GCD is the next step; it will be differentially tested
// against gcdLehmer before replacing it.

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
// O(M(n) log n) instead of O(n^2). This is the one that makes gcd(x-1, 2^p-1)
// practical at p > 1e8.
Nat gcdHalf(Nat a, Nat b);

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
extern thread_local std::function<void(size_t bitsNow, size_t bitsStart, u64 muls)> gGcdProgress;

// Measured cost exponent of gcdHalf (4x size -> ~6.9x time). Work done by the
// time the operands have shrunk from N to n is about 1 - (n/N)^GCD_WORK_EXPONENT,
// which is a far better progress estimate than the bit ratio: at the halfway
// point in bits, roughly 62% of the work is actually done, not 50%.
inline constexpr double GCD_WORK_EXPONENT = 1.39;

// Cap the worker threads used by the half-GCD's independent multiplications.
// 0 restores the default (all hardware threads).
void setGcdThreads(unsigned n);

inline Nat gcd(const Nat& a, const Nat& b) { return gcdHalf(a, b); }
