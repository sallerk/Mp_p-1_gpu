// Copyright (C) Mp_p-1_gpu
//
// Minimal non-negative multiprecision integers, base 2^64, little-endian limbs.
//
// This exists because P-1/P+1 needs gcd(x-1, 2^p-1) with p > 1e8, and there is
// no GMP on the target toolchain (MSVC only). Only what the GCD needs is
// implemented: compare, add, sub, shift, multiply (schoolbook / Karatsuba /
// Toom-Cook-3, by size) and truncating division (Knuth algorithm D).
//
// Invariant: `w` never has trailing zero limbs, so w.size() is the exact limb
// count and zero is represented by an empty vector.

#pragma once

#include "common.h"

#include <string>
#include <vector>

struct Nat {
  std::vector<u64> w;

  Nat() = default;
  explicit Nat(u64 v) { if (v) { w.push_back(v); } }

  void norm() { while (!w.empty() && w.back() == 0) { w.pop_back(); } }
  bool isZero() const { return w.empty(); }
  size_t size() const { return w.size(); }
  u64 operator[](size_t i) const { return i < w.size() ? w[i] : 0; }

  // Index of the highest set bit plus one; 0 for zero.
  size_t bits() const;
  bool bit(size_t i) const;

  bool isOne() const { return w.size() == 1 && w[0] == 1; }

  std::string hex() const;
  std::string dec() const;   // decimal, which is how factors are conventionally reported
};

// -1, 0, +1
int cmp(const Nat& a, const Nat& b);
inline bool operator<(const Nat& a, const Nat& b)  { return cmp(a, b) < 0; }
inline bool operator==(const Nat& a, const Nat& b) { return cmp(a, b) == 0; }

Nat add(const Nat& a, const Nat& b);
Nat sub(const Nat& a, const Nat& b);          // requires a >= b
Nat shl(const Nat& a, size_t n);
Nat shr(const Nat& a, size_t n);

Nat mul(const Nat& a, const Nat& b);          // schoolbook / Karatsuba / Toom-3 by size
// Each tier exposed so tests can cross-check it against the one below it.
Nat mulSchoolbook(const Nat& a, const Nat& b);
Nat mulKaratsuba(const Nat& a, const Nat& b);
Nat mulToom3(const Nat& a, const Nat& b);

// u = q*v + r with r < v. v must be non-zero.
void divrem(const Nat& u, const Nat& v, Nat& q, Nat& r);
Nat mod(const Nat& u, const Nat& v);

// 2^p - 1
Nat mersenne(u32 p);

// Decimal string -> Nat. Returns false on any non-digit.
bool fromDecimal(const std::string& s, Nat& out);

// a mod m for a single-limb m, without allocating. Trial division over
// millions of candidates goes through this.
u64 modU64(const Nat& a, u64 m);

// Modular arithmetic, for splitting and verifying a gcd result. These operate
// on small operands (a gcd of M_p is a few hundred bits at most), so the naive
// mul-then-reduce is fine.
Nat mulMod(const Nat& a, const Nat& b, const Nat& m);
Nat powMod(Nat base, const Nat& e, const Nat& m);

// Miller-Rabin. Deterministic for n < 3.3e24 with the bases used; a probable
// prime above that.
bool isProbablePrime(const Nat& n);

// The GPU residue is a packed little-endian bit array of 32-bit words.
Nat fromWords(const Words& words);
