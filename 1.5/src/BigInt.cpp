// Copyright (C) Mp_p-1_gpu
//
// See BigInt.h. Nothing here is clever; it is meant to be obviously correct,
// with Karatsuba and Toom-Cook-3/4 as the only asymptotic concessions. The GCD
// on top is where the real work happens.

#include "BigInt.h"
#include "Parallel.h"

#include <cassert>
#include <cstring>
#include <algorithm>

#ifdef _MSC_VER
#include <intrin.h>
#endif

using namespace std;

namespace {

// ---------------------------------------------------------------------------
// 64-bit primitives. MSVC has no __uint128_t, so these wrap the intrinsics.
// ---------------------------------------------------------------------------

inline u64 mul64(u64 a, u64 b, u64* hi) {
#ifdef _MSC_VER
  return _umul128(a, b, hi);
#else
  unsigned __int128 t = (unsigned __int128) a * b;
  *hi = u64(t >> 64);
  return u64(t);
#endif
}

// out = a + b + carryIn; returns carry out (0 or 1)
inline unsigned char addc(unsigned char carryIn, u64 a, u64 b, u64* out) {
#ifdef _MSC_VER
  return _addcarry_u64(carryIn, a, b, out);
#else
  unsigned __int128 t = (unsigned __int128) a + b + carryIn;
  *out = u64(t);
  return (unsigned char) (t >> 64);
#endif
}

// out = a - b - borrowIn; returns borrow out (0 or 1)
inline unsigned char subb(unsigned char borrowIn, u64 a, u64 b, u64* out) {
#ifdef _MSC_VER
  return _subborrow_u64(borrowIn, a, b, out);
#else
  unsigned __int128 t = (unsigned __int128) a - b - borrowIn;
  *out = u64(t);
  return (unsigned char) ((t >> 64) & 1);
#endif
}

inline int clz64(u64 x) {
  assert(x);
#ifdef _MSC_VER
  unsigned long i;
  _BitScanReverse64(&i, x);
  return 63 - int(i);
#else
  return __builtin_clzll(x);
#endif
}

// (hi:lo) / d, requires hi < d so the quotient fits in 64 bits.
inline u64 div128by64(u64 hi, u64 lo, u64 d, u64* rem) {
  assert(hi < d);
#if defined(_MSC_VER) && defined(_M_X64)
  return _udiv128(hi, lo, d, rem);
#else
  unsigned __int128 n = ((unsigned __int128) hi << 64) | lo;
  *rem = u64(n % d);
  return u64(n / d);
#endif
}

const size_t KARATSUBA_LIMBS = 40;   // below this, schoolbook wins
const size_t TOOM3_LIMBS = 8000;     // below this, Karatsuba wins -- tuned empirically

// mulToom4 is measurably faster than mulToom3 in isolated single-multiply
// benchmarks (1.2-1.4x from 8000 up to 1,000,000 limbs), but measurably
// SLOWER in the real recursive gcdHalf workload it actually needs to serve
// (production-scale gcdHalf: 229.5s with Toom-3 as the top tier vs 250s+
// with Toom-4 active, reproduced twice, both with and without Toom-4's own
// internal parallelism). The isolated benchmark measures one giant multiply;
// the real workload is many multiplies of many sizes across gcdHalf's
// recursion, where Toom-4's higher per-call overhead (7 sub-products and 10
// SNat combinations vs Toom-3's 5 and 6) apparently doesn't amortize the
// same way. Root cause not fully identified -- set high enough that mul()
// never reaches mulToom4 in practice, rather than ship a regression; kept
// implemented and differentially tested since it's correct, just not yet
// proven to help the workload that matters here.
const size_t TOOM4_LIMBS = 10000000;

} // namespace

// ---------------------------------------------------------------------------
// Basics
// ---------------------------------------------------------------------------

size_t Nat::bits() const {
  if (w.empty()) { return 0; }
  return w.size() * 64 - clz64(w.back());
}

bool Nat::bit(size_t i) const {
  size_t limb = i / 64;
  return limb < w.size() && ((w[limb] >> (i % 64)) & 1);
}

string Nat::hex() const {
  if (w.empty()) { return "0"; }
  string s;
  char buf[32];
  snprintf(buf, sizeof(buf), "%llx", (unsigned long long) w.back());
  s += buf;
  for (size_t i = w.size() - 1; i-- > 0; ) {
    snprintf(buf, sizeof(buf), "%016llx", (unsigned long long) w[i]);
    s += buf;
  }
  return s;
}

string Nat::dec() const {
  if (w.empty()) { return "0"; }
  // Repeatedly divide by 10^19, the largest power of ten fitting in a limb.
  const u64 CHUNK = 10000000000000000000ull;
  Nat cur = *this;
  const Nat d(CHUNK);
  string out;
  while (!cur.isZero()) {
    Nat q, r;
    divrem(cur, d, q, r);
    const u64 part = r.isZero() ? 0 : r.w[0];
    char buf[32];
    if (q.isZero()) { snprintf(buf, sizeof(buf), "%llu", (unsigned long long) part); }
    else            { snprintf(buf, sizeof(buf), "%019llu", (unsigned long long) part); }
    out = string(buf) + out;
    cur = std::move(q);
  }
  return out;
}

int cmp(const Nat& a, const Nat& b) {
  if (a.w.size() != b.w.size()) { return a.w.size() < b.w.size() ? -1 : 1; }
  for (size_t i = a.w.size(); i-- > 0; ) {
    if (a.w[i] != b.w[i]) { return a.w[i] < b.w[i] ? -1 : 1; }
  }
  return 0;
}

Nat add(const Nat& a, const Nat& b) {
  const size_t n = max(a.size(), b.size());
  Nat r;
  r.w.resize(n + 1);
  unsigned char c = 0;
  for (size_t i = 0; i < n; ++i) { c = addc(c, a[i], b[i], &r.w[i]); }
  r.w[n] = c;
  r.norm();
  return r;
}

Nat sub(const Nat& a, const Nat& b) {
  assert(cmp(a, b) >= 0);
  Nat r;
  r.w.resize(a.size());
  unsigned char bw = 0;
  for (size_t i = 0; i < a.size(); ++i) { bw = subb(bw, a[i], b[i], &r.w[i]); }
  assert(bw == 0);
  r.norm();
  return r;
}

Nat shl(const Nat& a, size_t n) {
  if (a.isZero()) { return Nat{}; }
  const size_t limbs = n / 64, bits = n % 64;
  Nat r;
  r.w.assign(a.size() + limbs + 1, 0);
  if (bits == 0) {
    for (size_t i = 0; i < a.size(); ++i) { r.w[i + limbs] = a.w[i]; }
  } else {
    u64 carry = 0;
    for (size_t i = 0; i < a.size(); ++i) {
      r.w[i + limbs] = (a.w[i] << bits) | carry;
      carry = a.w[i] >> (64 - bits);
    }
    r.w[a.size() + limbs] = carry;
  }
  r.norm();
  return r;
}

Nat shr(const Nat& a, size_t n) {
  const size_t limbs = n / 64, bits = n % 64;
  if (limbs >= a.size()) { return Nat{}; }
  Nat r;
  r.w.assign(a.size() - limbs, 0);
  if (bits == 0) {
    for (size_t i = 0; i < r.w.size(); ++i) { r.w[i] = a.w[i + limbs]; }
  } else {
    for (size_t i = 0; i < r.w.size(); ++i) {
      u64 lo = a.w[i + limbs] >> bits;
      u64 hi = (i + limbs + 1 < a.size()) ? (a.w[i + limbs + 1] << (64 - bits)) : 0;
      r.w[i] = lo | hi;
    }
  }
  r.norm();
  return r;
}

// ---------------------------------------------------------------------------
// Multiplication
// ---------------------------------------------------------------------------

Nat mulSchoolbook(const Nat& a, const Nat& b) {
  if (a.isZero() || b.isZero()) { return Nat{}; }
  Nat r;
  r.w.assign(a.size() + b.size(), 0);
  for (size_t i = 0; i < a.size(); ++i) {
    u64 carry = 0;
    const u64 ai = a.w[i];
    for (size_t j = 0; j < b.size(); ++j) {
      u64 hi, lo = mul64(ai, b.w[j], &hi);
      // r[i+j] += lo + carry, propagating into hi
      unsigned char c = 0;
      c = addc(0, r.w[i + j], lo, &r.w[i + j]);
      hi += c;
      c = addc(0, r.w[i + j], carry, &r.w[i + j]);
      hi += c;
      carry = hi;
    }
    // Flush the final carry; it cannot overflow because r is wide enough.
    size_t k = i + b.size();
    while (carry) {
      unsigned char c = addc(0, r.w[k], carry, &r.w[k]);
      carry = c;
      ++k;
    }
  }
  r.norm();
  return r;
}

namespace {

// Split at `k` limbs: a = lo + hi * 2^(64k)
void split(const Nat& a, size_t k, Nat& lo, Nat& hi) {
  if (a.size() <= k) { lo = a; hi = Nat{}; return; }
  lo.w.assign(a.w.begin(), a.w.begin() + k);
  hi.w.assign(a.w.begin() + k, a.w.end());
  lo.norm();
  hi.norm();
}

// ---------------------------------------------------------------------------
// Toom-Cook-3 support: signed magnitudes.
//
// Evaluating A(x)=a2*x^2+a1*x+a0 at x=-1 and combining the five evaluation
// products (see mulToom3) genuinely needs negative intermediates -- unlike
// Gcd.cpp's combine(), which can assume at most one side is negative, either
// side can be negative here. Nat itself stays non-negative-only (its
// documented contract in BigInt.h); this wraps it with a sign, kept local to
// this file. The five true output coefficients are always non-negative
// (product of non-negative-coefficient polynomials), so mulToom3 asserts
// that once interpolation is done -- a cheap, direct check on the formulas.
// ---------------------------------------------------------------------------

struct SNat { bool neg = false; Nat mag; };

// The only SNat constructor: forces neg=false whenever mag is zero, so "-0"
// can never appear. Every function below returns through this -- sAdd's
// cancelling branch and sMul's zero-operand case would otherwise both be
// able to produce a negative zero, which is numerically harmless but trips
// mulToom3's non-negative-coefficient assert for no real reason.
SNat mk(bool neg, Nat mag) { return {!mag.isZero() && neg, std::move(mag)}; }

SNat sOf(const Nat& n) { return mk(false, n); }

SNat sAdd(const SNat& x, const SNat& y) {
  if (x.mag.isZero()) { return y; }
  if (y.mag.isZero()) { return x; }
  if (x.neg == y.neg) { return mk(x.neg, add(x.mag, y.mag)); }
  if (cmp(x.mag, y.mag) >= 0) { return mk(x.neg, sub(x.mag, y.mag)); }
  return mk(y.neg, sub(y.mag, x.mag));
}

SNat sSub(const SNat& x, const SNat& y) { return sAdd(x, mk(!y.neg, y.mag)); }

SNat sMul(const SNat& x, const SNat& y) { return mk(x.neg != y.neg, mul(x.mag, y.mag)); }

SNat sMulSmall(const SNat& x, u64 k) { return mk(x.neg, mul(x.mag, Nat(k))); }

SNat sShl(const SNat& x, size_t n) { return mk(x.neg, shl(x.mag, n)); }

// Exact division by 2^n; the caller guarantees no remainder.
SNat sShr(const SNat& x, size_t n) { return mk(x.neg, shr(x.mag, n)); }

// Exact division by 3; the caller guarantees no remainder.
SNat sDiv3(const SNat& x) {
  Nat q, r;
  divrem(x.mag, Nat(3), q, r);
  assert(r.isZero());
  return mk(x.neg, q);
}

// Exact division by 5; the caller guarantees no remainder. Toom-4 only.
SNat sDiv5(const SNat& x) {
  Nat q, r;
  divrem(x.mag, Nat(5), q, r);
  assert(r.isZero());
  return mk(x.neg, q);
}

} // namespace

Nat mulKaratsuba(const Nat& a, const Nat& b) {
  const size_t k = (max(a.size(), b.size()) + 1) / 2;
  Nat a0, a1, b0, b1;
  split(a, k, a0, a1);
  split(b, k, b0, b1);

  // (a1*B + a0)(b1*B + b0) = a1b1*B^2 + ((a0+a1)(b0+b1) - a1b1 - a0b0)*B + a0b0
  Nat z0 = mul(a0, b0);
  Nat z2 = mul(a1, b1);
  Nat z1 = sub(sub(mul(add(a0, a1), add(b0, b1)), z2), z0);

  Nat r = add(z0, shl(z1, k * 64));
  r = add(r, shl(z2, 2 * k * 64));
  return r;
}

// Evaluate A(x)=a2x^2+a1x+a0, B(x)=b2x^2+b1x+b0 at x in {0,1,-1,2,inf}, then
// interpolate the degree-4 product's coefficients. Derived from first
// principles (see the plan this was implemented from) and checked
// numerically against a worked example before implementation.
Nat mulToom3(const Nat& a, const Nat& b) {
  const size_t k = (max(a.size(), b.size()) + 2) / 3;
  Nat a0, a1, a2, b0, b1, b2, rest;
  split(a, k, a0, rest); split(rest, k, a1, a2);
  split(b, k, b0, rest); split(rest, k, b1, b2);

  const SNat pa0 = sOf(a0), pa1 = sOf(a1), pa2 = sOf(a2);
  const SNat pb0 = sOf(b0), pb1 = sOf(b1), pb2 = sOf(b2);

  const SNat p1 = sAdd(sAdd(pa0, pa1), pa2);
  const SNat q1 = sAdd(sAdd(pb0, pb1), pb2);
  const SNat pm1 = sSub(sAdd(pa0, pa2), pa1);
  const SNat qm1 = sSub(sAdd(pb0, pb2), pb1);
  const SNat p2 = sAdd(sAdd(pa0, sShl(pa1, 1)), sShl(pa2, 2));
  const SNat q2 = sAdd(sAdd(pb0, sShl(pb1, 1)), sShl(pb2, 2));

  // The five evaluation-point products are independent of each other -- run
  // them across the shared thread budget (Parallel.h) the same way Gcd.cpp's
  // matCompose runs its 8 independent sub-multiplies. Only the raw magnitude
  // mul() is parallelized; the sign combination below is negligible next to
  // it, so it stays sequential.
  Nat m0, minf, m1, mm1, m2;
  const bool big = a.size() + b.size() >= PARALLEL_MIN_LIMBS;
  runParallel({
    [&]{ m0 = mul(a0, b0); }, [&]{ minf = mul(a2, b2); },
    [&]{ m1 = mul(p1.mag, q1.mag); }, [&]{ mm1 = mul(pm1.mag, qm1.mag); },
    [&]{ m2 = mul(p2.mag, q2.mag); },
  }, big);
  const SNat v0 = sOf(m0);
  const SNat vinf = sOf(minf);
  const SNat v1 = mk(p1.neg != q1.neg, m1);
  const SNat vm1 = mk(pm1.neg != qm1.neg, mm1);
  const SNat v2 = mk(p2.neg != q2.neg, m2);

  const SNat c0 = v0;
  const SNat c4 = vinf;
  const SNat c2 = sSub(sSub(sShr(sAdd(v1, vm1), 1), c0), c4);
  const SNat numerator3 =
      sSub(sSub(sSub(sAdd(v2, sMulSmall(c0, 3)), sMulSmall(c4, 12)), sMulSmall(v1, 3)), vm1);
  const SNat c3 = sDiv3(sShr(numerator3, 1));
  const SNat c1 = sSub(sShr(sSub(v1, vm1), 1), c3);

  assert(!c1.neg && !c2.neg && !c3.neg);

  Nat r = c0.mag;
  r = add(r, shl(c1.mag, k * 64));
  r = add(r, shl(c2.mag, 2 * k * 64));
  r = add(r, shl(c3.mag, 3 * k * 64));
  r = add(r, shl(c4.mag, 4 * k * 64));
  return r;
}

// Evaluate A(x)=a3x^3+a2x^2+a1x+a0, B(x) likewise, at x in
// {0,1,-1,2,-2,"1/2",inf} -- 7 points for the degree-6 product of two
// degree-3 splits. "1/2" is evaluated fraction-free as Ph=8a0+4a1+2a2+a3,
// Qh=8b0+4b1+2b2+b3, giving rh=Ph*Qh=64*p(1/2)q(1/2); this point set keeps
// every interpolation division small (2,3,4,5,8) rather than the larger
// divisors a "3" 7th point would need. Derived from first principles and
// checked numerically against a worked example (a=[1,2,3,4], b=[5,6,7,8],
// true c0..c6 = 5,16,34,60,61,52,32) before implementation, the same way
// mulToom3's formula was.
Nat mulToom4(const Nat& a, const Nat& b) {
  const size_t k = (max(a.size(), b.size()) + 3) / 4;
  Nat a0, a1, a2, a3, b0, b1, b2, b3, rest, rest2;
  split(a, k, a0, rest); split(rest, k, a1, rest2); split(rest2, k, a2, a3);
  split(b, k, b0, rest); split(rest, k, b1, rest2); split(rest2, k, b2, b3);

  const SNat pa0 = sOf(a0), pa1 = sOf(a1), pa2 = sOf(a2), pa3 = sOf(a3);
  const SNat pb0 = sOf(b0), pb1 = sOf(b1), pb2 = sOf(b2), pb3 = sOf(b3);

  const SNat p1 = sAdd(sAdd(pa0, pa1), sAdd(pa2, pa3));
  const SNat q1 = sAdd(sAdd(pb0, pb1), sAdd(pb2, pb3));
  const SNat pm1 = sSub(sAdd(pa0, pa2), sAdd(pa1, pa3));
  const SNat qm1 = sSub(sAdd(pb0, pb2), sAdd(pb1, pb3));
  const SNat p2 = sAdd(sAdd(pa0, sShl(pa1, 1)), sAdd(sShl(pa2, 2), sShl(pa3, 3)));
  const SNat q2 = sAdd(sAdd(pb0, sShl(pb1, 1)), sAdd(sShl(pb2, 2), sShl(pb3, 3)));
  const SNat pm2 = sSub(sAdd(pa0, sShl(pa2, 2)), sAdd(sShl(pa1, 1), sShl(pa3, 3)));
  const SNat qm2 = sSub(sAdd(pb0, sShl(pb2, 2)), sAdd(sShl(pb1, 1), sShl(pb3, 3)));
  const SNat ph = sAdd(sAdd(sShl(pa0, 3), sShl(pa1, 2)), sAdd(sShl(pa2, 1), pa3));
  const SNat qh = sAdd(sAdd(sShl(pb0, 3), sShl(pb1, 2)), sAdd(sShl(pb2, 1), pb3));

  // The seven evaluation-point products are independent -- parallelize them
  // the same way mulToom3 does (see Parallel.h). (Measurement showed making
  // these sequential instead is worse, not better -- see the TOOM4_LIMBS
  // comment below for why mulToom4 is not actually reached in practice.)
  Nat m0, minf, m1, mm1, m2, mm2, mh;
  const bool big = a.size() + b.size() >= PARALLEL_MIN_LIMBS;
  runParallel({
    [&]{ m0 = mul(a0, b0); }, [&]{ minf = mul(a3, b3); },
    [&]{ m1 = mul(p1.mag, q1.mag); }, [&]{ mm1 = mul(pm1.mag, qm1.mag); },
    [&]{ m2 = mul(p2.mag, q2.mag); }, [&]{ mm2 = mul(pm2.mag, qm2.mag); },
    [&]{ mh = mul(ph.mag, qh.mag); },
  }, big);
  const SNat r0 = sOf(m0);
  const SNat rinf = sOf(minf);
  const SNat r1 = mk(p1.neg != q1.neg, m1);
  const SNat rm1 = mk(pm1.neg != qm1.neg, mm1);
  const SNat r2 = mk(p2.neg != q2.neg, m2);
  const SNat rm2 = mk(pm2.neg != qm2.neg, mm2);
  const SNat rh = mk(ph.neg != qh.neg, mh);

  const SNat c0 = r0;
  const SNat c6 = rinf;

  const SNat A1 = sSub(sSub(r1, c0), c6);
  const SNat Am1 = sSub(sSub(rm1, c0), c6);
  const SNat A2 = sSub(sSub(r2, c0), sMulSmall(c6, 64));
  const SNat Am2 = sSub(sSub(rm2, c0), sMulSmall(c6, 64));
  const SNat Ah = sSub(sSub(rh, sMulSmall(c0, 64)), c6);

  const SNat S1 = sShr(sAdd(A1, Am1), 1);
  const SNat D1 = sShr(sSub(A1, Am1), 1);
  const SNat S2 = sShr(sAdd(A2, Am2), 3);
  const SNat D2 = sShr(sSub(A2, Am2), 2);

  const SNat c4 = sDiv3(sSub(S2, S1));
  const SNat c2 = sSub(S1, c4);

  const SNat D3 = sShr(sSub(sSub(Ah, sMulSmall(c2, 16)), sMulSmall(c4, 4)), 1);

  const SNat E1 = sDiv3(sSub(D2, D1));
  const SNat E2 = sDiv3(sSub(D3, D1));
  const SNat F = sDiv5(sSub(E2, E1));
  const SNat G = sSub(D1, E1);

  const SNat c5 = sDiv3(sSub(F, G));
  const SNat c1 = sAdd(c5, F);
  const SNat c3 = sSub(sSub(D1, c1), c5);

  assert(!c1.neg && !c2.neg && !c3.neg && !c4.neg && !c5.neg);

  Nat r = c0.mag;
  r = add(r, shl(c1.mag, k * 64));
  r = add(r, shl(c2.mag, 2 * k * 64));
  r = add(r, shl(c3.mag, 3 * k * 64));
  r = add(r, shl(c4.mag, 4 * k * 64));
  r = add(r, shl(c5.mag, 5 * k * 64));
  r = add(r, shl(c6.mag, 6 * k * 64));
  return r;
}

Nat mul(const Nat& a, const Nat& b) {
  if (a.isZero() || b.isZero()) { return Nat{}; }
  if (min(a.size(), b.size()) < KARATSUBA_LIMBS) { return mulSchoolbook(a, b); }
  if (min(a.size(), b.size()) < TOOM3_LIMBS) { return mulKaratsuba(a, b); }
  if (min(a.size(), b.size()) < TOOM4_LIMBS) { return mulToom3(a, b); }
  return mulToom4(a, b);
}

// ---------------------------------------------------------------------------
// Division: Knuth TAOCP vol 2, algorithm D
// ---------------------------------------------------------------------------

namespace {

// u[0..n) -= q * v[0..n). Returns the amount still to be taken out of u[n],
// i.e. the top limb of q*v plus any borrow left over from the low limbs.
u64 mulSub(u64* u, const u64* v, size_t n, u64 q) {
  u64 carry = 0;          // top half of the running product
  unsigned char bw = 0;   // borrow chain across the subtraction
  for (size_t i = 0; i < n; ++i) {
    u64 hi, lo = mul64(q, v[i], &hi);
    // Fold the previous product's high half into this limb.
    unsigned char c = addc(0, lo, carry, &lo);
    carry = hi + c;
    bw = subb(bw, u[i], lo, &u[i]);
  }
  return carry + bw;      // carry < 2^64-1 here, so this cannot wrap
}

// u[0..n) += v[0..n); returns carry out.
u64 addBack(u64* u, const u64* v, size_t n) {
  unsigned char c = 0;
  for (size_t i = 0; i < n; ++i) { c = addc(c, u[i], v[i], &u[i]); }
  return c;
}

} // namespace

void divrem(const Nat& u, const Nat& v, Nat& q, Nat& r) {
  assert(!v.isZero());

  if (cmp(u, v) < 0) { q = Nat{}; r = u; return; }

  // Single-limb divisor: plain long division.
  if (v.size() == 1) {
    const u64 d = v.w[0];
    q.w.assign(u.size(), 0);
    u64 rem = 0;
    for (size_t i = u.size(); i-- > 0; ) {
      // rem < d is the precondition of div128by64
      q.w[i] = div128by64(rem, u.w[i], d, &rem);
    }
    q.norm();
    r = Nat(rem);
    return;
  }

  // Normalise so the divisor's top bit is set; div128by64 needs that.
  const int s = clz64(v.w.back());
  Nat vn = shl(v, s);
  Nat un = shl(u, s);
  const size_t n = vn.size();          // shl by < 64 keeps the divisor's length
  assert(n == v.size() && n >= 2);
  // Algorithm D wants exactly one spare high limb on the dividend.
  un.w.resize(u.size() + 1, 0);
  const size_t m = u.size() - n;

  q.w.assign(m + 1, 0);

  for (size_t j = m + 1; j-- > 0; ) {
    const u64 top = un.w[j + n];
    u64 qhat, rhat;
    bool correct = true;
    if (top >= vn.w[n - 1]) {
      // The true digit is at or above the base; clamp to B-1 and let the
      // add-back below fix any remaining overshoot.
      qhat = ~u64(0);
      rhat = un.w[j + n - 1] + vn.w[n - 1];
      if (rhat < vn.w[n - 1]) { correct = false; }   // rhat overflowed
    } else {
      qhat = div128by64(top, un.w[j + n - 1], vn.w[n - 1], &rhat);
    }
    // Refine qhat using the next divisor limb; at most two decrements.
    while (correct) {
      u64 hi, lo = mul64(qhat, vn.w[n - 2], &hi);
      if (hi < rhat || (hi == rhat && lo <= un.w[j + n - 2])) { break; }
      --qhat;
      u64 next = rhat + vn.w[n - 1];
      if (next < rhat) { break; }        // overflow: qhat is now small enough
      rhat = next;
    }

    u64 borrow = mulSub(&un.w[j], vn.w.data(), n, qhat);
    unsigned char bw = subb(0, un.w[j + n], borrow, &un.w[j + n]);
    if (bw) {
      // qhat was one too large: give one back.
      --qhat;
      u64 c = addBack(&un.w[j], vn.w.data(), n);
      un.w[j + n] += c;
    }
    q.w[j] = qhat;
  }

  q.norm();
  un.w.resize(n);
  un.norm();
  r = shr(un, s);
}

Nat mod(const Nat& u, const Nat& v) {
  Nat q, r;
  divrem(u, v, q, r);
  return r;
}

// ---------------------------------------------------------------------------
// Conversions
// ---------------------------------------------------------------------------

u64 modU64(const Nat& a, u64 m) {
  assert(m);
  u64 rem = 0;
  for (size_t i = a.size(); i-- > 0; ) {
    // rem is a remainder mod m, so rem < m: the precondition of div128by64.
    div128by64(rem, a.w[i], m, &rem);
  }
  return rem;
}

Nat mulMod(const Nat& a, const Nat& b, const Nat& m) { return mod(mul(a, b), m); }

Nat powMod(Nat base, const Nat& e, const Nat& m) {
  Nat result(1);
  base = mod(base, m);
  for (size_t i = e.bits(); i-- > 0; ) {
    result = mulMod(result, result, m);
    if (e.bit(i)) { result = mulMod(result, base, m); }
  }
  return result;
}

bool isProbablePrime(const Nat& n) {
  if (n.isZero() || n.isOne()) { return false; }
  static const u64 SMALL[] = {2,3,5,7,11,13,17,19,23,29,31,37};
  for (u64 s : SMALL) {
    Nat sp(s);
    if (cmp(n, sp) == 0) { return true; }
    if (mod(n, sp).isZero()) { return false; }
  }

  // n - 1 = d * 2^r with d odd
  const Nat one(1);
  Nat d = sub(n, one);
  size_t r = 0;
  while (!d.bit(0)) { d = shr(d, 1); ++r; }

  const Nat nm1 = sub(n, one);
  for (u64 s : SMALL) {
    Nat x = powMod(Nat(s), d, n);
    if (x.isOne() || cmp(x, nm1) == 0) { continue; }
    bool witness = true;
    for (size_t i = 1; i < r; ++i) {
      x = mulMod(x, x, n);
      if (cmp(x, nm1) == 0) { witness = false; break; }
    }
    if (witness) { return false; }
  }
  return true;
}

Nat mersenne(u32 p) {
  Nat r;
  r.w.assign((p + 63) / 64, ~u64(0));
  const u32 tail = p % 64;
  if (tail) { r.w.back() = (u64(1) << tail) - 1; }
  r.norm();
  return r;
}

bool fromDecimal(const string& s, Nat& out) {
  if (s.empty()) { return false; }
  Nat r;
  const Nat ten(10);
  for (char ch : s) {
    if (ch < '0' || ch > '9') { return false; }
    r = add(mul(r, ten), Nat(u64(ch - '0')));
  }
  out = std::move(r);
  return true;
}

Nat fromWords(const Words& words) {
  Nat r;
  r.w.assign((words.size() + 1) / 2, 0);
  for (size_t i = 0; i < words.size(); ++i) {
    r.w[i / 2] |= u64(words[i]) << (32 * (i % 2));
  }
  r.norm();
  return r;
}
