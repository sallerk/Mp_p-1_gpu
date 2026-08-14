// Copyright (C) Mp_p-1_gpu
//
// See BigInt.h. Nothing here is clever; it is meant to be obviously correct,
// with Karatsuba as the only asymptotic concession. The GCD on top is where the
// real work happens.

#include "BigInt.h"

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

} // namespace

Nat mul(const Nat& a, const Nat& b) {
  if (a.isZero() || b.isZero()) { return Nat{}; }
  if (min(a.size(), b.size()) < KARATSUBA_LIMBS) { return mulSchoolbook(a, b); }

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
