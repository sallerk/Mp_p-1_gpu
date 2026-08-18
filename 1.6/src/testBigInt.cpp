// Copyright (C) Mp_p-1_gpu
//
// Gate G2: BigInt and GCD correctness.
//
// Structured so every claim has an independent check:
//   * Karatsuba is checked against schoolbook, which is short enough to read.
//   * divrem is checked against its defining identity u == q*v + r, r < v.
//   * gcdLehmer is checked against gcdEuclid, which is four lines.
//   * and both are checked against an identity that needs no reference at all:
//         gcd(2^a - 1, 2^b - 1) == 2^gcd(a,b) - 1
//     which gives exact expected answers at millions of bits.

#include "BigInt.h"
#include "Gcd.h"
#include "PM1.h"
#include "timeutil.h"

#include <algorithm>
#include <cstdio>
#include <string>

using namespace std;

namespace {

int failures = 0;
int checks = 0;

void check(bool ok, const string& what) {
  ++checks;
  if (!ok) {
    ++failures;
    printf("   FAIL  %s\n", what.c_str());
  }
}

// xorshift64*, so runs are reproducible.
struct Rng {
  u64 s;
  explicit Rng(u64 seed) : s{seed ? seed : 1} {}
  u64 next() {
    s ^= s >> 12; s ^= s << 25; s ^= s >> 27;
    return s * 0x2545F4914F6CDD1Dull;
  }
};

Nat randomNat(Rng& rng, size_t limbs) {
  Nat r;
  r.w.resize(limbs);
  for (size_t i = 0; i < limbs; ++i) { r.w[i] = rng.next(); }
  r.norm();
  return r;
}

u32 gcdU32(u32 a, u32 b) { while (b) { u32 t = a % b; a = b; b = t; } return a; }

} // namespace

int runBigIntTests() {
  printf("G2: BigInt and GCD\n\n");
  Rng rng{0x9E3779B97F4A7C15ull};

  // --- shifts round-trip ---------------------------------------------------
  printf("  shifts\n");
  for (int t = 0; t < 200; ++t) {
    Nat a = randomNat(rng, 1 + rng.next() % 12);
    size_t k = rng.next() % 200;
    check(shr(shl(a, k), k) == a, "shr(shl(a,k),k) == a");
  }

  // --- add / sub round-trip ------------------------------------------------
  printf("  add/sub\n");
  for (int t = 0; t < 200; ++t) {
    Nat a = randomNat(rng, 1 + rng.next() % 12);
    Nat b = randomNat(rng, 1 + rng.next() % 12);
    check(sub(add(a, b), b) == a, "sub(add(a,b),b) == a");
  }

  // --- Karatsuba vs schoolbook --------------------------------------------
  printf("  multiply (Karatsuba vs schoolbook)\n");
  for (int t = 0; t < 60; ++t) {
    Nat a = randomNat(rng, 1 + rng.next() % 160);
    Nat b = randomNat(rng, 1 + rng.next() % 160);
    check(mul(a, b) == mulSchoolbook(a, b), "mul == mulSchoolbook");
  }

  // --- Toom-3 vs Karatsuba --------------------------------------------------
  // Called directly (not through mul()'s size dispatch) so the edge cases
  // below -- especially the very unequal sizes -- actually exercise
  // mulToom3's own splitting/interpolation, regardless of where TOOM3_LIMBS
  // happens to be set.
  printf("  multiply (Toom-3 vs Karatsuba)\n");
  for (int t = 0; t < 60; ++t) {
    Nat a = randomNat(rng, 1 + rng.next() % 700);
    Nat b = randomNat(rng, 1 + rng.next() % 700);
    check(mulToom3(a, b) == mulKaratsuba(a, b), "mulToom3 == mulKaratsuba");
  }
  {
    // Unequal sizes, sizes not divisible by 3, and an operand smaller than a
    // third of the other's length (so a2 or b2 truncates to zero after the
    // 3-way split) -- exactly the cases a hand-derived interpolation is most
    // likely to get wrong.
    static const size_t pairs[][2] = {
      {301, 301}, {301, 302}, {301, 303}, {900, 1}, {1, 900}, {900, 100},
      {100, 900}, {450, 450}, {451, 449}, {600, 600}, {3, 900}, {900, 3},
    };
    for (auto& pr : pairs) {
      Nat a = randomNat(rng, pr[0]), b = randomNat(rng, pr[1]);
      check(mulToom3(a, b) == mulKaratsuba(a, b),
            "mulToom3 == mulKaratsuba (" + to_string(pr[0]) + "x" + to_string(pr[1]) + ")");
    }
  }

  // --- Toom-4 vs Toom-3 ------------------------------------------------------
  // Called directly (not through mul()'s size dispatch), same reasoning as
  // the Toom-3 block above.
  printf("  multiply (Toom-4 vs Toom-3)\n");
  for (int t = 0; t < 60; ++t) {
    Nat a = randomNat(rng, 1 + rng.next() % 900);
    Nat b = randomNat(rng, 1 + rng.next() % 900);
    check(mulToom4(a, b) == mulToom3(a, b), "mulToom4 == mulToom3");
  }
  {
    // Unequal sizes, sizes not divisible by 4, and an operand smaller than a
    // quarter of the other's length (so a2/a3 or b2/b3 truncate to zero after
    // the 4-way split) -- exactly the cases a hand-derived interpolation is
    // most likely to get wrong.
    static const size_t pairs[][2] = {
      {401, 401}, {401, 402}, {401, 403}, {401, 404}, {1200, 1}, {1, 1200},
      {1200, 100}, {100, 1200}, {600, 600}, {601, 599}, {4, 1200}, {1200, 4},
    };
    for (auto& pr : pairs) {
      Nat a = randomNat(rng, pr[0]), b = randomNat(rng, pr[1]);
      check(mulToom4(a, b) == mulToom3(a, b),
            "mulToom4 == mulToom3 (" + to_string(pr[0]) + "x" + to_string(pr[1]) + ")");
    }
  }

  // --- divrem defining identity -------------------------------------------
  printf("  divrem (u == q*v + r, r < v)\n");
  for (int t = 0; t < 300; ++t) {
    Nat u = randomNat(rng, 1 + rng.next() % 30);
    Nat v = randomNat(rng, 1 + rng.next() % 12);
    if (v.isZero()) { continue; }
    Nat q, r;
    divrem(u, v, q, r);
    check(add(mul(q, v), r) == u, "u == q*v + r");
    check(cmp(r, v) < 0, "r < v");
  }

  // --- gcdLehmer vs gcdEuclid ---------------------------------------------
  printf("  gcdLehmer vs gcdEuclid\n");
  for (int t = 0; t < 200; ++t) {
    Nat a = randomNat(rng, 1 + rng.next() % 20);
    Nat b = randomNat(rng, 1 + rng.next() % 20);
    if (a.isZero() || b.isZero()) { continue; }
    Nat g1 = gcdEuclid(a, b);
    Nat g2 = gcdLehmer(a, b);
    check(g1 == g2, "gcdLehmer == gcdEuclid (a=" + a.hex() + " b=" + b.hex() + ")");
  }

  // --- gcd with a shared planted factor -----------------------------------
  printf("  gcd with a planted common factor\n");
  for (int t = 0; t < 60; ++t) {
    Nat g = randomNat(rng, 1 + rng.next() % 6);
    Nat x = randomNat(rng, 1 + rng.next() % 10);
    Nat y = randomNat(rng, 1 + rng.next() % 10);
    if (g.isZero() || x.isZero() || y.isZero()) { continue; }
    Nat a = mul(g, x), b = mul(g, y);
    Nat got = gcdLehmer(a, b);
    // g divides the gcd, and the gcd divides both.
    check(mod(got, g).isZero() || mod(got, g).isZero(), "planted g divides gcd");
    check(mod(a, got).isZero() && mod(b, got).isZero(), "gcd divides both");
  }

  // --- the Mersenne identity, needing no reference -------------------------
  printf("  gcd(2^a-1, 2^b-1) == 2^gcd(a,b)-1\n");
  static const u32 pairs[][2] = {
    {12, 18}, {100, 60}, {521, 607}, {1279, 2203}, {9941, 11213},
    {86243, 110503}, {216091, 132049}, {859433, 1257787},
  };
  for (auto& pr : pairs) {
    Nat a = mersenne(pr[0]), b = mersenne(pr[1]);
    Nat want = mersenne(gcdU32(pr[0], pr[1]));
    Nat got = gcdLehmer(a, b);
    check(got == want, "mersenne gcd " + to_string(pr[0]) + "," + to_string(pr[1]));
  }

  // --- is Lehmer's fast path actually being taken? -------------------------
  // Without this, a gcdLehmer that always fell back to plain division would
  // pass every test above while delivering none of the speedup.
  printf("\n  Lehmer fast path vs Euclid fallback\n");
  {
    extern u64 gLehmerFastSteps, gLehmerFallbackDivisions;

    // NOTE on the choice of inputs. gcd(2^a-1, 2^b-1) is a poor performance
    // test even though it is an excellent correctness test: when a and b differ
    // in size the quotient is around 2^(a-b), i.e. astronomically large, and
    // Lehmer's leading-word loop only accelerates *small* quotients. Those
    // cases are division-bound by nature.
    //
    // The real workload is gcd(x-1, 2^p-1) with x-1 a residue of about p bits,
    // so both operands are the same size and the quotient sequence is the usual
    // random one (q == 1 roughly 41% of the time). That is what is measured
    // here.
    for (size_t limbs : {size_t(400), size_t(2000), size_t(6000), size_t(12000)}) {
      Nat a = randomNat(rng, limbs);
      Nat b = randomNat(rng, limbs);
      if (cmp(a, b) < 0) { std::swap(a, b); }

      gLehmerFastSteps = gLehmerFallbackDivisions = 0;
      Timer t1;
      Nat gl = gcdLehmer(a, b);
      const double lehmerSecs = t1.at();
      const u64 fast = gLehmerFastSteps, slow = gLehmerFallbackDivisions;

      Timer t2;
      Nat ge = gcdEuclid(a, b);
      const double euclidSecs = t2.at();

      check(gl == ge, "lehmer == euclid at " + to_string(limbs) + " limbs");
      printf("     %5zu limbs (%6zu bits): fast %6llu, fallback %5llu | "
             "lehmer %6.3fs euclid %6.3fs (%.1fx)\n",
             limbs, limbs * 64,
             (unsigned long long) fast, (unsigned long long) slow,
             lehmerSecs, euclidSecs,
             lehmerSecs > 0 ? euclidSecs / lehmerSecs : 0.0);
      check(fast > 0, "Lehmer fast path is exercised at " + to_string(limbs) + " limbs");
      check(fast > slow, "Lehmer fast path dominates at " + to_string(limbs) + " limbs");
    }
  }

  // --- half-GCD: correctness against the trusted Lehmer ---------------------
  printf("\n  gcdHalf vs gcdLehmer (differential)\n");
  for (int t = 0; t < 300; ++t) {
    size_t la = 1 + rng.next() % 400, lb = 1 + rng.next() % 400;
    Nat a = randomNat(rng, la), b = randomNat(rng, lb);
    if (a.isZero() || b.isZero()) { continue; }
    check(gcdHalf(a, b) == gcdLehmer(a, b), "gcdHalf == gcdLehmer");
  }
  for (auto& pr : pairs) {
    Nat a = mersenne(pr[0]), b = mersenne(pr[1]);
    Nat want = mersenne(gcdU32(pr[0], pr[1]));
    check(gcdHalf(a, b) == want, "gcdHalf mersenne " + to_string(pr[0]) + "," + to_string(pr[1]));
  }
  // Planted factor, large enough to exercise the recursion.
  for (int t = 0; t < 10; ++t) {
    Nat g = randomNat(rng, 200 + rng.next() % 100);
    Nat x = randomNat(rng, 500 + rng.next() % 300);
    Nat y = randomNat(rng, 500 + rng.next() % 300);
    if (g.isZero() || x.isZero() || y.isZero()) { continue; }
    Nat a = mul(g, x), b = mul(g, y);
    check(gcdHalf(a, b) == gcdLehmer(a, b), "gcdHalf == gcdLehmer (planted, large)");
  }

  // --- half-GCD: is it actually subquadratic? ------------------------------
  // Also reports the real mul-vs-div time split (gMulNanos/gDivNanos, see
  // Gcd.h) at each scale, up to the exact bit length of M82589933 -- the
  // production-scale question this exists to answer, not extrapolated from
  // the smaller tiers below it.
  printf("\n  gcdHalf scaling vs gcdLehmer\n");
  {
    double prev = 0;
    size_t prevLimbs = 0;
    for (size_t limbs : {size_t(2000), size_t(8000), size_t(32000), size_t(128000),
                          size_t(1290468)}) {
      Nat a = randomNat(rng, limbs), b = randomNat(rng, limbs);
      if (cmp(a, b) < 0) { std::swap(a, b); }

      gMulNanos = 0; gDivNanos = 0;
      Timer t1;
      Nat gh = gcdHalf(a, b);
      const double halfSecs = t1.at();
      const double mulSecs = double(gMulNanos.load()) / 1e9;
      const double divSecs = double(gDivNanos.load()) / 1e9;
      const double otherSecs = std::max(0.0, halfSecs - mulSecs - divSecs);

      // Lehmer is quadratic; only run it where that is still affordable.
      double lehSecs = -1;
      if (limbs <= 8000) {
        Timer t2;
        Nat gl = gcdLehmer(a, b);
        lehSecs = t2.at();
        check(gh == gl, "gcdHalf == gcdLehmer at " + to_string(limbs) + " limbs");
      }

      const double growth = (prev > 0) ? halfSecs / prev : 0.0;
      const double sizeRatio = (prevLimbs > 0) ? double(limbs) / prevLimbs : 0.0;
      printf("     %7zu limbs: half %8.3fs  [mul %7.3fs (%4.1f%%) div %7.3fs (%4.1f%%) other %7.3fs (%4.1f%%)]",
             limbs, halfSecs,
             mulSecs, halfSecs > 0 ? 100.0 * mulSecs / halfSecs : 0.0,
             divSecs, halfSecs > 0 ? 100.0 * divSecs / halfSecs : 0.0,
             otherSecs, halfSecs > 0 ? 100.0 * otherSecs / halfSecs : 0.0);
      if (lehSecs >= 0) { printf("  lehmer %7.3fs (%.1fx)", lehSecs, lehSecs / halfSecs); }
      if (growth > 0)   { printf("   [%.0fx size -> %.1fx time]", sizeRatio, growth); }
      printf("\n");
      prev = halfSecs;
      prevLimbs = limbs;
    }
  }

  // --- decimal round-trip --------------------------------------------------
  printf("\n  decimal conversion\n");
  for (int t = 0; t < 100; ++t) {
    Nat a = randomNat(rng, 1 + rng.next() % 8);
    Nat back;
    check(fromDecimal(a.dec(), back) && back == a, "dec/fromDecimal round-trip");
  }

  // --- primality -----------------------------------------------------------
  printf("  Miller-Rabin\n");
  {
    const char* primes[] = {"2", "3", "97", "65537", "2147483647",
                            "67280421310721",
                            "170141183460469231731687303715884105727"};
    for (const char* s : primes) {
      Nat n; fromDecimal(s, n);
      check(isProbablePrime(n), string("prime: ") + s);
    }
    const char* comps[] = {"1", "4", "91", "65536", "4294967297",
                           "78641950567367153678977"};
    for (const char* s : comps) {
      Nat n; fromDecimal(s, n);
      check(!isProbablePrime(n), string("composite: ") + s);
    }
  }

  // --- splitting a real composite gcd --------------------------------------
  // Exactly what a live run produced: M82589959 has two factors with small k,
  // so a single gcd returned their product and was reported as "a factor".
  printf("  splitting a composite gcd (from a real M82589959 run)\n");
  {
    Nat g;
    fromDecimal("78641950567367153678977", g);
    auto fs = splitFactors(g, 82589959);
    check(fs.size() == 2, "gcd splits into 2 factors");
    const char* want[] = {"23785908193", "3306241238689"};
    const u64 wantK[] = {144, 20016};
    for (size_t i = 0; i < fs.size() && i < 2; ++i) {
      check(fs[i].value.dec() == want[i], string("factor == ") + want[i]);
      check(fs[i].prime, "reported factor is prime");
      check(fs[i].dividesMp, "reported factor divides M_p");
      check(fs[i].k == wantK[i], "k matches");
      printf("     %-16s k=%-8llu %s\n", fs[i].value.dec().c_str(),
             (unsigned long long) fs[i].k,
             (fs[i].prime && fs[i].dividesMp) ? "prime, divides M_p" : "PROBLEM");
    }
  }

  printf("\nG2: %d/%d checks passed.\n", checks - failures, checks);
  return failures == 0 ? 0 : 1;
}
