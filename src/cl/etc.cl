// Copyright (C) Mihai Preda

#include "base.cl"

#if READRESIDUE

// Because the data "in" is stored transposed, and we want to read
// a number of logically successive values, we have a very bad read access pattern
KERNEL(32) readResidue(P(Word2) out, CP(Word2) in) {
  u32 me = get_local_id(0);
  u32 k = (ND - 16 + me) % ND;
  u32 y = k % BIG_HEIGHT;
  u32 x = k / BIG_HEIGHT;
  out[me] = in[WIDTH * y + x];
}
#endif

#if SUM64
KERNEL(64) sum64(global ulong* out, u32 sizeBytes, global ulong* in) {
  if (get_global_id(0) == 0) { out[0] = 0; }

  ulong sum = 0;
  for (i32 p = get_global_id(0); p < sizeBytes / sizeof(u64); p += get_global_size(0)) {
    sum += in[p];
  }
  u32 prev = atomic_add((global u32*)out, (u32) sum);
  u32 high = (sum + prev) >> 32;
  atomic_add(((global u32*)out) + 1, high);
}
#endif

#if ISEQUAL
// outEqual must be "true" on entry.
KERNEL(256) isEqual(global i64 *in1, global i64 *in2, P(int) outEqual) {
  for (i32 p = get_global_id(0); p < NWORDS * sizeof(Word) / sizeof(i64); p += get_global_size(0)) {
    if (in1[p] != in2[p]) {
      *outEqual = 0;
      return;
    }
  }
}
#endif

#if SUBSMALL
// Subtract a small constant from the residue's least significant word.
//
// P+1's Lucas ladder needs V_m*V_n - P with P a small integer seed. That is the
// same shape as the Lucas-Lehmer x^2-2 step, which the carry kernels implement
// by injecting -2 as the initial carry-in. Doing it here instead leaves carry.cl
// and carryfused.cl -- five FFT-type variants apiece, all currently passing
// every gate -- completely untouched.
//
// Word 0 is nudged directly and left slightly outside the balanced range; the
// next carry pass normalises it, exactly as the injected -2 is. For a seed of a
// few units against ~39-bit words this is far below the rounding margin.
//
// Layout: data is stored transposed, logical word k at
// WIDTH*(k % BIG_HEIGHT) + k / BIG_HEIGHT (see readResidue above), so logical
// word 0 lands at index 0 and needs no address arithmetic.
KERNEL(1) subSmall(P(Word2) io, u32 val) {
  if (get_global_id(0) == 0) { io[0].x -= (Word) val; }
}
#endif

#if SUBWORDS
// out = a - b, word by word, with NO carry propagation.
//
// Stage 2's accumulator step is acc *= (A_m - T_j). Both inputs are in the
// balanced representation, so each output word lands in roughly twice the usual
// range -- one extra bit -- and is deliberately not renormalised here. The
// difference is consumed immediately by a multiply, whose fftP reads it as plain
// integer data and whose carry pass normalises the product anyway. The
// alternative, a real carry pass on the difference, would cost an extra
// transform round trip per pair and add well over half again to stage 2.
//
// The price is that the multiply sees an input one bit larger than a squaring
// does, which raises its round-off error by about one bit. That headroom is
// measured, not assumed -- see the ROE figures in the stage-2 self-test.
//
// Layout-agnostic: the transposed storage is a permutation of the logical
// words, and subtraction is elementwise, so it commutes with any permutation.
KERNEL(256) subWords(P(Word2) out, CP(Word2) a, CP(Word2) b) {
  u32 p = get_global_id(0);
  Word2 x = a[p];
  Word2 y = b[p];
  out[p] = U2(x.x - y.x, x.y - y.y);
}
#endif

#if TEST_KERNEL
// Generate a small unused kernel so developers can look at how well individual macros assemble and optimize
kernel void testKernel(global int* in, global double* out) {
  const double TAB[8] = {M_PI/13, M_PI/17, M_PI, M_SQRT2, M_SQRT1_2, M_PI/7, M_PI*7, M_PI/15};

  int me = get_local_id(0);
  int p = me * in[me] % 8; // % 15;
  out[me] = TAB[p];
}
#endif
