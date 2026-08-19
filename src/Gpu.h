// Copyright (C) Mihai Preda

#pragma once

#include "Background.h"
#include "Buffer.h"
#include "Context.h"
#include "Queue.h"
#include "KernelCompiler.h"

#include "Saver.h"
#include "common.h"
#include "Kernel.h"
#include "Profile.h"
#include "GpuCommon.h"
#include "FFTConfig.h"

#include <vector>
#include <memory>
#include <filesystem>
#include <cmath>
#include <functional>

struct PRPResult;
struct Task;
struct Stage2Plan;

class Signal;
class ProofSet;

using TrigBuf = Buffer<double2>;
using TrigPtr = shared_ptr<TrigBuf>;

inline u64 residue(const Words& words) { return (u64(words[1]) << 32) | words[0]; }

struct PRPResult {
  bool isPrime{};
  u64 res64 = 0;
  u32 nErrors = 0;
  fs::path proofPath{};
  std::string res2048;
};

struct LLResult {
  bool isPrime;
  u64 res64;
};

struct ZAvg {
  double sum{};
  double n{};

  void update(double z, u32 inc) {
    sum += z * inc;
    n += inc;
  }

  double avg() { return sum / n; }
};

class RoeInfo {
public:
  RoeInfo() = default;
  RoeInfo(u32 n, double max, double mean, double sd) : N{n}, max{max}, mean{mean}, sd{sd} {
    // https://en.wikipedia.org/wiki/Gumbel_distribution
    gumbelBeta = sd * 0.779696801233676; // sqrt(6)/pi
    gumbelMiu = mean - gumbelBeta * 0.577215664901533; // Euler-Mascheroni
  }

  double z(double x = 0.5) const { return N ? (x - gumbelMiu) / gumbelBeta : 0.0; }

  double gumbelCDF(double x) const { return exp(-exp(-z(x))); }
  double gumbelRightCDF(double x) const { return -expm1(-exp(-z(x))); }

  std::string toString() const;

  u32 N{};
  double max{}, mean{}, sd{};
  double gumbelMiu{}, gumbelBeta{};
};

struct Weights {
  vector<double> weightsConstIF;
  vector<double> weightsIF;
  vector<u32> bitsCF;
};

class Gpu {
  Queue* queue;
  Background* background;

public:
  const Args& args;

private:
  std::unique_ptr<Saver<PRPState>> saver;

  u32 E;
  u32 N;

  FFTConfig fft;
  u32 WIDTH;
  u32 SMALL_H;
  u32 BIG_H;

  u32 hN, nW, nH;
  bool useLongCarry;
  u32 wantROE{};

  Profile profile{};

  KernelCompiler compiler;

  /* Kernels for FFT_FP64 or FFT_FP32 */
  Kernel kfftMidIn;
  Kernel kfftHin;
  Kernel ktailSquareZero;
  Kernel ktailSquare;
  Kernel ktailMul;
  Kernel ktailMulLow;
  Kernel kfftMidOut;
  Kernel kfftW;

  /* Kernels for NTT_GF31 */
  Kernel kfftMidInGF31;
  Kernel kfftHinGF31;
  Kernel ktailSquareZeroGF31;
  Kernel ktailSquareGF31;
  Kernel ktailMulGF31;
  Kernel ktailMulLowGF31;
  Kernel kfftMidOutGF31;
  Kernel kfftWGF31;

  /* Kernels for NTT_GF61 */
  Kernel kfftMidInGF61;
  Kernel kfftHinGF61;
  Kernel ktailSquareZeroGF61;
  Kernel ktailSquareGF61;
  Kernel ktailMulGF61;
  Kernel ktailMulLowGF61;
  Kernel kfftMidOutGF61;
  Kernel kfftWGF61;

  /* Kernels dealing with the FP data and product of NTT primes */
  Kernel kfftP;
  Kernel kCarryA;
  Kernel kCarryAROE;
  Kernel kCarryM;
  Kernel kCarryMROE;
  Kernel kCarryLL;
  Kernel kCarryFused;
  Kernel kCarryFusedROE;
  Kernel kCarryFusedMul;
  Kernel kCarryFusedMulROE;
  Kernel kCarryFusedLL;

  Kernel carryB;
  Kernel transpIn, transpOut;
  Kernel readResidue;
  Kernel kSubSmall;
  Kernel kSubWords;
  Kernel kernIsEqual;
  Kernel sum64;

  /* Weird test kernels */
  Kernel testTrig;
  Kernel testFFT4;
  Kernel testFFT14;
  Kernel testFFT15;
  Kernel testFFT;
  Kernel testTime;

  // Kernel testKernel;

  // Copy of some -use options needed for Kernel, Trig, and Weights initialization
  bool tail_single_wide;                // TailSquare processes one line at a time
  bool tail_single_kernel;              // TailSquare does not use a separate kernel for line zero
  u32 in_place;                         // Should GPU perform transform in-place. 1 = nVidia friendly memory layout, 2 = AMD friendly.
  u32 pad_size;                         // Pad size in bytes as specified on the command line or config.txt.  Maximum value is 512.

  // Twiddles: trigonometry constant buffers, used in FFTs.
  // The twiddles depend only on FFT config and do not depend on the exponent.
  // It is important to generate the height trigs before the width trigs because width trigs can be a subset of the height trigs
  TrigPtr bufTrigH;
  TrigPtr bufTrigM;
  TrigPtr bufTrigW;

  // Weights and the "bigWord bits" are only needed for FP64 and FP32 FFTs
  Weights weights;
  Buffer<double> bufConstWeights;
  Buffer<double> bufWeights;
  Buffer<u32> bufBits;  // bigWord bits aligned for CarryFused/fftP

  // "integer word" buffers. These are "small buffers": N x int.
  Buffer<Word> bufData;   // Main int buffer with the words.
  Buffer<Word> bufAux;    // Auxiliary int buffer, used in transposing data in/out and in check.
  Buffer<Word> bufCheck;  // Buffers used with the error check.

  // Carry buffers, used in carry and fusedCarry.
  Buffer<i64> bufCarry;  // Carry shuttle.
  Buffer<int> bufReady;  // Per-group ready flag for stairway carry propagation.

  // Small aux buffers.
  Buffer<Word> bufSmallOut;
  Buffer<u64> bufSumOut;
  Buffer<int> bufTrue;
  Buffer<float> bufROE; // The round-off error ("ROE"), one float element per iteration.
  Buffer<float> bufStatsCarry;

  u32 roePos{};   // The next position to write in the ROE stats buffer.
  u32 carryPos{}; // The next position to write in the Carry stats buffer.

  // The ROE positions originating from multiplications (as opposed to squarings).
  vector<u32> mulRoePos;

  // Auxilliary big buffers
  Buffer<double> buf1;
  Buffer<double> buf2;
  Buffer<double> buf3;

  unsigned statsBits;
  TimeInfo* timeBufVect;
  ZAvg zAvg;

  int NUM_CACHE_GROUPS = 3;

  void fftP(Buffer<double>& out, Buffer<double>& in) { fftP(out, reinterpret_cast<Buffer<Word>&>(in)); }
  void fftP(Buffer<double>& out, Buffer<Word>& in);
  void fftMidIn(Buffer<double>& out, Buffer<double>& in, int cache_group = 0);
  void fftMidOut(Buffer<double>& out, Buffer<double>& in, int cache_group = 0);
  void fftHin(Buffer<double>& out, Buffer<double>& in);
  void tailSquare(Buffer<double>& out, Buffer<double>& in, int cache_group = 0);
  void tailMul(Buffer<double>& out, Buffer<double>& in1, Buffer<double>& in2, int cache_group = 0);
  void tailMulLow(Buffer<double>& out, Buffer<double>& in1, Buffer<double>& in2, int cache_group = 0);
  void fftW(Buffer<double>& out, Buffer<double>& in, int cache_group = 0);
  void carryA(Buffer<double>& out, Buffer<double>& in) { carryA(reinterpret_cast<Buffer<Word>&>(out), in); }
  void carryA(Buffer<Word>& out, Buffer<double>& in);
  void carryM(Buffer<Word>& out, Buffer<double>& in);
  void carryLL(Buffer<Word>& out, Buffer<double>& in);
  void carryFused(Buffer<double>& out, Buffer<double>& in);
  void carryFusedMul(Buffer<double>& out, Buffer<double>& in);
  void carryFusedLL(Buffer<double>& out, Buffer<double>& in);

  vector<Word> readWords(Buffer<Word> &buf);
  void writeWords(Buffer<Word>& buf, vector<Word> &words);

  vector<Word> readOut(Buffer<Word> &buf);
  void writeIn(Buffer<Word>& buf, vector<Word>&& words);

  enum LEAD_TYPE {LEAD_NONE = 0, LEAD_WIDTH = 1, LEAD_MIDDLE = 2};

  void square(Buffer<Word>& out, Buffer<Word>& in, enum LEAD_TYPE leadIn, enum LEAD_TYPE leadOut, bool doMul3 = false, bool doLL = false);
  void square(Buffer<Word>& io) { square(io, io, LEAD_NONE, LEAD_NONE, false, false); }
  void squareCERT(Buffer<Word>& io, enum LEAD_TYPE leadIn, enum LEAD_TYPE leadOut) { square(io, io, leadIn, leadOut, false, false); }
  void squareLL(Buffer<Word>& io, enum LEAD_TYPE leadIn, enum LEAD_TYPE leadOut) { square(io, io, leadIn, leadOut, false, true); }

  u32 squareLoop(Buffer<Word>& out, Buffer<Word>& in, u32 from, u32 to, bool doTailMul3);
  u32 squareLoop(Buffer<Word>& io, u32 from, u32 to) { return squareLoop(io, io, from, to, false); }

  bool isEqual(Buffer<Word>& bufCheck, Buffer<Word>& bufAux);
  u64 bufResidue(Buffer<Word>& buf);
  
  vector<u32> writeBase(const vector<u32> &v);
  
  void exponentiate(Buffer<Word>& bufInOut, u64 exp, Buffer<double>& buf1, Buffer<double>& buf2, Buffer<double>& buf3);

  // out = base^exp mod (2^E-1) for a machine-word exponent. Square-and-multiply
  // straight off the bits; stage 2 only uses it for setup, where the total is a
  // few thousand multiplies against millions in the main loop, so there is no
  // reason to be cleverer.
  void powSmall(Buffer<Word>& out, Buffer<Word>& base, u64 exp);

  // io := io * one mod (2^E-1), where `one` holds the constant 1 -- forces a
  // full carry pass through the ordinary modMul/carry path, using no kernel
  // beyond what already exists (the caller owns `one` as a small persistent
  // scratch buffer, the same way stage2()'s own T/A/S/C buffers are owned by
  // the call that needs them, rather than adding a permanent class member for
  // a value only P+1 stage 2 needs).
  //
  // Used on every kSubWords result in P+1 stage 2's block-step and T-table
  // build that becomes PERSISTED, reused state (unlike P-1 stage 2's diff,
  // consumed by a single immediate multiply and never touched again).
  //
  // The precaution turned out not to be load-bearing, and it is worth writing
  // down why rather than deleting the reasoning: the naive worry was that an
  // uncarried "one extra bit" width would compound across steps, since each
  // step's kSubWords result becomes the NEXT step's input. It does not,
  // because every such result is produced by kSubWords(modMul(...), prior)
  // -- one operand is always fresh out of a modMul or squareLL, whose own
  // carry pass narrows its output regardless of the OTHER operand's width. So
  // the wide operand entering any given kSubWords is always exactly ONE
  // application wide, never the sum of many, which is the same "one extra
  // bit" margin P-1's own kSubWords already relies on (measured safe there).
  // Verified directly, not just argued: --selftest=pp1stage2 runs the same
  // plan with and without these calls at a realistic exponent and up to
  // ~100000 accumulator slots, and the two accumulators come out identical.
  // Left in anyway, since it costs a small fraction of stage 2's dominant
  // per-slot cost and the alternative failure mode -- a silently wrong
  // accumulator -- is worse than the cost of a margin that turned out not to
  // be needed.
  void renormalize(Buffer<Word>& io, Buffer<Word>& one);

  void writeState(u32 k, const vector<u32>& check, u32 blockSize);

  // does either carrryFused() or the expanded version depending on useLongCarry
  void doCarry(Buffer<double>& out, Buffer<double>& in, Buffer<Word>& tmp);

  void mul(Buffer<Word>& ioA, Buffer<double>& inB, Buffer<double>& tmp1, Buffer<double>& tmp2, bool mul3 = false);
  void mul(Buffer<Word>& io, Buffer<double>& inB);

  void modMul(Buffer<Word>& ioA, Buffer<Word>& inB, bool mul3 = false);
  void modMul(Buffer<Word>& ioA, Buffer<Word>& inB, enum LEAD_TYPE leadInB, bool mul3 = false);

  fs::path saveProof(const Args& args, const ProofSet& proofSet);
  std::pair<RoeInfo, RoeInfo> readROE();
  RoeInfo readCarryStats();

  u32 updateCarryPos(u32 bit);

  PRPState loadPRP(Saver<PRPState>& saver);

  vector<Word> readChecked(Buffer<Word>& buf);

  // void measureTransferSpeed();

  static void doDiv9(u32 E, Words& words);
  static bool equals9(const Words& words);
  void selftestTrig();

public:
  Gpu(Queue* q, GpuCommon shared, FFTConfig fft, u32 E, const vector<KeyVal>& extraConf, bool logFftSize);
  static unique_ptr<Gpu> make(Queue* q, u32 E, GpuCommon shared, FFTConfig fft,
                              const vector<KeyVal>& extraConf = {}, bool logFftSize = true);

  ~Gpu();

  PRPResult isPrimePRP(const Task& task);
  LLResult isPrimeLL(const Task& task);
  array<u64, 4> isCERT(const Task& task);

  double timePRP(int quick = 7);

  tuple<bool, u64, RoeInfo, RoeInfo> measureROE(bool quick);
  tuple<bool, RoeInfo> measureCarry();

  Saver<PRPState> *getSaver();

  void writeIn(Buffer<Word>& buf, const vector<u32> &words);
  
  u64 dataResidue()  { return bufResidue(bufData); }
  u64 checkResidue() { return bufResidue(bufCheck); }

  bool doCheck(u32 blockSize);

  void logTimeKernels();

  Words readAndCompress(Buffer<Word>& buf);
  vector<u32> readCheck();
  vector<u32> readData();


  u32 getFFTSize() { return N; }

  // Bytes for ONE stage-2 T_j buffer, so the memory budget can size the pairing
  // table exactly instead of estimating it. An earlier version guessed
  // exponent/18 words of 4 bytes and was wrong twice over -- Word is 64-bit, and
  // these NTTs carry ~39 bits/word, not 18.
  u64 stage2BufferBytes() const { return u64(N) * sizeof(Word); }

  // return A^h * B
  Words expMul(const Words& A, u64 h, const Words& B, bool doSquareB);

  // return A^h * B^2
  Words expMul2(const Words& A, u64 h, const Words& B);

  // A:= A^h * B
  void expMul(Buffer<Word>& A, u64 h, Buffer<Word>& B);

  // return A^(2^n)
  Words expExp2(const Words& A, u32 n);

  // P-1 stage 1: 3^E mod (2^p - 1), E as little-endian 64-bit limbs.
  // progress(done, total) is called every reportEvery squarings; returning
  // false stops the ladder early (used for Ctrl-C / checkpointing).
  // resumeFrom/resumeBit continue a checkpointed ladder: the residue already
  // accounts for every bit above resumeBit. saveEvery/save checkpoint as it
  // goes, handing back (residue, next bit still to process).
  Words powBase3(const vector<u64>& expLimbs, u32 reportEvery,
                 const std::function<bool(u64, u64)>& progress,
                 const Words* resumeFrom = nullptr, u64 resumeBit = 0,
                 u32 saveEvery = 0,
                 const std::function<void(const Words&, u64)>& save = {});

  // P+1 stage 1: V_exp(seed, 1) mod (2^E - 1), the Lucas sequence with Q = 1.
  //
  // Montgomery ladder on the pair (V_k, V_{k+1}), whose difference is always
  // V_1 = seed. Every step is a squaring or a multiply followed by subtracting
  // a SMALL CONSTANT (2 or seed) -- no full-residue subtraction is ever needed,
  // which is what makes this implementable without a general subtract.
  //   V_2k   = V_k^2 - 2
  //   V_2k+1 = V_k * V_k+1 - seed
  // Cost is one squaring plus one multiply per bit, so about 2x P-1 for the
  // same B1.
  Words lucasV(u32 seed, const vector<u64>& expLimbs, u32 reportEvery,
               const std::function<bool(u64, u64)>& progress,
               const Words* resumeA = nullptr, const Words* resumeB = nullptr,
               u64 resumeBit = 0, u32 saveEvery = 0,
               const std::function<void(const Words&, const Words&, u64)>& save = {});

  // base^exp mod (2^E - 1) for an arbitrary base -- used to extend a completed
  // stage 1 to a larger B1 without redoing it. Same resume/checkpoint contract
  // as powBase3.
  Words powResidue(const Words& base, const vector<u64>& expLimbs,
                   u32 reportEvery,
                   const std::function<bool(u64, u64)>& progress,
                   const Words* resumeFrom = nullptr, u64 resumeBit = 0,
                   u32 saveEvery = 0,
                   const std::function<void(const Words&, u64)>& save = {});
  // P-1 stage 2. `x` is the stage-1 residue; the result is the accumulated
  // product, whose gcd with M_p carries any factor q = 2kp+1 whose k is
  // B1-smooth apart from a single prime in (B1, B2].
  //
  // Walks the plan's open (m, j) slots doing acc *= (A_m - T_j), with
  // A_m = x^((m*D)^2) carried forward by A *= S, S *= C rather than recomputed.
  // progress(done, total) is called every reportEvery multiplies and counts the
  // recurrence as well as the accumulator, so the rate and ETA are honest;
  // returning false stops early and returns the partial accumulator.
  //
  // If mulRoeOut is given, round-off error is collected over the run. The
  // difference A_m - T_j is not carry-normalised, so these multiplies see inputs
  // about one bit wider than a squaring does -- this is how that cost is
  // measured rather than assumed.
  // normalizeDiff is for verification only. It sends each difference through a
  // GPU->CPU->GPU round trip, which carries it back into the balanced range
  // without changing its value. Running a plan both ways and comparing the
  // accumulators bit for bit is what proves the un-normalised operand is
  // harmless at a given exponent; it is far too slow for production.
  // Where a stage-2 walk is, so it can be checkpointed and resumed. The plan is
  // not part of it: the plan is a pure function of (b1, b2, d, w), so only the
  // position and the three live residues need saving.
  struct Stage2Pos {
    u64 m = 0;
    u64 jIdx = 0;
    u64 done = 0;
    Words acc, a, s;
  };

  // Fills T[i] with x^(jset[i]^2) for the whole plan, by a second-difference
  // chain stepping j by 2 rather than one exponentiation per entry:
  //
  //   x^((j+2)^2) = x^(j^2) * x^(4j+4)      x^(4(j+2)+4) = x^(4j+4) * x^8
  //
  // Two multiplies per step of 2 in j -- ONE per unit of j, and independent of
  // how many j are actually stored. d is even, so every j coprime to d is odd
  // and the chain lands on all of them; jset[0] is always 1, where it starts.
  //
  // This is where P-1 borrows from the Lucas-trace stage 2 (P+1's here, and
  // PrMers' "V-trace" for P-1): the win in those is not the pairing, which is
  // no better than this one's, but that a LINEAR index lets the baby table be
  // built incrementally. x^(j^2) is quadratic in j, so it needs a second
  // difference rather than a first, and that is the only difference -- no
  // modular inverse required, which is what a true V-trace port would cost
  // here (x + x^-1 needs an extended gcd on a p-bit number, dearer than the
  // whole of stage 2 it would be speeding up).
  //
  // T is sized by the caller and owns the buffers; xBuf must already hold x.
  void buildStage2Table(vector<Buffer<Word>>& T, Buffer<Word>& xBuf,
                        const Stage2Plan& plan);

  // Builds the table with buildStage2Table and checks every entry against a
  // direct powSmall(x, j*j), returning the number that differ -- 0 means the
  // chain reproduces independent exponentiation exactly. Direct exponentiation
  // is what shipped through 1.6, so this tests the new code against a version
  // that has real production mileage, not against a fresh reference. Selftest
  // only: it holds the whole table plus a scratch pair, so keep the plan small
  // enough that J residues fit.
  u32 stage2TableCheck(const Words& x, const Stage2Plan& plan);

  // accSeed starts the accumulator at a value other than 1. That is all a B2
  // extension needs: the accumulator is a plain product over slots, so walking
  // (b2old, b2new] on top of the finished product for (b1, b2old] gives exactly
  // the product for (b1, b2new]. Mutually exclusive with `resume`, where the
  // seed is already folded into the saved accumulator.
  Words stage2(const Words& x, const Stage2Plan& plan, u32 reportEvery,
               const std::function<bool(u64, u64)>& progress,
               RoeInfo* mulRoeOut = nullptr, bool normalizeDiff = false,
               const Stage2Pos* resume = nullptr, u32 saveEvery = 0,
               const std::function<void(const Stage2Pos&)>& save = {},
               Stage2Pos* stoppedAt = nullptr, const Words* accSeed = nullptr);

  // V_n(base, 1) mod (2^E-1), for an ARBITRARY full-residue base (unlike
  // lucasV, whose `seed` is a small u32 baked in via a dedicated small-constant
  // kernel) and a moderate index n -- tens of bits, not millions: this is only
  // ever called a handful of times per P+1 stage-2 run, to seed the block
  // recurrence (V_D, V_(mFirst*D), V_((mFirst-1)*D)), never in a hot per-bit
  // loop. No checkpoint/resume -- cheap enough to recompute from scratch on
  // every resume, same precedent as stage2()'s own per-run constant C.
  //
  // Same Montgomery-ladder-pair structure as lucasV (V_2k=A^2-2 via squareLL,
  // V_2k+1=A*B-base via modMul+kSubWords), but kSubWords's result here can be
  // squared on a later bit -- a precision regime lucasV's kSubSmall path never
  // exercises (kSubSmall's own comment notes its perturbation is a few units,
  // far below the rounding margin; kSubWords leaves a full extra bit). Every
  // kSubWords result here is carried back to a normal residue via
  // renormalize() before it can be squared. pp1Stage2's own verification
  // found the analogous worry unfounded THERE (see Gpu::renormalize) -- but
  // that measurement covers modMul with one wide operand, not squaring a
  // value against itself, and this ladder has not been separately measured
  // with the precaution removed. At O(log n) <= ~32 steps the cost of keeping
  // it is negligible, so it stays rather than extending an unverified finding
  // to a case it was not measured on.
  //
  // n == 0 returns V_0 == 2 directly. This is NOT a rare edge case:
  // buildStage2Plan's own invariant (w*d/2 <= b1) forces mFirst >= 1, making
  // (mFirst-1)*D == 0 the common case at realistic bounds, not a corner one.
  Words lucasVResidue(const Words& base, u64 n);

  // P+1 stage 2. `y1` is the completed stage-1 residue V_E(seed,1). Finds a
  // factor q of M_p whenever q+1 is B1-smooth apart from a single prime in
  // (B1, B2], by an approach that does not carry over from P-1's: P-1 works
  // with an EXPLICIT ring element (x^n, freely divisible); P+1's group element
  // is never explicit, only its Lucas-V trace V_n(y1,1) is computable in
  // Z/NZ. Writing beta = alpha^E (alpha the unknown "P+1 unit" for a hidden
  // factor q, ord(alpha) | q+1) and W_n := V_n(y1,1), so W_n mod q =
  // beta^n+beta^-n, direct expansion gives, for any a, b:
  //
  //   W_a - W_b = beta^-a * (beta^(a-b) - 1) * (beta^(a+b) - 1)
  //
  // beta^-a is a unit and q is prime, so W_a == W_b (mod q) iff ord(beta)
  // divides (a-b) OR (a+b) -- the SAME "catches either candidate" property
  // P-1 manufactures by squaring (mD)^2-j^2=(mD-j)(mD+j), except V's built-in
  // evenness (V_-n=V_n) gives it for free on LINEAR indices. So, reusing
  // Stage2Plan's (m,j) pairing geometry completely unchanged (it is pure
  // sieve combinatorics, independent of the algebra that tests a slot):
  //
  //   A_m = V_(m*D)(y1,1)   T_j = V_j(y1,1)   diff = A_m - T_j
  //
  // and gcd(prod(diff over open slots), M_p) carries q whenever ord(beta) is
  // exactly the prime assigned to that slot.
  //
  // A_m is carried forward by a SLIDING PAIR (Stage2Pos::a = A_curr,
  // Stage2Pos::s = A_prev), stepped via the 3-term identity
  // V_(a+b)=V_a*V_b-V_(a-b) with a=mD, b=D: A_next = A_curr*V_D - A_prev. The
  // T-table is built by the same identity walked one index at a time
  // (T_(n+1) = y1*T_n - T_(n-1)), needing no per-entry ladder since P+1's
  // indices are linear rather than squared. BOTH of these persist and reuse a
  // kSubWords result across later steps, unlike P-1's single-use diff -- see
  // Gpu::renormalize's comment for the full account, including why that
  // turned out not to matter here (every persisted value is always exactly
  // one kSubWords application wide, never a compounded sum of many, because a
  // modMul separates every pair of successive kSubWords calls) and why the
  // precaution stays anyway.
  //
  // skipRenormalize is for verification ONLY -- never set true in production.
  // It omits the renormalize() calls, so a selftest can measure whether the
  // precaution is load-bearing rather than merely cautious. As it happens it
  // is not, at least up to ~100000 accumulator slots at a realistic exponent
  // (--selftest=pp1stage2): the "correct" (default) and "broken" accumulators
  // come out bit-identical. Kept anyway; see Gpu::renormalize.
  Words pp1Stage2(const Words& y1, const Stage2Plan& plan, u32 reportEvery,
                  const std::function<bool(u64, u64)>& progress,
                  RoeInfo* mulRoeOut = nullptr, bool normalizeDiff = false,
                  const Stage2Pos* resume = nullptr, u32 saveEvery = 0,
                  const std::function<void(const Stage2Pos&)>& save = {},
                  Stage2Pos* stoppedAt = nullptr, const Words* accSeed = nullptr,
                  bool skipRenormalize = false);

  vector<Buffer<Word>> makeBufVector(u32 size);

  void clear(bool isPRP);

private:
  u32 getProofPower(u32 k);
  void doBigLog(u32 k, u64 res, bool checkOK, float secsPerIt, u32 nIters, u32 nErrors);
};

// Compute the size of an FFT/NTT data buffer depending on the FFT/NTT float/prime.  Size is returned in units of sizeof(double).
// Data buffers require extra space for padding.  We can probably tighten up the amount of extra memory allocated.
// The worst case seems to be !INPLACE, MIDDLE=4, PAD_SIZE=512.

#define MID_ADJUST(size,M,pad)                  ((pad == 0 || M != 4) ? (size) : (size) * 5/4)
#define PAD_ADJUST(N,M,inplace,pad)             (inplace ? 3*N/2 : MID_ADJUST(pad == 0 ? N : pad <= 128 ? 9*N/8 : pad <= 256 ? 5*N/4 : 3*N/2, M, pad))
#define FP64_DATA_SIZE(W,M,H,inplace,pad)       PAD_ADJUST(W*M*H*2, M, inplace, pad)
#define FP32_DATA_SIZE(W,M,H,inplace,pad)       PAD_ADJUST(W*M*H*2, M, inplace, pad) * sizeof(float) / sizeof(double)
#define GF31_DATA_SIZE(W,M,H,inplace,pad)       PAD_ADJUST(W*M*H*2, M, inplace, pad) * sizeof(uint) / sizeof(double)
#define GF61_DATA_SIZE(W,M,H,inplace,pad)       PAD_ADJUST(W*M*H*2, M, inplace, pad) * sizeof(ulong) / sizeof(double)
#define TOTAL_DATA_SIZE(fft,W,M,H,inplace,pad)  (int)fft.FFT_FP64 * FP64_DATA_SIZE(W,M,H,inplace,pad) + (int)fft.FFT_FP32 * FP32_DATA_SIZE(W,M,H,inplace,pad) + \
                                                (int)fft.NTT_GF31 * GF31_DATA_SIZE(W,M,H,inplace,pad) + (int)fft.NTT_GF61 * GF61_DATA_SIZE(W,M,H,inplace,pad)
