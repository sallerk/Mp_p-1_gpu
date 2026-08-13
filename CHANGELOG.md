# Changelog

## 1.3

**P+1 gets its own B1 model.** P+1 used to silently borrow whatever B1 P-1's
Dickman-rho model chose — a different smoothness target (`q+1`, not `q-1`),
optimised for the wrong quantity. `pp1Prob()` is the same rho-integration as
`pm1Prob()`, but removes only 1 bit (the known factor of 2 in
`q+1 = 2(kp+1)`) instead of `log2(exponent)+1`: unlike P-1's `q-1 = 2kp`, `p`
never divides `q+1` (`kp+1 == 1 mod p` for any k), so it isn't a free factor
here. `choosePP1B1()` picks the B1 minimising expected work across the
configured seed list, folding in that each seed pays its own full ladder (one
squaring plus one multiply per bit, no fused-carry shortcut) and its own gcd,
and that only ~50% of seeds land in the residue class P+1's mechanism
actually needs (`1 - 0.5^n` chance across n seeds). B2 and the pairing shape
are still borrowed from P-1's own choice — P+1 has no cost model of its own
for stage 2 yet.

**This is a real behavioural change, not just a bugfix**: P+1's auto-chosen
B1 is now typically much smaller than before (the quantity that has to be
B1-smooth for P+1 is about `p` times larger than P-1's, and Dickman's rho
falls fast, so paying for a bigger B1 buys much less here). A `method = both`
job now runs P-1 and P+1 at two different B1s, shown separately in the
startup summary and in `--bounds`.

**Also fixed**, found by inspecting a real `method = both` run: the banner
said "P-1 factoring" regardless of method; the `M<exponent>, trial-factored
to...` summary line was missing its indent; "estimated success" was P-1's
own model shown even when P-1 wasn't running (a pure `method = pp1` job
never runs P-1 at all); and the P+1 and P-1 sections ran together with no
separator. P-1 and P+1 now get their own labelled `estimated ... success`
lines, and the P-1 section gets a `P-1 attempt` header (matching P+1's own
`P+1 attempt N of M`) with a blank line before it when both methods run.

**`--tune` fixes and visibility**, found by running it for real against a
small (~3.3M) exponent: `tune.txt` could end up silently empty for any
exponent small enough that the smallest available FFT shape already exceeded
the size ceiling `--tune` uses to skip needlessly-oversized candidates --
that ceiling now retries once, bounded to 10x, if the strict pass finds
nothing. Separately, `TAIL_TRIGS32`/`TABMUL_CHAIN32` tested against an FP32
NTT shape that never went through the same broken-baseline check the main
shape gets, so on hardware where that shape fails, every value came back an
identical failure and "Best ..." was picked from noise -- it now gets the
same verify-and-substitute treatment. `--tune` also now prints `[i/N]
Finding best PARAM` progress through its ~20 kernel-option searches and
`[shape i/N]` through its FFT-shape sweep, and every logged timing shows
`FAILED (wrong residue -- not used)` inline instead of a bare sentinel
number, with a `WARNING` line if an entire parameter's search came back all
failures. `-h`'s help text was also corrected to match: two working
`--selftest=` modes (`extend`, `pp1`) were missing from its list, STAGES
never mentioned P+1 at all, and the `noconfig`/`--bench`/`minexp=`/`maxexp=`
descriptions no longer matched current behaviour.

The message `--tune` prints after finishing (`fft = auto` already picks and
verifies the same recommended entry automatically, so pinning it is optional,
not required) now says so instead of reading like an instruction, and its
warning about broken entries no longer blames "the kernel-option search" by
default or points at the wrong file to delete -- and no longer conflates
shapes that are structurally too imprecise for the exponent (an expected,
harmless mismatch) with shapes that actually compute wrong results.

## 1.2

**P+1 stage 2.** P+1 could only find a factor q whose q+1 was entirely
B1-smooth. It can now also catch q where q+1 is B1-smooth apart from one
additional prime in (B1, B2] — the same near-miss stage 2 has always caught
for P-1.

P-1's exact pairing trick doesn't transfer: P-1 works with an explicit ring
element (x^n, freely divisible); P+1's group element is never explicit, only
its Lucas-V trace is computable in Z/NZ. Real implementations (GMP-ECM) handle
this with a much heavier machine (Chebyshev/RLP resultants via FFT) built for
asymptotic speed. This instead uses a simpler construction, derived from first
principles and checked three ways (literature cross-check, an independent
numerical check on a real prime, and an algebraic proof) — see
`1.2/README.md` and the comment above `Gpu::pp1Stage2`. In short: V's built-in
evenness (V_-n = V_n) gives P+1 the same "catches either candidate" property
P-1 manufactures by squaring exponents, but for free, on linear indices — so
**`Stage2Plan`, the sieve and pairing-assignment code, is reused completely
unchanged**. Only the arithmetic underneath a slot is new.

A precision question came up during design — do the block-step and T-table
build need explicit re-carrying, since (unlike P-1's single-use per-slot
difference) they persist state across many steps? Verified empirically
(`--selftest=pp1stage2`, up to ~100000 accumulator slots at a realistic
exponent, with a deliberately-broken control arm) that the answer is no: the
recurrence always interposes an ordinary carrying multiply between successive
uncarried subtracts, so the width never compounds past what P-1's own
single-use case already relies on being safe. The precaution is kept anyway —
cheap, and the failure mode of removing it and being wrong would be a silent
bad answer, not a crash — but the reasoning is written down rather than just
assumed, including the fact that the original worry turned out to rest on an
incomplete mental model of the recurrence.

Shares P-1's B1 and B2/pairing-shape choice (P+1 already had no B1 model of
its own; this is the same category of approximation, not a new one). Per-seed
checkpointing and interrupted-walk resume are supported; B2 extension for P+1
is explicitly out of scope for this version — the checkpoint format has room
for it later without another format bump.

Also fixed: a pure `method = pp1` job's B2 was never actually computed —
`chooseBounds` was only ever told a stage 2 might be wanted when P-1 itself
was active, so a P+1-only run always saw stage 2 as "not worth it" regardless
of `b2` or `stages`. And P+1's progress lines, previously hardcoded to a
literal `[n/3]`, now use the same phase-counting machinery P-1's lines do, so
a stage-2 run reports `[n/5]` correctly instead of always claiming 3 phases.

**New self-test** `--selftest=pp1stage2`: the underlying algebraic identity
checked directly (no GPU); the new Lucas-ladder primitive against an
independent CPU reference (ordinary exponentiation in a quadratic ring,
sharing no code with the GPU ladder); the full walk against the same kind of
reference; the precision question above; interrupt/resume; checkpoint
round-trip and rejection rules; and a real factor of M862907 whose missing
prime is provably outside B1's reach and inside B2's.

## 1.1

**B2 extension.** Raising B2 on an exponent whose stage 2 already finished now
walks only the new range instead of starting over, the way raising B1 already
reused a completed stage 1. Going from B2 = 4M to 60M costs the gap, not the
whole range.

The accumulator is a product over slots, so `acc(b1,b2new] == acc(b1,b2old] *
acc(b2old,b2new]`; the earlier run's pairing shape (D, w) does not have to match
this one's. Consequently:

- the stage-2 checkpoint is no longer deleted when a run finishes — it is
  rewritten as a completed record holding the final accumulator, and the record
  it supersedes is removed once the wider one is safely on disk;
- re-running an exponent whose stage 2 is already complete for that B2 skips
  straight to the gcd;
- checkpoint format v2: adds `complete`, `fromB2` and a res64 of the stage-1
  residue. Every one of them is compared on load, so an accumulator can never be
  resumed or seeded into a computation it does not belong to. Stage-2
  checkpoints written by 1.0 are rejected with a message; stage-1 checkpoints
  are unaffected.

Turn it off with `extend = 0`, which now covers both stages.

**New self-test** `--selftest=b2extend`: the seeded accumulator against exact CPU
arithmetic, a real factor of M86255591 whose missing prime lies in the gap (and
which is provably *not* found before the extension), and the rejection rules for
foreign, stale or unfinished accumulators. `--selftest=stage2plan` additionally
checks that two half-ranges cover exactly the primes one whole range would.

## 1.0

First release. GPU P-1 with both stages, P+1 stage 1, automatic bound selection,
checkpoint/resume, and PrimeNet-format JSON results.

Ships as a prebuilt 64-bit Windows binary — statically linked CRT, OpenCL
resolved from the driver at run time, so there is no toolchain or SDK to install.
