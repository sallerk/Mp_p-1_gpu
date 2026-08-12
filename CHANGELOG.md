# Changelog

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
