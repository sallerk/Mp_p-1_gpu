# Changelog

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
