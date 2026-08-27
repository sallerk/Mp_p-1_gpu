# Changelog

## 1.9.2

**Four more ways a worktodo.txt line could produce a wrong reported value,
found by auditing the file against 1.9.1 rather than against 1.9.**

**A line longer than 1023 bytes was silently split into two entries.** Both
the load and the consume loop read with `fgets(buf, 1024, f)`, which returns
the tail of a longer line as though it were the next line. A comment of 2009
bytes whose tail happened to be digits queued **M1000003 -- an exponent that
appeared nowhere in the file** -- and wrote its result to results.txt under
it. Consuming that phantom was worse: the truncated head went back to disk
without the newline it never had, gluing the next real entry onto the end of
a comment and dropping it from the queue with nothing said. The realistic
trigger is a `Pminus1=` line whose `"known_factors"` list runs past the
buffer. Both loops now share one whole-line reader, so their line numbering
and what consume copies through cannot disagree.

**Known factors were never checked against M_p.** They were validated as
decimal digits, then echoed verbatim into results.txt as `"known-factors"` --
a claim about M_p that PrimeNet reads -- and that echo was the *only* thing
the program ever did with them, so a typo had nothing downstream that could
notice. `Pfactor=1,2,1000003,-1,70,1.5,"123456789"` was accepted without
complaint. Each entry is now verified with `2^p == 1 (mod f)`, the same test
already applied to every factor this program reports on its own account.

**A known factor found again was reported as a discovery.** A listed factor's
`k` is B1-smooth by construction -- that is why it is on file -- so it comes
out of the stage-1 gcd of every run whose bounds reach it. That ended the
job as "factor found", skipped the stage 2 the assignment had asked for, and
submitted a status `F` line for a factor PrimeNet has held for years. A real
case: M4444091 with its two smallest factors declared known used to stop at
stage 1 and report both; it now reports them as rediscovered, runs stage 2,
and finds **636358278473**, which the old behavior never reached at all. If
nothing new turns up the result is `NF`, with the known factors still echoed.

**An assigned B1 below 105 aborted the job after stage 1 had finished.**
A `Pminus1=` line always asks for a stage 2, and the smallest pairing shape
(D=210, w=1) cannot walk primes at or below `w*D/2`; `buildStage2Plan` throws
-- from inside stage 2, an hour in, and again on every restart, because the
entry is still queued. The floor is now a named constant (`STAGE2_MIN_B1`,
asserted by the stage-2 self-test to track `STAGE2_D_CANDIDATES[0]`) and is
refused at parse time for an assignment, and before any work for a `b1`
pinned by hand in config.txt.

Three new self-test sections (`--selftest=worktodo` O, P, Q) cover the buffer
boundary in both loops, eight known-factor lists accepted or refused on
divisibility, and the B1 floor. Sections F, G and I had to be fixed to land
them: their known-factor fixtures were written against exponent 86243, which
is a Mersenne **prime** and therefore has no factors at all. Invented values,
invisible while nothing checked them. They now use real factors of M4444091.

## 1.9.1

**1.9's exponent fix was incomplete, and this completes it.** 1.9 guarded
`parseKbnc`, which serves `Pfactor=`/`Pminus1=` lines. A **bare exponent**
takes a different branch and still wrapped: `9993000001` on its own line
became **M1403065409** and ran. The corpus added in 1.9 only ever exercised
assignment lines, which is exactly why the gap survived the release that was
meant to close it.

**An out-of-range B2 silently became a run that did nothing.**
`Pminus1=1,2,1000099,-1,100000,1e18` produced `B1 = 0`, no stage 2,
"estimated P-1 success 0.000%", and completed -- writing a *no factor found*
result for an assignment that had asked for real work. The cause was not the
parser: `chooseBounds` discards every candidate above the machine-word cap
(stage 2 forms `(m*D)^2` as a machine word), the candidate list came back
empty, and `if (cands.empty()) { return best; }` handed the caller a
default-constructed `Bounds` -- b1 = 0. Nothing checked. Reachable from
config.txt's own `b2` as well as from an assignment, so it is fixed in both
places: the cap now lives in one spot (`STAGE2_B2_MAX`, Bounds.h) and is
validated at the parser, while `chooseBounds` and `choosePP1Bounds` throw
instead of returning zeros.

**Numeric fields no longer go through `strtod` where the value is used rather
than checked.** `strtod` accepts far more than the format specifies, and every
one of these was a live silent-wrong-answer path:

- `1000099.7` rounded to a **different exponent**, M1000100, and ran.
- `0x10` became M16; `1e9` became M1000000000.
- `nan` defeated every range test -- NaN compares false against everything, so
  `n < 3`, `n > 4294967295`, `b1 < 2` and `b2 < 2` all missed it -- and
  arrived as **M0**, failing much later with "No FFT config passed
  verification for M0" rather than naming the real problem.
- `b = 1.6` and `b = 2.4` were both accepted as base 2, because the test was
  `u64(b + 0.5) != 2`.

The exponent now takes a strict decimal-integer parse, and the remaining
range tests are written in positive form (`!(x >= lo && x <= hi)`) so NaN
fails them instead of slipping through.

**The exponent must now be prime, which nothing had ever checked.**
`config.txt` has always said so. For composite `p`, every `d | p` gives an
algebraic factor `2^d - 1` of `2^p - 1`, so a run against a mistyped exponent
could "find" a factor and write it to `results.txt` as though it were a
discovery.

**Worktodo parse errors are no longer discarded by `--bounds` and `--tune`.**
The peek that reads the queue for those modes threw the error away one frame
from where it would have helped, so a malformed line surfaced as the far less
useful "needs an exponent in worktodo.txt".

**New self-test section N: 18 inputs that must be refused, each with a
diagnostic.** Every entry is something that was once accepted and quietly
turned into a different job. It lists the bare and assignment shapes side by
side on purpose -- covering only one of them is what let the 1.9 gap through.
A silently dropped line counts as a failure here, not just a wrong one.

## 1.9

**Fixed a silent wrong-answer bug: a worktodo exponent above 2^32-1 wrapped
into a different exponent and ran.** `parseKbnc` bounded the exponent from
below (`n < 3`) but never from above, then did `exponent = u32(n + 0.5)`.
So `Pfactor=1,2,9993000001,-1,96,2` became **M1403065409** -- not an error,
not a warning, just a full P-1 run against a Mersenne number nobody asked
about, with a `results.txt` line to match. Anything above the limit is now
refused by the parser with the exponent quoted back. The limit is
representational rather than a policy choice: `WorktodoEntry::exponent` and
`Config::exponent` are both `u32`, so such an exponent cannot be held at all.
Found while widening the worktodo corpus below, not by a bug report -- which
is the argument for corpus tests over hand-picked fixtures.

**Stage 1 + stage 2 against known factors, as a new self-test section
(`--selftest=pm1`, G3 section D).** The suite had two ends covered and a gap
in the middle: G3 proved stage 1 reaches a factor whose `k` is B1-smooth, and
M7 proved the stage-2 engine matches a CPU reference, but nothing ran *the
whole method at published bounds and got the published factor back*. Eight
vectors now do, spanning 256K- to 8M-word transforms and exponents from 1e6
to 3e8. Every one was verified offline before being written down: the factor
divides `2^p-1`, and `(f-1)/(2p)` is B1-powersmooth apart from exactly one
prime in `(B1, B2]` -- the property that makes stage 2, specifically, the
thing that finds it. Two vectors carry a second factor the same bounds also
reach, so the expected value is the product.

Four of the eight are drawn per run, so the gate stays bounded (~6 minutes)
while every vector gets exercised across repeated runs. The selection seed is
printed and `MP_SELFTEST_SEED` pins it: a random subset you cannot replay is
a bad test, because the run that fails is the one you need to repeat.

One constraint worth recording, since it shapes which bounds are usable at
all: `buildStage2Plan` refuses a shape with `w*D/2 > B1`, so at the floor
shape (D=210, w=1) no vector with B1 below 105 can run a stage 2. Vectors
whose natural bound is smaller are raised to 1000, which is harmless -- a
larger B1 only covers more, and each vector's stage-2 prime stays well above
it.

**A 60-entry `Pfactor=` corpus in the worktodo self-test (section M).**
A PrimeNet-shaped sweep from 1e6 to 1e10 and TF 67 to 96 bits. Section F
already covered the *shape* of a `Pfactor=` line (assignment ID, known
factors); this covers its *range*, which is where the parser was actually
wrong. The five over-range entries are the point of the corpus rather than an
afterthought -- they assert refusal, and would have caught the wraparound
above the day it was written.

## 1.8.1

**A `Pminus1=` assignment's own B1/B2 are now used whenever present, by
default.** 1.8 shipped this gated behind a `bounds_source = auto |
assignment` config key that defaulted to *ignoring* the assignment and
recomputing bounds instead -- which meant a `Pminus1=` line from AutoPrimeNet
silently didn't do what it said unless the user knew to flip a switch. That
key is gone: an assignment's own bounds now always win, falling back to
config.txt's own b1/b2 (auto by default) only for a `Pfactor=` line or a bare
exponent, which have none to honor. This changes what a submitted result
under the default config *means* for anyone using the AutoPrimeNet interop
added in 1.8 -- the reason for the version bump.

**Fixed a real crash: `gcdHalf(x, x)` -- two exactly equal operands -- could
hit a division-by-zero assertion.** Found while auditing the self-test suite
for coverage gaps (see below): `reduceRec` can legitimately drive the second
operand to exactly zero (the `a == b` case reduces to `(x, 0)` in one step),
and the caller's fallback division didn't check for that before dividing by
it. No prior test had ever called `gcdHalf` with two equal large operands --
every existing differential/planted-factor test builds its two operands from
independently random values, which are never equal in practice. `gcdHalf` is
not the default gcd (GMP's `mpz_gcd` is, since 1.8) so this had no
production impact, but it's a real bug in the reversion path, not just a
test gap.

**Selftest audit, on request: found and closed 6 coverage gaps.** Most
notably, the bounds-precedence logic above had zero automated coverage
before this release -- it had already changed direction twice in one
session with nothing but manual smoke tests catching regressions. Pulled it
out into a testable `resolveBounds()` (Worktodo.h/cpp) with its own test
section. Also: a vacuous self-test assertion (`check(X.isZero() ||
X.isZero(), ...)`, an OR of the identical expression) fixed to the real
check it was clearly meant to be; new zero/equal/word-boundary-size gcd
edge cases across all four gcd tiers (this is the test that caught the
`gcdHalf` bug above); a correctness check that `gcd_threads` cannot silently
corrupt results on the default (single-threaded `gcdGmp`) path even though
it no longer has any effect on it; and two more malformed-input fixtures for
the `Pfactor=`/`Pminus1=` worktodo parser (a non-decimal `known_factors`
entry, and a smaller, more typo-shaped ambiguous `sieve_depth`/`B2_start`
value).

**Clarified three stale/unclear config.txt comments**, reported directly by
a user reading the file: `gcd_threads`'s comment still described the old
multi-threaded `gcdHalf` behavior as current, when the default gcd has
ignored it since 1.8; `wait_for_work`/`wait_poll_seconds` never stated what
values were accepted; and the AUTOPRIMENET INTEROP block wrongly called
`Pfactor=` "trial-factoring-flavored work" -- this program does no trial
factoring at all (that's Prime95's separate, unsupported `Factor=` keyword).
Both `Pfactor=` and `Pminus1=` are P-1 work; they differ only in whether
PrimeNet already chose B1/B2 (`Pminus1=`) or the client computes its own
from `how_far_factored`/`tests_saved` (`Pfactor=` -- the same role this
program's own `factored_to`/`bias` already play for `auto`).

**`checkpoint_seconds` was mislabeled the same way, more seriously: it sat
directly under the "CPU THREADS FOR PHASE 3 (the gcd)" heading with no
comment of its own**, making it look like it checkpoints the gcd, which it
cannot -- that phase is one opaque call with no intermediate state to save.
It actually governs stage 1's squaring ladder and stage 2's pairing walk
only. Given its own heading now, saying so. Also lowered its default from
300s to 90s, on request.

## 1.8

**worktodo.txt now accepts `Pfactor=`/`Pminus1=` assignment lines, the shape
AutoPrimeNet (github.com/tdulcet/AutoPrimeNet, successor to gpuowl's
`primenet.py`) and Prime95 itself write.** This program still never talks to
PrimeNet -- no networking of any kind was added -- but it can now sit in the
same directory as AutoPrimeNet and read/write the same two files. A bare
exponent line still works exactly as before; a `Pfactor=`/`Pminus1=` line
additionally carries its own k,b,n,c (validated as a Mersenne number: k=1,
b=2, c=-1), an optional assignment ID, and optional known factors, all
verified against Prime95's own C parser (`commonc.c`'s `parseWorkToDoLine`)
rather than assumed. The assignment ID and known factors are echoed back into
`results.txt` as `"aid"`/`"known-factors"`, matching Prime95's own submission
JSON keys, so AutoPrimeNet's upload step stays consistent with what PrimeNet
expects back.

**A `Pminus1=` assignment's own B1/B2 are used whenever present.** A
`Pfactor=` line or a bare exponent has none to honor, and falls back to
config.txt's own `b1`/`b2` (auto by default) exactly as before. (An earlier
draft of this feature gated this behind a `bounds_source = auto |
assignment` config key, defaulting to ignoring the assignment's own bounds;
that key has been removed in favor of always preferring them -- the whole
point of running AutoPrimeNet alongside this program is to act on what it
fetches.)

A `Pfactor=`'s `how_far_factored`/`Pminus1=`'s `sieve_depth` now feeds
`factored_to`'s existing auto-default (an explicit `factored_to=` in
config.txt still wins over either). `Pfactor=`'s `tests_saved` is parsed and
printed as FYI only -- it is a GIMPS PRP-credit concept with no equivalent in
this program's own cost model (`bias` means "value of a factor vs. a
composite result," not "tests saved"), so no heuristic maps it into `bias`.

**`Pminus1=`'s optional `B2_start` (stage 2 already partly covered
elsewhere) is parsed but not honored.** The existing `extend` reuse path
needs a local checkpoint whose exponent/B1/residue all match -- there is no
local artifact to seed an externally declared `B2_start` with. A job
carrying one now prints a clear warning and walks the full `(B1,B2]` range:
always correct, just not optimal. Out of scope for the same reason: using an
assignment's known factors to work against the cofactor rather than `M_p`
directly -- that is a real numeric-algorithm change, not a parsing one.

**New config.txt key, `wait_for_work = yes | no` (plus `wait_poll_seconds`,
default 5).** Off by default, preserving the existing "empty queue exits
cleanly" behavior a script or a one-shot run expects. On, an empty
worktodo.txt is rechecked on an interval instead of exiting -- for running
unattended alongside AutoPrimeNet, so a queue that runs dry for a moment
doesn't need a relaunch by hand. Ctrl-C during the wait exits the same clean
way an empty queue always has.

**The production gcd(x-1, 2^p-1) call now runs through GMP instead of this
project's own hand-rolled half-GCD.** Measured at the exact scale this call
runs at in a real run (p ~ 82.5M, 1,290,468 limbs): GMP's `mpz_gcd` finished
in 12.5s, single-threaded, against the previous implementation's 133-141s
across up to 20 threads -- an ~11x wall-clock win, over 200x per core. Same
algorithm family underneath (Lehmer plus a subquadratic HGCD, same idea as
this project's own `gcdHalf`); the gap is entirely GMP's hand-tuned assembly
multiply, including an FFT tier this project has no equivalent of. Statically
linked (`x64-windows-static` via vcpkg -- see `build.bat`), so `Mp_p-1_gpu.exe`
stays the single file every prior release was; no DLL is shipped. GMP is
GPLv3-compatible (LGPLv3+/GPLv2+, licensee's choice; see `ATTRIBUTION.md`) --
Mp_p-1_gpu is already GPLv3 in full. One real trade-off: `mpz_gcd` has no
progress callback, so unlike the previous implementation this specific call
can no longer be interrupted once started -- Ctrl-C is still honored
immediately before it starts, so the uninterruptible window shrank from
gcdHalf's ~140s worst case to gcdGmp's ~12s, rather than being eliminated.
The previous implementation (`gcdHalf`, `gcdLehmer`, `gcdEuclid` in `Gcd.cpp`)
is untouched and still fully self-tested as the reversion path if this is
ever pulled.

## 1.7

**Stage 2's `T_j` table is built by a recurrence instead of an exponentiation
per entry.** The table holds `x^(j^2)` for every `j` in the plan -- the values
each accumulator multiply is taken against. It was built with one
`powSmall(x, j*j)` per entry, roughly `2*log2(j)` multiplies each. It is now a
second-difference chain stepping `j` by 2:

```
x^((j+2)^2) = x^(j^2) * x^(4j+4)        x^(4(j+2)+4) = x^(4j+4) * x^8
```

which costs two multiplies per step of 2 in `j` -- one per unit of `j`, and
independent of how many entries are kept. `D` is even, so every `j` coprime to
`D` is odd and the chain lands on all of them.

This changes no result: the accumulator comes out bit for bit identical and the
multiply count is unchanged. It is setup cost that is no longer paid. The win
therefore scales with table size, and table size is capped by GPU memory, so it
is largest where the table is largest. Measured, same machine, same bounds,
same shape, 1.6 against 1.7:

| exponent | shape | `T_j` buffers | stage 2, 1.6 | stage 2, 1.7 | |
|---|---|--:|--:|--:|--:|
| M5378909 | D=2310 w=9 | 2160 | 33s | 24s | -27% |
| M82589933 | D=420 w=5 | 240 | 2m20s | 2m14s | -4.3% |

Both runs matched 1.6's accumulator exactly (`acc res64` `446bd2d4ded89bd7` and
`5f12b67e85c663c3`) on an unchanged multiply count -- the shape a setup-only
change should have. The wavefront exponent gains least because GPU memory caps
its table at 240 entries; the 7-digit one holds 2160.

**The saving does reach the end of a run.** Worth spelling out, because the
overlap 1.6 introduced makes it easy to get backwards: the stage-1 gcd runs
alongside stage 2, so it is tempting to compare the two and conclude the
shorter one is free. That is the wrong comparison. The final gcd needs the
stage-2 accumulator, so it cannot start until stage 2 ends, which puts
stage 1, stage 2 and the final gcd in series. A full 1.7 run of M82589933 at
B2 = 2,000,000:

```
setup       24s   |  timing FFT candidates (21.9s), init, launch
stage 1   1m33s   |
stage 2   2m07s   |  stage-1 gcd (3m18s) runs here, and now OUTLASTS
final gcd 3m03s   |  stage 2 -- see below
                     24s + 1m33 + 2m07 + 3m03 = 7m07s, as observed
```

`setup` is listed because otherwise the phases do not sum to the run: 1.7
spends ~22s timing FFT candidates before it starts, against ~8s in 1.6, which
verifies only the first candidate that works. The stage-1 gcd's 3m18s is *not*
a term -- it overlaps stage 2 and the final gcd, which together outlast it.

So stage 2 is on the critical path and a second saved there is a second saved
overall.

What the overlap buys is the *stage-1* gcd -- but by the end of 1.7 it no
longer buys all of it. Under 1.6 that gcd and stage 2 take about the same time
(2m57s against 2m47s) and end together, leaving the final gcd the CPU to
itself. 1.7's stage 2 finishes a minute earlier while the stage-1 gcd still
needs ~3m, so the final gcd now overlaps it by around 70s and the two contend
for the same worker threads. That is why 1.7's final gcd measures *slower*
than 1.6's on identical gcd code. Making stage 2 faster moved the bottleneck,
which is worth knowing before optimising it further: the next second saved in
stage 2 is no longer a full second saved overall.

Verified by a new differential check in `--selftest=stage2`: build the table
both ways and compare every entry on the GPU. Direct exponentiation is what
shipped through 1.6, so this pins the new code against a version with real
production mileage rather than against a reference written alongside it.

**A tune.txt built for one exponent could pick a transform 4x too big for
another, and did.** `--tune` is documented as a once-per-GPU step, so its
output routinely gets used on exponents it was never measured at. Every
tune.txt entry passes `maxExp` for any *smaller* exponent, and entries are
offered ahead of the untuned fallback in the order their cost was measured --
at the exponent they were tuned for, which says nothing about this one. Tune
at ~2e7, then run M5378909, and the 1048576-word transform (5.13 bits/word) is
offered before the 262144-word one that actually fits (20.52 bits/word). Both
verify, so the first wins: **416 us/it against 90**, a four-times-too-large
transform for 4.6x the time, silently, on the workflow the README recommends.

`minBpw` already rejected a transform too *small* to be correct -- the same
oversizing in its extreme form, where the engine refuses to run at all. This
adds the missing upper bound: candidates more than 2x the smallest usable
transform are demoted behind everything of a sane size. Within that 2x band
the tuned order still rules, because that is the range where shape tuning
genuinely decides. Demoted rather than dropped, so an oversized transform that
works still beats failing outright.

Verified both directions: M5378909 with a tune.txt from a larger exponent now
picks 262144 words instead of 1048576, and M82589933 with its own tune.txt
still picks the tuned 1:512:8:256:202 exactly as before.

**A run now times its transform candidates instead of taking the first one
that works.** Candidate order came from
tune.txt, whose costs were measured at whatever exponent the tune targeted; the
oversizing fix above bounds how wrong that order can be, but inside the
remaining 2x band it is still a guess -- and the table below shows a 4.9x
spread inside that band. Selection took the first candidate that passed its
correctness check and stopped looking.

It now times two or three candidates and keeps the fastest. The budget gates
*starting* a measurement, not finishing one, and a candidate cannot be cut
short, so a run overshoots rather than stopping dead -- see the ordering entry
below for what that settles at and why rejected transforms can push it further.

There is no exponent floor on this. There was at first -- 8 digits, on the
reasoning that a shorter job cannot afford the wait -- but what scales with the
exponent is the cost of looking, not the value of it. A candidate at 262144
words builds and times in about two seconds against ten at 2097152, so the same
comparison costs 7s at M786433 and 23s at M99700031. The floor also turned out
to be actively harmful; see the ordering entry.

How wide the field actually is, measured on M13466917 with no tune.txt at all,
every candidate 524288 or 1048576 words and so all inside the size band the
oversizing fix allows:

| shape | us/it |
|---|--:|
| 2:256:4:256:101 | 222 |
| 4:256:4:256:101 | 244 |
| 1:256:4:256:101 | 333 |
| 256:6:256:101 | 668 |
| 256:7:256:101 | 716 |
| 512:3:256:101 | 740 |
| 256:3:512:101 | 831 |
| 1K:2:256:101 | 1094 |

Every one of these passes its correctness check, so under the old rule the
first in candidate order won outright -- a 4.9x spread decided by an ordering
that never measured this exponent.

**The candidate ORDER now comes from the bpw table, and the timing is long
enough to trust.** The search above only helps if the candidates worth timing
are the ones it reaches, and if the numbers it takes are real. Neither held.

Candidates arrived in `FFTShape::allShapes()` catalog order, which has nothing
to do with speed, and only the first eight are ever timed. Measured cold on an
RTX 3070 with no tune.txt, that cost real time two different ways:

- The fast shape sat too late to reach. M55000013 chose 687 us/it on its first
  run, 650 on the second, 630 on the third, converging only as the verdict
  cache filled in. Four of twelve exponents behaved this way, all between 52M
  and 83M.
- Worse, M20000003 settled on `4:256:4:256` having already timed
  `1:256:4:256` and ranked it slower. Measured over 8000 iterations the two are
  306 and 169 us/it. It converged on a shape **81% slower** than one in its own
  candidate list, and no amount of re-running fixed it.

The order now ranks on precision headroom, `maxBpw - E/size`, smallest first,
after sorting by transform size. Every family trades speed for bits per word,
so the cheapest arithmetic that still clears the exponent is the one to try
first: here M31^2*M61 (~32 bpw) beats M31*M61 (~40) beats M31^2*M31*M61 (~48),
and plain FP64 is 4.3x off the pace at M65000011 because consumer nVidia runs
FP64 at 1/64 rate. Headroom expresses that without naming any of it, which is
the point -- given a card with fast FP64 the same rule reaches for FP64 exactly
where it fits. Size has to lead, or a tighter-fitting shape one size up
outranks the right answer: at M20000003 that is the M61 shape at 1048576 words
(216 us/it) jumping the 524288-word M31*M61 (169).

Shapes also appear once per usable variant, and for every NTT family the bpw
table is variant-independent, so `:101` and `:202` of one shape tie on both
keys and land adjacent. They also run at the same speed -- 641.4 vs 641.8,
887.4 vs 886.7, 937.2 vs 936.0 us/it, inside 0.2% each time -- so every shape
now gets one variant tried before any shape gets a second.

The measurement was the other half. `verifyOne` timed 400 iterations, about
0.3s of GPU work, which is not a measurement of the transform so much as of
where the clocks happened to be. At M65000011 it put `2:512:8:256:101` at 891
us/it and `2:256:16:256:101` at 667, where careful runs give 609 and 621 --
wrong by 46%, and in the wrong order, which is the part that costs time. It now
runs 3000 iterations. That is affordable only because the ordering fixed what
gets timed: the budget buys accuracy on a couple of good candidates instead of
noise on many. Most of a candidate's 6-12s is building the transform, not
timing it, so candidates rather than iterations are what the budget really
buys, and two comparisons proved sufficient at every exponent tested.

Rejections no longer starve the comparison, either. A candidate that fails its
correctness check costs as much wall time as one that is timed, and they come
in runs -- M51900019 rejects five straight, 36 seconds of them, which used to
exhaust the budget one measurement after the first shape that worked, leaving
the search to "choose" from a field of one. The budget may now only stop the
search once two candidates have actually been compared.

Cold-start first-run pick, against exhaustive interleaved head-to-heads of
every safe shape at each exponent:

| | before | after |
|---|--:|--:|
| exactly optimal | 8 of 13 | **11 of 13** |
| within 2% | 8 of 13 | **13 of 13** |
| worst case | +81% | **+2.0%** |

The rule pays most where the families are closest to each other in precision
and furthest apart in speed. At M5378909 the old catalog order took
`2:256:2:256` at 148.7 us/it where M61's `3:256:2:256` does 52.2 -- the
ordering makes that exponent **2.85x faster**.

It also has a failure mode, found by running the rule against the full range
rather than the wavefront: it assumes the cheapest arithmetic is the one with
the lowest bits-per-word ceiling, and FP64 is the exception -- lowest ceiling
of all at ~19 bpw, and 4.3x off the pace here, because consumer nVidia runs
FP64 at 1/64 rate. FP64 only becomes usable below ~19 bpw, which at the
smallest transform means an exponent under about 5.1e6 -- precisely the region
the 8-digit floor left unmeasured. At M786433 that ordered FP64 first and took
it unmeasured: 234 us/it where M61 does 148.8, 57% given away. Removing the
floor fixes it outright, because two candidates get compared and FP64 loses on
the numbers rather than on an assumption. On a card with fast FP64 it would win
there, which is the point of measuring instead of asserting.

The residue is a limit worth knowing: run-to-run GPU clock variation is around
10%, so a 3000-iteration measurement cannot reliably separate shapes closer
than about 5%. At M65000011 the top two are 2% apart and it takes either. What
the search now reliably avoids is the large mistake, which is what mattered.

**The tune suggestion is now only printed where a tune could pay for itself.**
A run with no usable tune.txt used to close with "Consider: Mp_p-1_gpu.exe
--tune quick=10,minexp=...,maxexp=..." at every exponent. --tune runs for
15-90 minutes; a 6- or 7-digit job finishes in well under one, so at that size
the advice was recommending an hour of setup to speed up something already
over. It is suppressed below 8 digits. The factual line above it, that no
tune.txt entry covers this exponent, stays -- it explains why the candidates
being timed are untuned ones.

Selection costs 7-23s across the range from M786433 to M120000007, scaling
with transform size. The exception is M51900019, where four transforms fail
their correctness check before a working one appears and push it to 50s -- the
case where spending the time is most clearly right.

**`fft-verified.txt` gained a field.** Each line now ends with the measurement
scale (`q5`). Costs taken over different iteration counts are not comparable,
and mixing them is the exact mistake the change above exists to stop, so an
entry written by an older build is kept for its verdict and its cost dropped.
Those exponents re-time once.

The budget covers only candidates that must actually be measured. `fft-verified.txt`
already recorded a cost next to each verdict; those are now read back, so a
repeat run of the same exponent on the same GPU and driver re-picks the same
winner instantly. Delete the file to force a fresh search. A cost is read as
optional -- a missing or zero one means "unknown" and is re-measured rather
than treated as free, which would otherwise win every comparison.

**`--tune`'s kernel options now record the exponent they were measured for, and
are ignored for any other.** `Mp_p-1_gpu-tune-config.txt` holds the winners of
the tune's `-use` option search -- `MODM31`, `TABMUL_CHAIN61`, `IN_WG` and
about a dozen more. The loader applied every line in it to every job it ever
ran. It had no alternative: the file recorded no exponent, no shape and no
transform size, so there was nothing to check a job against even in principle.

These settings are not neutral when misapplied. A 15-option file replayed
against a transform it was not measured on cost about 32% per squaring -- worse
than running with no tuning at all. The file is also opened for *append*, so
tuning twice leaves two blocks in it and the later keys win by parse order,
which has nothing to do with which one fits.

`--tune` now writes a `# tuned-for exponents <lo>-<hi>` line above each block,
and only a block whose range covers the exponent about to run is applied;
anything else leaves the stock defaults in place, and says so. A bare `--tune`
targets the worktodo exponent, so that range is normally a single exponent.
Blocks written before this -- carrying no range -- are ignored and reported,
rather than being trusted by default.

The decision also moved to where the exponent is known. Options used to be read
once at startup, before the worktodo loop, so every entry in a queue ran under
one option set no matter how far apart the exponents were -- a queue holding a
5-million and an 80-million exponent got the same kernels for both. They are now
re-decided per entry. `--tune`, `--selftest` and `--bench` are
unaffected in the sense that matters: a tune still measures from the stock
baseline, a selftest now tests the engine rather than one tuning, and `--bench`
-- which does have exactly one exponent -- applies a tuning measured for it.

Tagging the block keeps a tune-config off exponents it was not measured for,
but on its own it says nothing about whether the tune was right for the
exponent it *does* claim -- and it was not, because the option search measured
against a shape hard-coded in `tune.cpp` rather than the one that exponent
selects. That is the next entry.

**`--tune` now picks the transform first and trains the kernel options on it.**
It ran the other way round: the option search measured against a shape
hard-coded in `tune.cpp` -- `FFT64 512:15:512`, or `FFT3161 512:8:512` for the
NTT side -- at that shape's own top exponent, whatever the tune was aimed at.
So the winners of ~50 measurements described a 7.5M-word transform at around
150M, and were written to a file a run then applies to whatever transform its
own exponent selects. That mismatch is what cost ~32% per squaring, and the
exponent tag added above only kept it off *other* exponents; it did nothing
about the tune being wrong for its own.

`--tune` now calls `chooseVerifiedFFT`, the same selector a real run uses, at
the exponent the tune targets, and the option search trains on whatever comes
back. The hard-coded shapes are gone, and so is the FP64-versus-NTT probe that
timed two fixed shapes at a fixed `141000001` to guess which arithmetic this
GPU prefers -- picking the transform for the target exponent answers that by
measuring the shapes that are actually candidates.

Which options get searched now follows from that shape's own arithmetic. An
option like `MODM31` used to be measured against a substitute shape that had
GF31 even when the target transform had none, and the winner still went into a
file that applies it to every job. Options for arithmetic the chosen transform
does not contain are simply not searched: tuning M5378909 runs 14 steps rather
than 22, because its transform is a pure M61 NTT with no GF31, FP32 or FP64
work in it at all.

**This exposed a second bug, and it is the more interesting one.** `quick` is an
iteration count, and against a 7.5M-word transform even its shortest setting --
400 iterations -- is nearly three seconds of work. Pointed at the transform an
exponent actually uses, that assumption collapses: at 262144 words an iteration
takes about 50us, so 400 iterations is **20 milliseconds** of timing, which
measures launch overhead and clock state rather than the kernel. One such run
reported `3:256:2:256:101` at 57.3, 173.0, 220.6, 372.8, 419.5, 446.2 and 458.7
us/it across adjacent option values whose true cost is about 52 -- then wrote
the noise out as winners. The resulting config made M5378909 **3x slower** than
no tuning at all.

So the measurement scale now comes from the transform rather than a constant:
`--tune` probes the chosen shape once and uses the shortest run whose timed
window still clears 1.5 seconds. A `quick=` argument is honoured as a ceiling
on speed only -- it can ask for more accuracy than that, never less. On the
same shape the spread across those option values falls from 8x to 5%.

Measured after, config present on one round and parked on the next so drift
cannot land on one arm:

| exponent | transform chosen | tuned vs stock |
|---|---|--:|
| M5378909 | `3:256:2:256:101` | 51.7 vs 51.6 us/it -- neutral |
| M82589933 | `1:512:8:256:101` | 852 vs 839 us/it -- neutral |

Both neutral, and that is the finding: **on this GPU these kernel options do not
produce a win this measurement can see.** Every value in a block lands inside a
few percent of every other, which is why the previous behaviour -- writing the
best of them down as an active setting -- was recording coin tosses. The
transforms chosen are the ones an exhaustive head-to-head of every safe shape
confirms fastest, so what the tune contributes here is the shape, not the
options.

One interim measurement said M82589933 came out 2.4% faster tuned, consistently
across four rounds. It did not replicate: repeated later under sustained load,
with the same config and the same shape, the same comparison gives 3 rounds of 6
ahead and a mean of 0.45%, scatter +/-4%. The absolute level had moved from
691-711 to 800-893 us/it between the two, the GPU having warmed up. Recorded
because it is a good illustration of the trap this whole entry is about: four
consistent rounds looked like a result and were a thermal state.

**An option now has to beat the current setting by more than the measurement
noise to be recommended.** `configsUpdate`'s margins were constants -- 0.003 at
most call sites, 0.000 at six of them -- so anything 1% ahead was written as an
active `-use` line, while run-to-run variation on one shape is around 2-3%. Half
those lines were coin tosses.

`--tune` now measures the floor instead of assuming it: after fixing the
transform and the iteration count it times that exact configuration four times
over and takes the spread, which on `1:512:8:256:101` came to 3.4% (821.2 to
848.9 us/it). Each call site keeps its own literal as a lower bound, so the
threshold is `max(what the author wanted, what the hardware can resolve)`.

The estimate is taken the way the option search itself measures -- repeat runs
inside the one process -- because that is the noise those measurements are
actually subject to, and it re-measures per tune, so a quiet machine admits
more settings than a loaded one.

Nothing is discarded by this. A sub-threshold win still gets written, into the
commented "slightly faster" block that already says it wants timing over a
longer duration before being trusted; the change only moves the line between
recommended and suggested to where the measurement can support it. On
M82589933 that moved five of twelve settings across.

**Three things `--tune` got wrong once it started choosing its own transform,
all reported from a real run.**

*It advised you to run a tune, while running one.* Calling the production FFT
selector brought its "no tune.txt entry covers M... / Consider: --tune ..."
notice along with it, which is sound advice to a job and nonsense to a tune.
The selector now takes a flag saying who is calling.

*It compared three candidates and gave up.* The selector stops early on purpose
-- a job is waiting on it, so it buys a good-enough answer in about twenty
seconds. A tune has been asked for exactly this measurement and has minutes to
spend, so reporting "3 compared, 9 left untimed" is not a tune, it is a slower
guess. For a tune the budget and the comparison cap are both lifted and the
candidate ceiling goes from 8 to 24, so every sane-sized shape gets timed.

*Ctrl-C did not stop it.* `Gpu`'s stop checks go through
`Signal::stopRequested()`, which reports a SIGINT handler this program never
installs -- `Signal`'s constructor is what installs one, and nothing anywhere
constructs a `Signal`. So every one of those checks has always been dead. It
went unnoticed while `timePRP` was only reached for a few seconds of
verification; `--tune` spends its entire run inside it. What does work is
`gInterrupted`, set by the console control handler `main` installs, and which
the P-1 job path has used all along. `timePRP` now honours it, `verifyOne`
rethrows rather than swallowing a stop as "that candidate is broken" (which
would have cost one Ctrl-C per candidate), and `main` reports an interrupted
tune as interrupted rather than as a failure, keeping whatever was already
written.

**Not done: a "V-trace"/Pair95 stage 2.**
[PrMers](https://github.com/cherubrock-seb/PrMers) runs its P-1 stage 2 on the
Lucas trace `V_n = H^n + H^-n`, where `V_(kD) - V_j` covers `kD-j` and `kD+j`
in one scalar, with a "Pair95" extension adding irregular offsets `unit`,
`unit+D`, `unit+3D`, `unit+7D`, ... This was scoped for 1.7 and rejected on the
numbers:

- **The pairing is the same trade at the same price.** Pair95's `L` levels give
  each prime `L` candidate partners for `L*phi(D)/2` table entries; the
  existing `stage2_w` window gives `w` candidates for `w*phi(D)/2` entries.
  Same candidates per unit of memory. `--selftest=stage2plan` shows the rate is
  a function of candidates per prime and barely of `D` at all (1.15-1.21 at
  w=1, 1.47-1.54 at w=5, 1.56 at w=9), and that it flattens hard -- w=1 to w=3
  buys 0.2 primes per multiply, w=7 to w=9 buys 0.03. There is no pairing
  headroom here to go and get.
- **The port would cost more than it saves.** V-trace for P-1 needs
  `V_1 = x + x^-1 mod M_p`, an extended gcd on a p-bit number. The plain gcd
  already runs for minutes at a wavefront exponent, so the inverse alone would
  exceed the whole stage 2 it was meant to speed up. P+1's stage 2 here is
  already a Lucas-trace walk, because a P+1 group element is only ever
  available as a trace -- there the construction is free and natural.

What did transfer is the part that is not about pairing at all: a linear index
lets the baby table be built incrementally. `x^(j^2)` is quadratic in `j`, so it
takes a second difference rather than a first -- and needs no inverse. That is
the change above, and it is where 1.7's speedup comes from.

## 1.6

**Resuming an interrupted stage 1 reprocessed one exponent bit, corrupting the
residue.** `powBase3`, `lucasV`, and `powResidue` all restarted a resumed walk
at `startBit = resumeBit + 1`. `resumeBit` is saved AFTER that bit's `square()`
already ran -- it is the bit the interrupted run *just finished*, already
folded into the saved residue, not "the next bit to process" the `+ 1`
assumed. Resuming past it reprocessed that same bit a second time, silently
producing the wrong `x` for the rest of stage 1 -- and therefore the wrong
`gcd(x-1, M_p)`, meaning a resumed run could miss a real factor without any
error or warning. This is not new in 1.6: the mechanism is unchanged since it
was written, so any interrupted-and-resumed stage-1 walk in any released
version -- Ctrl-C, a crash, a reboot, anything that left a partial
checkpoint and got restarted -- has been affected. Also ported to 1.5, whose
released binary has the identical bug; see that section below.

Found by writing the first real test for this exact scenario -- interrupt a
real ladder mid-walk, resume it, and check the result against an
uninterrupted run bit for bit. The existing coverage never exercised this:
the B1-extension tests only resume a *completed* stage 1 onto a higher B1,
which never touches this code path at all. Fixed by starting the resumed
loop at `resumeBit` itself.

**Stage 1 now checkpoints immediately on interrupt, instead of only
periodically.** Stage 2 already captured its exact position the moment
Ctrl-C was noticed; stage 1's squaring and Lucas ladders did not -- they only
saved on a `checkpoint_seconds` timer or on completion, so an interrupt
between two scheduled saves discarded up to a full checkpoint interval of GPU
work (redone, not lost, but wasted). Fixed by forcing a full carry on the
iteration that precedes every progress check when checkpointing is on (P+1's
ladder needed no such forcing -- its state is always fully carried), so an
interrupt always has a valid residue to save right away. Writing the first
real test for this -- which had to actually interrupt and resume a walk
rather than just extend a completed one -- is what surfaced the resume bug
above.

**Ctrl-C did not stop a run during the gcd phases.** `gGcdProgress` was
void-returning, so nothing inside `gcdHalf` could ever be told to stop. The
gcd is the longest single phase of a run (minutes at production scale) and
every GPU phase already had an equivalent bool-returning progress callback
that broke its loop on `false` -- the gcd was the one place Ctrl-C was
silently ignored until it finished on its own, at which point the job looked
complete and (with a stage 2 configured) was removed from `worktodo.txt`
without ever writing a result. Found live: a `worktodo.txt` that lost its
queued exponent with no matching `results.txt` entry after Ctrl-C during a
gcd. `gGcdProgress` now returns bool; returning false throws a `GcdAborted`
that unwinds to whichever caller installed the hook, which sets that result's
existing (previously dead for this case) `interrupted` field. Also ported to
1.5, whose released binary has the same bug -- see that section below.

**The stage-1 gcd now runs while the GPU does stage 2.** Not a faster gcd --
the same gcd, moved off the critical path. It was CPU-only work with the card
sitting idle, and at M82589933 that is ~184s of a 613s run during which the
GPU does nothing at all.

Stage 2 never needed the gcd's *answer*. It continues from the stage-1
residue, which is ready the moment stage 1 ends; it only needs to know
whether a factor turned up, because a factor in hand makes stage 2 pointless.
So the driver starts the gcd on another thread and goes straight into stage 2.
Measured on the same job: **613s -> 470s**, 137s saved, with the stage-2
accumulator residue bit-identical to 1.4 and 1.5.

It is a speculation, and the odds are on the record -- the bounds model prints
its own stage-1 success estimate, typically a few percent. In the ~97% of jobs
where the gcd finds nothing the overlap is free. In the rest the stage-2 work
done so far is discarded, but a factor was found and the job ends successfully
anyway: losing that race is the good outcome. The saving is bounded by stage
2's own duration, so it is largest exactly where it matters -- big exponents,
where both phases are long.

Nothing is reported until the gcd lands. Announcing "no factor found with B1"
while it is still running would be a guess rather than a result, so the
stage-1 report moved after the join, and the overlapped gcd runs silent so its
progress line and stage 2's cannot interleave on one terminal.

**A latent data race, made reachable by the above, is fixed.** `gGcdProgress`
was a plain global: stage 2's gcd installs a callback capturing its own locals
by reference, and with two gcds now genuinely concurrent the other one could
invoke it from the wrong thread. That was safe for exactly as long as only one
gcd could be running, which was true until this release. The hook and its
bookkeeping counters are `thread_local` now.

**The schoolbook base case got its threshold corrected and its inner loop
tightened.** The post-1.5 profile put ~37% of ALL samples inside one ~62-byte
span -- `mulSchoolbook`'s inner multiply-accumulate, reached ~740M times as
the base case of every other tier's recursion.

`KARATSUBA_LIMBS` had been 40 since before Toom-3, Toom-4 and the
thread-local allocator existed, and re-measuring found the real crossover is
48: at 40 limbs schoolbook is still the faster of the two (1.579us against
Karatsuba's 1.765us), and Karatsuba only takes the lead at 48 (2.500us
against 2.302us).

The inner loop itself was read-modify-writing the accumulator word twice per
iteration -- once to add the product's low half, once to add the incoming
carry -- two dependent round-trips to the same address. Folding the carry
into the low half in registers first leaves one. Identical arithmetic, half
the accumulator traffic: production gcd **143.0s -> 139.3s** (two runs
agreeing within 0.25s), aggregate multiply CPU-time 325s -> 314s.

## 1.5

**The gcd is about 2.4x faster.** Nothing about running the program changed --
no new options, no config keys, no file formats. This is entirely about the
phase that prints `gcd CPU` and leaves the GPU idle, which at 100M-digit
exponents was taking longer than the GPU work it followed. At the production
size used for measurement throughout (M82589933, a 1,290,468-limb gcd), one
`gcdHalf` call went from **338s to 143s**.

Four changes got there, and the order matters, because two of them only paid
off because of the others.

**Toom-Cook-3's sub-products now run in parallel.** `Gcd.cpp` already ran its
matrix sub-multiplies across cores; `mulToom3`'s own five independent
evaluation-point products did not. They do now, through a thread-budget
dispatcher shared between `Gcd.cpp` and `BigInt.cpp` (new `Parallel.h`) --
shared rather than per-file because these call sites nest, and two independent
budgets would oversubscribe the machine. 338s -> 226s.

**The gcd was allocator-bound, which was not what anyone assumed.** Every
`Nat` operation allocates a backing vector, and at production scale the gcd
does that hundreds of millions of times. A CPU sampling profile (~1.8M
samples, aggregated by module) put **~45% of all process CPU time inside the
Windows heap**, against ~53% in the program's own arithmetic. `Nat`'s
allocator is now a thread-local, lock-free free list (new `Pool.h`), and
`Parallel.cpp`'s workers became long-lived so those per-thread pools survive
between tasks instead of being born cold. 226-233s -> **148.8s**. A
re-profile of the result put the heap at **3.7%**, down from ~45%.

**Toom-Cook-4 was added, and it is switched on -- but only after the
allocator fix.** This is the part worth reading. Toom-4 was implemented,
verified against Toom-3 by differential testing, and then deliberately left
*dormant*, because switching it on made the real gcd slower (250s+ against a
229.5s baseline) even though every isolated single-multiply benchmark said it
was 1.2-1.4x faster. The profile explained the contradiction: Toom-4 trades
arithmetic for temporaries, so it delivered exactly what it promised --
12-13% fewer arithmetic samples -- and lost anyway, because in an
allocator-bound workload the temporaries cost more than the arithmetic they
saved. With the thread-local pool in place the trade turns favourable and the
advantage finally shows up as wall time: 148.8s -> **143.0s**, aggregate
multiply CPU-time 373s -> 325s.

**Both Toom tiers stopped recomputing shared partial sums.** `mulToom4`'s
`p1`/`pm1` and `p2`/`pm2` are plus/minus pairs built from the same two partial
sums, each computed independently twice; `mulToom3`'s `p1`/`pm1` had the same
pattern once. Sharing them cut combination-phase time 7-11%.

Two things were tried and **reverted**, which is worth recording because the
first one was right about the target and wrong about everything else. A
pooled allocator using one *global* free list with a mutex per size bucket
was 3.4x *slower* -- a contended global lock is strictly worse than the
per-thread caching the Windows heap already does, across ~19 threads. Same
idea, same target, wrong sharing model; the thread-local version above is
what that attempt should have been. Toom-Cook-6/8 were scoped and declined:
the Toom-3 -> Toom-4 step returned only 3.9% of wall time for a 13%
arithmetic cut, and Toom-6/8 need 11 and 15 evaluation points against
Toom-4's 7, so the per-call cost moves the wrong way faster than the
asymptotics improve.

No behavior outside `BigInt`/`Gcd` changed. The full selftest regression
(`gcd`, `exponent`, `stage2plan`, `bounds`, `worktodo`, `engine`, `pm1`,
`extend`, `pp1`, `stage2`, `b2extend`, `pp1stage2`) passes unchanged, with
Toom-4's own differential tests added to gate G2 (1986 checks, up from 1914).

**Post-release fix (source only, same folder):** the released 1.5 binary does
not honor Ctrl-C during the gcd phases -- see the 1.6 entry above for the
bug and fix, which also applies here since `1.5/src` was patched in place.

**Post-release fix (source only, same folder):** the released 1.5 binary
reprocesses one exponent bit when resuming an interrupted stage 1, corrupting
the residue -- see the 1.6 entry above for the bug and fix (`powBase3`,
`lucasV`, `powResidue` in `Gpu.cpp`), which also applies here. This affects
any run that was interrupted mid-stage-1 (Ctrl-C, a crash, a reboot) and then
resumed from the checkpoint; a run that completed in one sitting, or was
interrupted only during stage 2 or the gcd, is not affected.

## 1.4

**Exponents now come from `worktodo.txt`, and it can hold a queue.** Every
prior version took exactly one exponent per run, from `exponent = ` in
`config.txt`. That key is gone: exponents live in `worktodo.txt` instead
(`worktodo_file` in `config.txt`, mirroring `results_file`'s own pattern),
one per line, `#` comments and blank lines allowed same as `config.txt`. A
run works the queue's first line, appends the result to `results.txt`,
removes that line, and moves to the next -- all in one process, no relaunch
needed between exponents. An old `config.txt` with an `exponent = ` line
still in it fails with a clear migration message rather than being silently
ignored.

Interruption and failure are handled the way the rest of the program already
handles them, not with new retry logic: a Ctrl-C mid-job leaves that entry's
line in place, so relaunching resumes it from its checkpoint (unchanged
behavior) and then continues down the queue from there. A job that fails
(an actual error, not "no factor found" -- e.g. an exponent too small for
any FFT shape) stops the whole queue at that entry rather than skipping past
it, the same thing that already happened to a single-exponent run failing;
the bad line stays put rather than being silently dropped or retried
forever. `--bounds` and `--tune` read the queue's first entry, without
consuming it, to scope themselves the way they used to scope to `config.txt`'s
`exponent`.

The job driver itself was restructured to make this possible: the
FFT-selection-through-stage-2 logic that used to be the entire body of
`runMain` after device bring-up is now its own function, `runOneJob`, called
once per queue entry from a loop around it. The extraction was mechanical --
no logic changed, only where it lives -- which is what let a per-entry
failure keep working the same way it already did (an uncaught exception
still reaches `runMain`'s own top-level `catch`, stopping the run) without
writing any new error-handling code.

New `--selftest=worktodo` gate (no GPU): round-trip parsing, exact-entry
removal, malformed-line and missing-file handling, and specifically a file
whose last line has no trailing newline -- routine for a hand-edited file,
and the reason the new `Worktodo.cpp` reads the file itself rather than
reusing `File`'s existing line iterator (used by `tune.txt`'s reader and by
a vendored, previously-unused `deleteLine()` helper in `fs.cpp` that looks
like it was meant for exactly this): that iterator throws if a line does not
end in `\n`, which is correct for a machine-written file but wrong for one a
human edits by hand.

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
actually needs (`1 - 0.5^n` chance across n seeds).

**P+1 now gets its own B2 too.** `choosePP1Bounds()` (renamed from
`choosePP1B1`) chooses B1 and B2 jointly, the same tolerance-banded scan
`chooseBounds` already runs for P-1, instead of taking P-1's B2 as fixed
input. That coupling had a real effect: because P+1 had to walk a stage-2
range sized for *P-1's* bounds, raising `bounds_tolerance` could push P+1's
B1 up mainly to shrink an expensive borrowed-B2 tail it never asked for —
not because P+1's own economics justified it (measured case: tolerance 0.3
recommended P+1 spend 116.6% of a PRP test's worth of work for a 0.288%
chance, a far worse trade than P-1's own 75.1%/4.25% at the same tolerance).
The stage-2 *cost* constants (`mulsPerPrime`, `s2MulCost`) are still reused
from P-1's own measurement rather than remeasured — `Gpu::pp1Stage2` reuses
P-1's exact pairing geometry and `modMul`-based accumulator recurrence, so
the same per-multiply cost applies regardless of which method is walking it
(now checked live in `--selftest=pp1stage2`, mirroring the check M6b already
does for P-1's own ratio). The stage-2 pairing *shape* (`stage2_d`,
`stage2_w`, T-table sizing) also stays shared between both methods — a
GPU-memory budget decision, not a cost-model gap.

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

**Runs now leave a log file.** `log()` (used throughout `--tune` and, as of
this release, every job's outcome) has always written to stdout only:
`initLog()` — the call that opens its file half — was declared and defined
but never actually invoked anywhere in this program, upstream included, so
nothing was ever recorded once a console window closed. `initLog()` is now
called at startup, appending to `Mp_p-1_gpu.log`, and every failure/outcome
message (factor found, no factor, interrupted, and the top-level `FAILED:
...` a crash exits with) goes through `log()` rather than a plain `printf()`
so it lands there too — routine per-iteration progress stays console-only.
The file is explicitly flushed after every write, unlike the buffered
default a `FILE*` gets from `fopen`: the whole point of the file is
surviving a crash or a forced kill, which an unflushed buffer would not.

**P+1 stage 2 now leaves a completed checkpoint, like P-1's already did.**
Found by re-running an already-finished exponent: P+1's stage 2 deleted its
checkpoint the moment it finished, so a later run at the identical bounds had
nothing to resume from and re-walked the whole `(B1, B2]` range from scratch
every time -- P+1's own stage 1 (and P-1's stage 1 and stage 2, since the 1.1
B2-extension work) already skipped straight to the gcd in this situation,
so this was a real inconsistency, not by design. `Pp1Stage2Save.h`'s
`complete` field already existed for this (kept for schema parity when B2
extension was scoped out of 1.2) but nothing ever set it. `runPP1Stage2` now
writes a completed record instead of deleting the checkpoint, and checks for
one before building a plan at all. This is a narrower fix than full B2
extension (seeding a *larger* B2's walk from a completed *smaller* one,
still out of scope -- see `Pp1Stage2Save.h`): only the exact-match case,
re-running the same `(exponent, b1, b2, seed)`, is affected. Raising B2 on a
finished P+1 stage 2 still re-walks the whole range.

**Audited a real run end to end** (an 8-digit exponent, M28733813, plus
earlier real runs on M786613 and M7873417 -- the latter found a genuine
3-factor composite, `236202511 * 2629453581823 * 12346147729361`) --
cross-checking every log line, `results.txt` entry, and checkpoint file
against what the code should have produced. The bounds model checked out
exactly: `--bounds` against the same config independently reproduces both
methods' chosen B1/B2 to the exact value the real run used. Two real
reporting defects turned up in the process, both fixed:

- **A multi-factor report was missing context depending on which method
  found it.** P-1's stage 1 printed "the gcd was a product of N factors"
  before a multi-factor announcement; the `reportFactors()` helper --
  despite being documented as "the shared helper for all factor
  announcements" and used by P+1's stage 1, P+1's stage 2, and P-1's stage
  2 -- never had that preamble, because P-1's stage 1 kept its own separate
  inline copy of the announcement loop instead of calling the helper. Seen
  live: P+1 announced M7873417's 3-factor composite with no preamble; P-1,
  re-finding the identical 3 factors 90 seconds later, announced them WITH
  it. Moved the preamble into `reportFactors()` and pointed P-1's stage 1 at
  it, so every call site is consistent and the duplicated loop is gone.
- **"appended to results.txt" didn't reliably mean a line was appended.**
  On P-1's side it printed unconditionally right after stage 1, including
  the common case (stage 2 about to run) where nothing has been written
  yet -- so a "no factor, stage 2 next" job showed the confirmation twice
  for one real write. On P+1's side, for `method = both`, the same message
  was gated behind the pp1-only early return and so never printed at all,
  even though P+1's own write did happen. Both fixed by moving the log line
  inside the branch that actually calls `writeResultJson` (or the
  `reportFactors` that wraps it).

Also surfaced, not a code bug but worth knowing: `results.txt`'s
`timestamp` is UTC (`gmtime_s`, matching what mersenne.org wants) while
`Mp_p-1_gpu.log`'s and the console's timestamps are local time
(`localtime`) -- the same event can legitimately show several hours apart
between the two, which is easy to mistake for a bug during exactly this
kind of audit. And `method = both` does not skip P-1 once P+1 already found
a factor -- P-1 still runs its own full stage 1 (and stage 2, if stage 1
doesn't find it) at full cost, confirmed live by the M7873417 run re-finding
the same factor a second time. Both are now documented in `README.md`'s
Scope section rather than being surprises; the `method = both` gap is also
tracked as outstanding work (`tasks/todo.md`) since the existing
within-method "stop, we're done" rule could in principle extend across
methods too.

**`config.txt`'s own comments had gone stale.** The `method` section still
described P+1 as borrowing P-1's B1 *and* B2 with no bounds model of its
own -- true when it was written, false since the B1 (and, later, B2) work
above landed, and never updated because a config template doesn't sit next
to the code that made it stale. Rewritten to describe the current,
independent-bounds behavior, and `pp1_seeds` gained an explicit comment that
it is a literal comma-separated list of seed VALUES, not a count --
`pp1_seeds = 3` requests exactly one seed (value 3), not three. That gap
was real enough to trip up a config written fresh during this same
project's own testing.

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
