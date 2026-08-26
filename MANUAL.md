# Mp_p-1_gpu — manual

GPU **P-1** and **P+1** factoring of Mersenne numbers `M_p = 2^p - 1`, both with
two stages. Windows / MSVC. Nothing to install to run it: no CUDA or OpenCL
SDK, and GMP (used for the CPU gcd step since 1.8) is linked statically, so
there's no DLL for it either.

> **Licence: GPLv3.** Derived from [gpuowl / PRPLL](https://github.com/preda/gpuowl)
> by Mihai Preda and George Woltman — see [ATTRIBUTION.md](ATTRIBUTION.md).

> This is the manual for the current source in this repository. What changed
> between releases is in [CHANGELOG.md](CHANGELOG.md); every earlier version is
> on the [Releases page](https://github.com/sallerk/Mp_p-1_gpu/releases), each
> at its own tag.

> **1.7 builds stage 2's `T_j` table incrementally.** The table is the
> `x^(j^2)` values every accumulator multiply is taken against. It used to be
> built with one exponentiation per entry, about `2*log2(j)` multiplies each;
> it is now a second-difference chain costing one multiply per unit of `j`
> regardless of how many entries are kept. Nothing about the result changes --
> the accumulator comes out bit for bit identical -- it is purely the setup
> that got cheaper. The win scales with how large the table is, so it is
> biggest exactly where the table is biggest: 27% off stage 2 at a 7-digit
> exponent with 2160 buffers, a few percent at a GIMPS-wavefront exponent
> where GPU memory only allows a few hundred.
>
> That reaches the end of a run: the final gcd needs the stage-2 accumulator,
> so stage 1, stage 2 and the final gcd run in series, and a second saved in
> stage 2 is a second saved overall. The stage-1 gcd runs alongside stage 2 --
> that is what the overlap covers, not stage 2 itself. As of 1.7 it no longer
> covers all of it: stage 2 is now short enough that the stage-1 gcd outlasts
> it and overlaps the final gcd, and the two compete for the same worker
> threads. See CHANGELOG.md.
>
> Carried over from 1.6, which fixed a real correctness bug in
> checkpoint/resume: if you have ever resumed a stage-1 run under 1.5 or
> earlier, its result should be treated as unverified.

---

## Download

Take the zip from the
[latest release](https://github.com/sallerk/Mp_p-1_gpu/releases/latest), unzip
it anywhere, and run it. There is nothing to install:

- **no Visual Studio and no compiler** — the C runtime is linked statically
- **no OpenCL SDK** — OpenCL is loaded from your GPU driver at run time
- **no GMP install** — GMP (used for the gcd step, since 1.8) is linked
  statically into the exe, not a separate DLL or runtime dependency
- no CUDA, no Python

The only requirement is a GPU driver with OpenCL support, which the AMD, nVidia
and Intel desktop drivers all install by default. 64-bit Windows; the binary
targets the baseline instruction set, so any x64 CPU will run it.

The executable is not code-signed, so the first launch will probably raise
SmartScreen — *"Windows protected your PC"* → **More info** → **Run anyway**.
If you would rather not run a stranger's binary, build it yourself; it takes
about a minute.

## Build from source (optional)

1. Install Visual Studio 2019 or 2022 with the **Desktop development with C++**
   workload. This includes a copy of [vcpkg](https://github.com/microsoft/vcpkg)
   already — nothing extra to install for step 2.
2. Install GMP, once:
   ```
   "%VCPKG_ROOT%\vcpkg.exe" install gmp:x64-windows-static
   ```
   (`VCPKG_ROOT` is set automatically by Visual Studio's build environment,
   which `build.bat` also uses — run this from a "Developer Command Prompt for
   VS", or just try `build.bat` first: it prints this exact command, with the
   real path filled in, if GMP isn't found. If you have an older Visual
   Studio without a bundled vcpkg, clone
   [vcpkg](https://github.com/microsoft/vcpkg) yourself, run
   `bootstrap-vcpkg.bat`, then the `install` line above, and either set
   `VCPKG_ROOT` to that folder or place it at `Mp_p-1_gpu\vcpkg`.)
3. Run `build.bat`. It produces `Mp_p-1_gpu.exe`.

Python is optional — it regenerates `src/bundle.cpp` from the OpenCL sources,
and the generated file is committed so the build works without it.

## Quick start

Add an exponent to `worktodo.txt` (one per line — see that file's own header
comment for the queue format):

```
81679223
```

Edit `config.txt` for everything about *how* to work it:

```
factored_to = 76         # bits already trial-factored, from mersenne.org
b1 = auto
b2 = auto
```

Then, once for this exponent:

```
Mp_p-1_gpu.exe --tune noconfig
```

(Optional — a run without it picks and verifies its own transform. See
[Tuning](#tuning) for why it is per exponent and not per GPU.)

and to run the job:

```
Mp_p-1_gpu.exe
```

## config.txt

Every setting is documented inline in the file. The ones you will actually
touch:

| key | meaning |
|---|---|
| `worktodo_file` | the exponent queue, one per line — `worktodo.txt` by default |
| `factored_to` | trial-factoring depth already done, in bits |
| `b1`, `b2` | `auto`, or explicit bounds |
| `stages` | `auto`, `both`, or `1` for stage 1 only |
| `method` | `pm1` (default), `pp1`, or `both` |
| `pp1_seeds` | comma-separated bases to try for P+1, default `3, 5, 7` |
| `bias` | what a factor is worth to you, in units of one PRP test |
| `bounds_tolerance` | how much extra work to accept for a better chance |
| `gcd_threads` | CPU threads for the gcd phase — no longer consulted by the default gcd (GMP's `mpz_gcd`, single-threaded, since 1.8); still respected by this program's own `gcdHalf` implementation, kept as the reversion path |
| `username`, `computer_name` | written into results for PrimeNet |
| `wait_for_work` | `no` (default) exits on an empty queue; `yes` waits and rechecks it — for running alongside AutoPrimeNet |
| `wait_poll_seconds` | how often to recheck, when `wait_for_work = yes` — default 5 |

## worktodo.txt — the exponent queue

As of 1.4, exponents no longer live in `config.txt` — they live in
`worktodo.txt` (or whatever `worktodo_file` points at), one per line,
processed top to bottom. `#` starts a comment, blank lines are ignored, same
as `config.txt`. A bare exponent is this program's own shorthand: every such
entry runs under `config.txt`'s settings — `method`, `b1`/`b2`, `bias`,
`stages`, and everything else — only the exponent (and, for `factored_to =
auto`, its default) varies per entry.

A line can also be a `Pfactor=` or `Pminus1=` assignment — the shape
[AutoPrimeNet](https://github.com/tdulcet/AutoPrimeNet) (the maintained
successor to gpuowl's `primenet.py`) and Prime95 itself write. This program
never talks to PrimeNet itself — no networking of any kind was added — but
AutoPrimeNet runs as a separate program in the same directory, fetching
assignments as these lines and uploading new `results.txt` lines, so the two
tools share a queue without this program needing to know PrimeNet exists.
Such a line carries its own `k,b,n,c` (validated as a Mersenne number — `k=1,
b=2, c=-1` — anything else is a hard error, since this program factors
nothing else) and, for `Pminus1=`, its own `B1`/`B2`:

```
Pminus1=1,2,86243,-1,50000,3000000
```

That B1/B2 is used as-is whenever present. A `Pfactor=` line or a bare
exponent has none to honor, and falls back to `config.txt`'s own `b1`/`b2`
(auto by default) exactly as before. Either way, the entry's exponent,
trial-factoring depth (`Pfactor=`'s `how_far_factored` or `Pminus1=`'s
`sieve_depth`, when present — it wins over `factored_to = auto`'s built-in
default, though an *explicit* `factored_to = N` in config.txt still wins
over that), assignment ID, and known factors are always taken from the line.
An assignment ID and known factors, when present, are echoed back into
`results.txt` (see **Results**, below) so AutoPrimeNet's own upload step
sees what it expects.

One thing an assignment can carry that this program cannot use yet:
`Pminus1=`'s optional `B2_start`, meaning stage 2 was already partly walked
elsewhere. There is no local checkpoint to seed that partial progress from
(the `extend` reuse below needs one it can verify matches), so such a job
prints a warning and walks the full `(B1,B2]` range instead — always
correct, just not the minimal amount of work.

A finished job's line is removed once its result is written to
`results.txt`, and the next line starts — in the same run, no relaunch
needed. What "finished" means, and what happens otherwise:

- **found a factor, or ran to completion with none** — line removed, next
  entry starts.
- **interrupted** (Ctrl-C, or the process killed) — line left in place.
  Relaunching resumes that exponent from its checkpoint exactly as before,
  then continues down the queue from there.
- **failed** (an error — not "no factor found"; e.g. an exponent too small
  for any FFT shape) — the whole queue stops there and the line is left in
  place, so you can fix or remove it before continuing. A bad entry does not
  get silently skipped or silently retried forever.

`--bounds` and `--tune` read the queue's first entry (without consuming it)
to know which exponent to scope themselves to, the same role `config.txt`'s
`exponent` key used to play.

## Command line

| | |
|---|---|
| *(no arguments)* | run the job in `config.txt` |
| `--config <file>` | use a different job file |
| `--bounds` | show the B1/B2 trade-off for your exponent |
| `--tune[=opts]` | measure FFT shapes and kernel options for your GPU, targeting the worktodo exponent (see [Tuning](#tuning)) |
| `--bench` | time every FFT that can hold an exponent |
| `--list-devices` | list OpenCL GPUs |
| `--selftest[=which]` | run self-checks |
| `-d <n>` | use GPU *n* |
| `-fft <spec>` | force an FFT configuration |
| `-h` | full help |

## Tuning

`--tune` writes two files, and they are not the same kind of thing.

**`tune.txt`** — the cost of each FFT *shape* it managed to time. A run reads
it as a ranking of candidates, verifies them at your exponent, and can fall
through to untuned shapes if every tuned entry fails or is the wrong size. It
is advisory, and a stale one costs you a slower transform at worst.

**`Mp_p-1_gpu-tune-config.txt`** — the winners of the `-use` *kernel option*
search: `MODM31`, `TABMUL_CHAIN61`, `IN_WG`, and about a dozen more. These are
not advisory. They change how the kernels are compiled, and they are specific
to the transform they were measured on: one such file replayed against a
transform it was not measured for cost about 32% per squaring, worse than no
tuning at all.

So each block in it is tagged with the exponent range the tune targeted:

```
# New settings based on a -tune run.
# tuned-for exponents 82589933-82589933
   -use IN_WG=64,IN_SIZEX=8,OUT_WG=256,...
```

A run applies a block only when its range covers the exponent being worked, and
otherwise uses the stock defaults and says so. A bare `--tune` targets the first
exponent in `worktodo.txt`, so that range is normally one exponent — which is
why the recommended workflow tunes per exponent, not once per GPU. Pass
`minexp=`/`maxexp=` to tune for a range on purpose.

The file is *appended* to, so tuning several exponents accumulates blocks in it
and each is used for its own. Blocks written by a version before this existed
carry no range; they are ignored and reported, because there is no way to know
what they were measured on. Re-run `--tune` to replace them.

A queue is handled per entry: the options are re-decided for each exponent in
`worktodo.txt` as it comes up, not once when the program starts.

`--tune` itself, and `--selftest`, always measure from the stock defaults — a
tune has to be reproducible, and a selftest should test the engine rather than
one tuning. `--bench` does apply a tuning matching its exponent, since the
number it reports is meant to be the number a real run of that exponent gets.

`--tune` picks the transform before it tunes anything: it runs the same
selector a real run uses, at the exponent it is targeting, and measures the
options against whatever that returns. So the recorded range describes where
the numbers actually came from. Options for arithmetic the chosen transform
does not contain are not searched at all, which is why tuning a small exponent
takes fewer steps than a wavefront one.

It also sizes its own measurements. `quick=` is an iteration count, and the
right number of iterations depends entirely on the transform — 400 of them is
three seconds on a 7.5M-word FFT and twenty milliseconds on a 262144-word one.
`--tune` probes the chosen shape and uses the shortest run whose timed window
still clears 1.5 seconds, so a `quick=` you pass is treated as a ceiling on
speed: it can ask for more accuracy, never less.

And it measures how much it can trust itself. Before searching, it times the
chosen configuration several times over and takes the spread as a noise floor;
an option has to beat the current setting by more than that to be written as an
active `-use` line. Anything closer goes to the commented "slightly faster"
block instead, which is what that block is for. Expect a loaded or hot GPU to
admit fewer settings than an idle one — that is the threshold doing its job,
not the tune failing.

A tune times **every** candidate transform, not the two or three a job start
settles for — it has been asked for the measurement, so it takes the time.
Expect the FFT selection at the head of a tune to run for minutes rather than
the ~20 seconds a job start spends.

`Ctrl-C` stops a tune. Whatever was already written to `tune.txt` and
`Mp_p-1_gpu-tune-config.txt` is kept; the option search appends its results
before the shape sweep begins, so a tune stopped late still leaves something
usable.

**What to expect from it.** On the machine this was developed on, the kernel
options do not produce a win these measurements can resolve: tuned and stock
come out level at both a small exponent and a wavefront one. What `--tune`
reliably contributes is `tune.txt` — the shape ranking — and that only saves
the FFT search some work, since the search finds the same transform without it.
Treat a tune as optional, and treat an empty active block as an honest answer
rather than a failure.

### Picking the transform

With `fft = auto` (the default), a run ranks candidates from `tune.txt`, drops
any that cannot hold the exponent or that are more than twice the smallest
usable size, and verifies them in order against a correctness check.

It does not stop at the first one that passes: it times two or three candidates
and takes the fastest. Shapes that all pass the correctness check have measured
4.9x apart at one exponent. It never stops before two have actually been
compared, so neither a slow first candidate nor a run of rejected ones can
leave a job on an unexamined transform.

The wait scales with the transform, since that is what costs: about 7 seconds
at the smallest exponent the engine accepts, 15–25 at a wavefront one, and more
when candidates fail their correctness check and have to be skipped.

Untuned candidates are ordered by **precision headroom**: smallest transform
first, and within one size, the cheapest arithmetic whose bits-per-word ceiling
still clears the exponent. That is an ordering, not a decision — FP64 is the
one family it ranks wrongly on this hardware (lowest ceiling, but 4.3x slower
than the NTTs here), and it loses on the measurement instead. Every family here buys precision with speed, so the
tightest fit is the one to try first, and one variant of every shape is tried
before any shape gets a second (the `:101` and `:202` of a shape measure within
0.2% of each other). Entries from `tune.txt` keep their measured order and are
tried ahead of all of it.

A caveat on precision: GPU clocks vary about 10% run to run, so the search
cannot reliably separate two shapes closer than roughly 5%. It is built to
avoid the large mistake, not to split hairs — before this ordering existed it
could settle on a transform 81% slower than one it had already timed.

Verdicts *and timings* are remembered in `fft-verified.txt`, keyed by exponent,
GPU and driver, so a repeat run of the same exponent re-picks the same winner
instantly. Delete that file to search again. Each line records the measurement
scale it was taken at; costs from a different scale are not comparable, so an
entry from an older build keeps its verdict but is re-timed. `-fft <spec>` skips selection
entirely (the spec is still verified unless `verify_fft = 0`).

## Reading the output

A run reports its phases as `[2/5 stage 1 GPU]` — five phases when a stage 2
follows, three without:

```
[1/5 exponent CPU]   building the stage-1 exponent E
[2/5 stage 1 GPU]    the squaring ladder
[3/5 gcd CPU]        gcd(x-1, M_p)      GPU idle here; this is normal
[4/5 stage 2 GPU]    the pairing walk
[5/5 gcd CPU]        gcd(acc, M_p)
```

Progress lines rewrite in place on a terminal, and are written once a minute
when redirected to a file.

The gcd phases report a percentage that jumps straight from 0% to roughly
**62%**, then fills in the rest more smoothly. That is not a stall: the
half-GCD halves its operand's bit length each pass, and each halving costs
roughly half of what the one before it did, so the percentage is weighted by
work (`1 - (bits left / bits start)^1.39`, the exponent measured from the
algorithm's own cost curve), not by bits. The very first pass -- on the
full-size operand, the single most expensive step in the whole gcd -- reports
no progress at all while it runs, then lands on `1 - 0.5^1.39 ≈ 61.8%` the
instant it finishes. This is expected on every gcd, stage 1 or stage 2, P-1 or
P+1 -- they all share the same `gcdWithProgress` reporting.

## Results

`results.txt` gets one JSON object per line, the format
[mersenne.org](https://www.mersenne.org/manual_result/) accepts:

```json
{"status":"F","exponent":81679223,"worktype":"P-1","b1":2000000,"b2":60000000,
 "factors":["..."],"program":{"name":"Mp_p-1_gpu","version":"1.8"},
 "timestamp":"2026-07-28 17:29:54","user":"...","computer":"...",
 "aid":"...","known-factors":["..."]}
```

`status` is `F` (factor), `NF` (no factor), or `C` (a divisor that could not be
split into primes — recorded for you, not submittable). One line per job, with
both bounds when stage 2 ran. `aid` and `known-factors` appear only when the
job came from a `Pfactor=`/`Pminus1=` worktodo entry that carried them (see
**worktodo.txt**, above) — present or absent the same way `user`/`computer`
already are, and using the same JSON keys Prime95's own submission does, so
AutoPrimeNet's upload step can match the line back to its assignment.

Its `timestamp` is **UTC**, deliberately — it is the field mersenne.org
actually wants. The console and `Mp_p-1_gpu.log` timestamps are your machine's
**local** time instead, since that is what a person watching the run wants. Do
not be surprised if a line in the log and the `results.txt` entry it caused
carry timestamps several hours apart; both are correct, just in different
zones.

## Choosing bounds

`b1 = auto` weighs the chance of a factor against the work. Run `--bounds` to
see the trade-off and what `bounds_tolerance` buys:

```
tolerance   bounds                P(factor)   work (PRP=1)
0.00        B1=0.20M B2=6M           2.363%      0.022
0.05        B1=2.00M B2=60M          5.059%      0.100   <-- default
0.20        B1=4.00M B2=240M         6.720%      0.273
```

Raise `bounds_tolerance` to spend more for a better chance; raise `bias` if a
factor is worth more to you than one PRP test.

P-1 and P+1 have their own, independent B1 models (P+1 targets `q+1`, a
different — and harder to hit — smoothness condition than P-1's `q-1`; see
the P+1 section below), so a `method = both` job with `b1 = auto` runs each
method at its own B1, shown separately at startup and in `--bounds`. P+1's
auto-chosen B1 is typically much smaller than P-1's at the same exponent.

## Resuming, and raising the bounds

Runs checkpoint automatically and resume when restarted. Ctrl-C is safe.

Completed work is reused rather than repeated, so raising a bound costs only the
difference:

- **raise `b1`** — a completed smaller B1 is extended to the new one; the saving
  is exactly `oldB1 / newB1` of the squarings.
- **raise `b2`** — a completed stage 2 leaves its accumulator behind, and the
  next run walks only `(oldB2, newB2]`. The pairing shape may differ between the
  two; nothing about the accumulator requires it to match.

So the sensible way to work an exponent is to start at modest bounds and raise
them if nothing turns up. Re-running an exponent that is already finished at
those bounds does no GPU work at all.

Both are on by default; `extend = 0` in `config.txt` turns them off.

Checkpoints are named `pm1_<exponent>_b1_<B1>.save` (stage 1) and
`pm1_<exponent>_s2_<B1>_<B2>.save` (stage 2). Deleting them costs work, never
correctness — every field that could make a residue mean something else is
recorded and checked, so a file that does not match the job is refused, not
adapted.

## P+1

`method = pp1` (or `both`, to try P+1 then fall through to P-1) runs P+1
instead of, or alongside, P-1. P+1 finds a factor q of `M_p` when `q+1` is
B1-smooth — a different target from P-1's `q-1`, so it catches factors P-1
cannot, and vice versa, for the seeds tried (`pp1_seeds`, default `3, 5, 7`).

P+1 has its own B1 and B2 model (`pp1Prob`/`choosePP1Bounds` in
`src/Bounds.cpp`) — both used to silently borrow P-1's own choice, which
optimises for the wrong smoothness target and, for B2, could tie P+1's B1 to
an inflated range it never asked for; now `b1 = auto` and `b2 = auto` pick
each method's bounds independently. Stage 2 catches `q+1` that is B1-smooth
apart from one prime in `(B1, B2]`, exactly like P-1's stage 2, and reuses
P-1's own measured stage-2 cost constants (same pairing geometry, same
underlying multiply — see the comment above `choosePP1Bounds` in
`src/Bounds.h`), but still shares P-1's `stage2_d`/`stage2_w` pairing *shape*
— a GPU-memory budget decision (one T-table sized for whichever method needs
the smaller B1), not a cost-model gap. `stages` still switches stage 2 off
for both methods at once. Checkpointing and
interrupted-walk resume work the same way as P-1's, per seed:
`pp1_<exponent>_b1_<B1>_s<seed>.save` (stage 1) and
`pp1_<exponent>_s2_<B1>_<B2>_s<seed>.save` (stage 2). Raising B2 on a
completed P+1 stage 2 does **not** yet reuse the earlier walk the way P-1's
does — see Scope below.

The pairing math is different from P-1's under the hood (P-1's exact trick
needs an explicit ring element, which P+1 never has — only a Lucas-sequence
trace), derived and verified independently rather than ported; see the
comment above `Gpu::pp1Stage2` in `src/Gpu.h` for the full derivation.

## How stage 2 pairs primes

Stage 2 catches a factor `q` whose `q-1` is B1-smooth apart from one prime in
`(B1, B2]`. Done naively that is one accumulator multiply per prime. Montgomery
pairing gets two: pick a highly composite `D`, write each prime as
`m*D ± j`, and note that with

```
A_m = x^((m*D)^2)      T_j = x^(j^2)
A_m - T_j  is divisible by both  m*D-j  and  m*D+j
```

a single multiply by `A_m - T_j` covers *both* candidates whenever both happen
to be prime. `stage2_w` widens `j` from `D/2` to `w*D/2`, which gives each
prime about `w` candidate partners instead of exactly one, turning slot
assignment into a matching problem worth solving. That costs `w` times as many
`T_j` buffers, and GPU memory is what caps it.

The pairing rate that buys, measured by `--selftest=stage2plan`:

| candidates per prime | `T_j` buffers | primes per multiply |
|---|---|--:|
| `w=1` | `phi(D)/2` | 1.15–1.21 |
| `w=3` | `3*phi(D)/2` | 1.35–1.43 |
| `w=5` | `5*phi(D)/2` | 1.47–1.54 |
| `w=7` | `7*phi(D)/2` | 1.53 |
| `w=9` | `9*phi(D)/2` | 1.56 |

Two things worth reading off it. The rate depends on `w` and barely on `D` —
it is candidates per prime that matters, not block size. And it flattens fast:
`w=1→3` buys 0.2 primes per multiply, `w=7→9` buys 0.03. The ceiling is 2.0
and it is not reachable, because plenty of primes have no prime partner to
pair with at any window.

**On PrMers' "V-trace"/"Pair95".** [PrMers](https://github.com/cherubrock-seb/PrMers)
runs a different-looking stage 2, built on the Lucas trace `V_n = H^n + H^-n`,
where `V_(kD) - V_j` covers `kD-j` and `kD+j` in one scalar. Pair95 extends it
with irregular offsets `unit`, `unit+D`, `unit+3D`, `unit+7D`, … This was
scoped for 1.7 and deliberately not ported, for reasons worth recording:

- **The pairing is the same trade, priced the same.** Pair95's `L` levels give
  each prime `L` candidate partners for `L*phi(D)/2` table entries. The `w`
  window above gives `w` candidates for `w*phi(D)/2` entries — the same
  candidates per unit of memory, on the same curve, and symmetric about the
  block rather than one-sided.
- **A real port would cost more than it saves.** V-trace for P-1 needs
  `V_1 = x + x^-1 mod M_p`, and that inverse is an extended gcd on a p-bit
  number. The plain gcd already runs for minutes at a wavefront exponent; the
  extended one would exceed the whole of the stage 2 it was meant to speed up.
  P+1's stage 2 here is already a Lucas-trace walk, because P+1's group
  element is *only* ever available as a trace — there the inverse is free, and
  the same construction is the natural one.

What did carry over is the part that has nothing to do with pairing: a linear
index lets the baby table be built by a recurrence instead of an
exponentiation per entry. `x^(j^2)` is quadratic in `j`, so it needs a second
difference rather than a first, but it works the same way and needs no
inverse — which is what 1.7 does, and where its speedup comes from.

## Stage-2 memory, and the largest exponent that fits

`stage2()` allocates `J + 9` full-residue buffers: the `T_j` table
(`J = phi(D)*w/2`), six scratch (`x`, `acc`, `A`, `S`, `C`, `diff`), and three
more, transient, for the table-building chain. One residue is
`8 * fft_words` bytes. **Nothing in it is sized by B1 or B2** -- a longer walk
changes the multiply count, not the allocation. Measured at a fixed exponent,
FFT and shape, B2 = 2M against B2 = 80M (40x the primes): 701 MB against
711 MB, which is sampling noise.

B1 does bound the *shape*, and only downward: `chooseStage2Shape` rejects any
`(D, w)` with `w*D/2 > B1`, since every prime in `(B1, B2]` has to exceed the
window half-width. At B1 = 200 the only shape left is the floor; by B1 = 1000
the constraint stops binding and the memory budget takes over. Below B1 = 105
even the floor violates it and stage 2 throws -- not reachable with any
auto-chosen bound, but worth knowing if you pin a tiny `b1` by hand.

That floor is `D=210 w=1`, 24 buffers: `pickStage2Shape` falls back to it when
nothing fits the budget, so it is also the point at which stage 2 stops being
able to run at all. Its cost, and the largest exponent each transform size can
hold:

| FFT words | one residue | max exponent | stage-2 floor |
|--:|--:|--:|--:|
| 262,144 | 2 MB | 13,046,875 | 203 MB |
| 524,288 | 4 MB | 26,015,625 | 276 MB |
| 1,048,576 | 8 MB | 51,953,125 | 422 MB |
| 2,097,152 | 16 MB | 99,785,156 | **701 MB** |
| 4,194,304 | 32 MB | 198,515,625 | **1258 MB** |
| 8,388,608 | 64 MB | 394,240,722 | **2400 MB** |
| 16,777,216 | 128 MB | 746,905,516 | **4940 MB** |
| 33,554,432 | 256 MB | 1,399,804,687 | 9.5 GB |
| 67,108,864 | 512 MB | 2,650,781,250 | 18.8 GB |

Bold entries are measured (RTX 3070, 1.8.1, two shapes per size so the table
cost and the fixed overhead separate cleanly); the rest follow the fit those
four give, `130 MB + (8*buffers + 100) * fft_words`, which reproduces them to
within 3%. Max exponents come from the transform tables via `--bounds`.

Read off as a card-size limit:

| GPU RAM | headless | with a desktop on it (~1.4 GB) |
|---|--:|--:|
| 2 GB | 198 M | 52 M |
| 4 GB | 394 M | 394 M |
| 6 GB | 747 M | 394 M |
| 8 GB | 747 M | 747 M |
| 10 GB | 1.40 G | 747 M |
| 11-16 GB | 1.40 G | 1.40 G |
| 24 GB and up | 2.65 G | 2.65 G |

Four things this table is not:

- **It is not a speed target.** The floor is where stage 2 still *runs*. At a
  GIMPS-wavefront exponent the budget picks 144 buffers or more, and running at
  24 costs roughly a quarter of stage-2 throughput. Pin `stage2_d` /
  `stage2_w` if you want the floor deliberately.
- **There is no graceful degradation past it.** When even 24 buffers exceed the
  budget, `pickStage2Shape` still returns `D=210 w=1` -- there is no "stage 2
  declined, not enough memory" path, so the allocation simply fails.
- **Verification can bite before memory does.** At 16,777,216 words on this
  card only `1:1K:8:1K:202` passed its correctness check; the seven other shapes
  tried at that size, across all four arithmetic types, failed. Near the top of
  a size class the usable limit can be below the memory limit.
- **It assumes the smallest fitting transform.** A real run picks by measured
  cost among verified configs and can land a size class up, doubling every
  number here. For M100000223 it picked the smallest (4,194,304 words).

## Self-tests

```
Mp_p-1_gpu.exe --selftest
```

`gcd`, `exponent`, `stage2plan`, `bounds` and `worktodo` need no GPU. `engine`,
`pm1`, `pp1`, `extend`, `stage2`, `b2extend` and `pp1stage2` exercise the GPU
against exact CPU arithmetic and against known factors of real Mersenne
numbers.

## Scope and limitations

- **Windows / MSVC only.** No Makefile, no Linux build, no CI.
- **No PrimeNet automation by this program itself** — no networking of any
  kind, so no auto-fetch and no auto-submit. As of 1.8 `worktodo.txt`/
  `results.txt` understand assignment IDs and known factors, the shape
  [AutoPrimeNet](https://github.com/tdulcet/AutoPrimeNet) writes and reads —
  see **worktodo.txt** and **Results**, above — so running AutoPrimeNet
  alongside this program as a separate process covers the fetch/submit half.
  Without it, results are written for you to upload manually, same as always.
- **`Pminus1=`'s `B2_start` is parsed but not honored** (see **worktodo.txt**,
  above) — a job carrying one prints a warning and walks the full `(B1,B2]`
  range rather than only the uncovered part.
- **The gcd phase can no longer be interrupted once started, only before.**
  Since 1.8 it runs through GMP's `mpz_gcd`, which has no progress callback
  (the previous implementation's Ctrl-C support relied on one). Ctrl-C is
  still honored right up to the moment the call starts; the window this
  narrows is small in practice — GMP's own worst case at real production
  scale measured ~12s, against the previous implementation's ~140s.
- **Minimum exponent is 786613.** Every FFT shape needs at least 3 bits/word
  (`FFTShape::minBpw()`); the smallest shape in the catalog is 256x256x2
  (262144 words), so nothing below this exponent has any usable shape at
  all — `--bounds`, `--tune`, and a real run all just fail with "no FFT fits
  this exponent".
- **Maximum exponent is 4294967295.** Exponents are held in a `u32`
  (`WorktodoEntry::exponent`, `Config::exponent`), so anything larger cannot
  be represented at all. The worktodo parser refuses such a line with the
  exponent quoted back; before 1.9 it wrapped silently and ran a different
  exponent, and 1.9 fixed only the `Pfactor=`/`Pminus1=` path -- a bare
  exponent still wrapped until 1.9.1. Both shapes are covered by the
  refusal tests now. In practice the transform catalog runs out well
  before the type does.
- **The exponent must be prime.** `M_p = 2^p - 1` is only a GIMPS candidate
  for prime `p`, and for composite `p` every `d | p` gives an algebraic
  factor `2^d - 1` -- so a run against a mistyped composite exponent can
  "find" a factor and report it as though it were a discovery. Enforced at
  the parser since 1.9.1; before that nothing checked it, despite
  config.txt having always said so.
- **P+1 is a secondary mode.** Its yield per unit work is well below P-1's, so
  P-1 is the default; run it in earnest only once P-1 has been tried.
- **`method = both` does not stop early across methods.** Finding a factor in
  P+1 skips that method's own remaining seeds and its stage 2 (a factor
  already in hand makes more searching wasted work) — but P-1 still runs
  afterward regardless, at full cost. The equivalent "stop, we're done"
  check exists *within* each method, not *between* them yet.
- **P+1 has its own B1 and B2 model now, but still shares P-1's pairing
  shape** (`stage2_d`, `stage2_w`) — a GPU-memory budget decision (one
  shared T-table) rather than a cost-model gap.
- **No B2 extension for P+1.** P-1's stage 2 reuses a completed smaller-B2
  walk when you raise B2 (see Resuming above); P+1's does not yet — raising B2
  re-walks the whole range from B1 rather than only the new part. An
  identical re-run at the *same* B2 is not affected by this — like P-1's, a
  completed P+1 stage 2 leaves its accumulator behind and a later run at the
  same bounds skips straight to the gcd. The checkpoint format has room for
  true extension later.
- **Bound selection ignores work already done.** Raising B2 on a finished stage 2
  is cheap for P-1 (not yet for P+1, per above), but `auto` still picks bounds
  as though nothing had been computed yet either way. Set `b1` and `b2`
  explicitly if you are deliberately working an exponent upwards.
- Not a drop-in replacement for gpuowl: this does factoring only, no PRP or LL.

## Licence

GPLv3 — see [LICENSE](LICENSE) and [ATTRIBUTION.md](ATTRIBUTION.md).
The upstream copyright notices in `src/` are required by GPLv3 §5 and must not
be removed.

The corresponding source for a released binary is the tagged commit it was built
from, in this repository — which satisfies GPLv3 §6. The zip carries `LICENSE`
and `ATTRIBUTION.md` with it.

## Name

"M_p, P-1, GPU". It implements P-1 and P+1 — it is **not** an ECM
implementation.
