# Mp_p-1_gpu 1.3

GPU **P-1** and **P+1** factoring of Mersenne numbers `M_p = 2^p - 1`, both with
two stages. Windows / MSVC, no external dependencies — no GMP, no CUDA or
OpenCL SDK.

> **Licence: GPLv3.** Derived from [gpuowl / PRPLL](https://github.com/preda/gpuowl)
> by Mihai Preda and George Woltman — see [../ATTRIBUTION.md](../ATTRIBUTION.md).

> This is the 1.3 release, kept exactly as it shipped. What changed since is in
> [../CHANGELOG.md](../CHANGELOG.md); the current version is
> [../1.7](../1.7/README.md).

---

## Download

Take `Mp_p-1_gpu-1.3-win64.zip` from the
[Releases page](https://github.com/sallerk/Mp_p-1_gpu/releases), unzip it
anywhere, and run it. There is nothing to install:

- **no Visual Studio and no compiler** — the C runtime is linked statically
- **no OpenCL SDK** — OpenCL is loaded from your GPU driver at run time
- no CUDA, no GMP, no Python

The only requirement is a GPU driver with OpenCL support, which the AMD, nVidia
and Intel desktop drivers all install by default. 64-bit Windows; the binary
targets the baseline instruction set, so any x64 CPU will run it.

The executable is not code-signed, so the first launch will probably raise
SmartScreen — *"Windows protected your PC"* → **More info** → **Run anyway**.
If you would rather not run a stranger's binary, build it yourself; it takes
about a minute.

## Build from source (optional)

1. Install Visual Studio 2019 or 2022 with the **Desktop development with C++**
   workload.
2. Run `build.bat`. It produces `Mp_p-1_gpu.exe`.

Nothing else is required. Python is optional — it regenerates `src/bundle.cpp`
from the OpenCL sources, and the generated file is committed so the build works
without it.

## Quick start

Edit `config.txt`:

```
exponent = 81679223      # the p in M_p = 2^p - 1; must be prime
factored_to = 76         # bits already trial-factored, from mersenne.org
b1 = auto
b2 = auto
```

Then, once per GPU:

```
Mp_p-1_gpu.exe --tune noconfig
```

and to run the job:

```
Mp_p-1_gpu.exe
```

## config.txt

Every setting is documented inline in the file. The ones you will actually
touch:

| key | meaning |
|---|---|
| `exponent` | the Mersenne exponent to factor |
| `factored_to` | trial-factoring depth already done, in bits |
| `b1`, `b2` | `auto`, or explicit bounds |
| `stages` | `auto`, `both`, or `1` for stage 1 only |
| `method` | `pm1` (default), `pp1`, or `both` |
| `pp1_seeds` | comma-separated bases to try for P+1, default `3, 5, 7` |
| `bias` | what a factor is worth to you, in units of one PRP test |
| `bounds_tolerance` | how much extra work to accept for a better chance |
| `gcd_threads` | CPU threads for the gcd phase |
| `username`, `computer_name` | written into results for PrimeNet |

## Command line

| | |
|---|---|
| *(no arguments)* | run the job in `config.txt` |
| `--config <file>` | use a different job file |
| `--bounds` | show the B1/B2 trade-off for your exponent |
| `--tune[=opts]` | measure FFT configurations for your GPU |
| `--bench` | time every FFT that can hold an exponent |
| `--list-devices` | list OpenCL GPUs |
| `--selftest[=which]` | run self-checks |
| `-d <n>` | use GPU *n* |
| `-fft <spec>` | force an FFT configuration |
| `-h` | full help |

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
 "factors":["..."],"program":{"name":"Mp_p-1_gpu","version":"1.3"},
 "timestamp":"2026-07-28 17:29:54","user":"...","computer":"..."}
```

`status` is `F` (factor), `NF` (no factor), or `C` (a divisor that could not be
split into primes — recorded for you, not submittable). One line per job, with
both bounds when stage 2 ran.

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

## Self-tests

```
Mp_p-1_gpu.exe --selftest
```

`gcd`, `exponent`, `stage2plan` and `bounds` need no GPU. `engine`, `pm1`,
`pp1`, `extend`, `stage2`, `b2extend` and `pp1stage2` exercise the GPU against
exact CPU arithmetic and against known factors of real Mersenne numbers.

## Scope and limitations

- **Windows / MSVC only.** No Makefile, no Linux build, no CI.
- **One exponent per run**, from `config.txt`. No worktodo queue, no PrimeNet
  automation — results are written for you to upload manually.
- **Minimum exponent is 786613.** Every FFT shape needs at least 3 bits/word
  (`FFTShape::minBpw()`); the smallest shape in the catalog is 256x256x2
  (262144 words), so nothing below this exponent has any usable shape at
  all — `--bounds`, `--tune`, and a real run all just fail with "no FFT fits
  this exponent".
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

GPLv3 — see [../LICENSE](../LICENSE) and [../ATTRIBUTION.md](../ATTRIBUTION.md).
The upstream copyright notices in `src/` are required by GPLv3 §5 and must not
be removed.

The corresponding source for a released binary is the tagged commit it was built
from, in this repository — which satisfies GPLv3 §6. The zip carries `LICENSE`
and `ATTRIBUTION.md` with it.

## Name

"M_p, P-1, GPU". It implements P-1 and P+1 — it is **not** an ECM
implementation.
