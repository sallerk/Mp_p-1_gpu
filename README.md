# Mp_p-1_gpu 1.1

GPU **P-1** factoring of Mersenne numbers `M_p = 2^p - 1`, with both stages, and P+1 at currently stage 1 only.
Windows / MSVC, no external dependencies — no GMP, no CUDA or OpenCL SDK.

> **Licence: GPLv3.** Derived from [gpuowl / PRPLL](https://github.com/preda/gpuowl)
> by Mihai Preda and George Woltman — see [ATTRIBUTION.md](ATTRIBUTION.md).

---

## Download

Take `Mp_p-1_gpu-1.1-win64.zip` from the
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

## Results

`results.txt` gets one JSON object per line, the format
[mersenne.org](https://www.mersenne.org/manual_result/) accepts:

```json
{"status":"F","exponent":81679223,"worktype":"P-1","b1":2000000,"b2":60000000,
 "factors":["..."],"program":{"name":"Mp_p-1_gpu","version":"1.0"},
 "timestamp":"2026-07-28 17:29:54","user":"...","computer":"..."}
```

`status` is `F` (factor), `NF` (no factor), or `C` (a divisor that could not be
split into primes — recorded for you, not submittable). One line per job, with
both bounds when stage 2 ran.

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

## Self-tests

```
Mp_p-1_gpu.exe --selftest
```

`gcd`, `exponent`, `stage2plan` and `bounds` need no GPU. `engine`, `pm1`,
`pp1`, `extend`, `stage2` and `b2extend` exercise the GPU against exact CPU
arithmetic and against known factors of real Mersenne numbers.

## Scope and limitations

- **Windows / MSVC only.** No Makefile, no Linux build, no CI.
- **One exponent per run**, from `config.txt`. No worktodo queue, no PrimeNet
  automation — results are written for you to upload manually.
- **P+1 is stage 1 only**, and is a secondary mode. Prime95 has a more complete
  P+1 (both stages, CPU); this exists because no GIMPS *GPU* tool offers P+1 at
  all. Its yield per unit work is well below P-1's, so P-1 is the default.
- **P+1 currently borrows P-1's B1**, which is chosen for a different smoothness
  target.
- **Bound selection ignores work already done.** Raising B2 on a finished stage 2
  is cheap, but `auto` still picks bounds as though nothing had been computed
  yet. Set `b1` and `b2` explicitly if you are deliberately working an exponent
  upwards.
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
