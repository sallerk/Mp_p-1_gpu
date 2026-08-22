# Mp_p-1_gpu

GPU **P-1** and **P+1** factoring of Mersenne numbers `M_p = 2^p - 1`, both with
two stages. Windows / MSVC. Nothing to install to run it: no CUDA or OpenCL
SDK, and GMP (used for the CPU gcd step since 1.8) is linked statically, so
there's no DLL for it either.

> **Licence: GPLv3.** Derived from [gpuowl / PRPLL](https://github.com/preda/gpuowl)
> by Mihai Preda and George Woltman — see [ATTRIBUTION.md](ATTRIBUTION.md).

---

## Download

**[Get the latest release](https://github.com/sallerk/Mp_p-1_gpu/releases/latest)**
— a prebuilt 64-bit Windows binary. Nothing to install: no Visual Studio, no
CUDA Toolkit, no OpenCL SDK. The only requirement is `OpenCL.dll`, which ships
with your GPU driver. Unzip it anywhere and run it.

The executable is not code-signed, so the first launch will probably raise
SmartScreen — *More info* → *Run anyway*. If you would rather not run a
stranger's binary, building it yourself takes about a minute (see Quick start).

Every version before it is on the
[Releases page](https://github.com/sallerk/Mp_p-1_gpu/releases), newest first,
each with its own notes on what changed. Those older releases are source only —
their binaries are no longer distributed — but every one is still a tag in this
repository, so its exact source is a checkout away (`git checkout v1.4`), and
each builds the same way. [CHANGELOG.md](CHANGELOG.md) is the same history in
one file.

GitHub's green *Code → Download ZIP* button gives you the **whole repository**,
every version at once.

## Speed

Full P-1 runs (stage 1 to B1, stage 2 to B2) at four exponents, same GPU
(NVIDIA RTX 3070), same bounds throughout (B1=100,000, B2=2,000,000). No
checkpoint, no tuning cache. 1.8's M82589933 row is two runs agreeing within
2%; every other row is a single clean run (see the note below the table for
what that trades off).

| exponent | tool | version | setup | stage 1 | stage 2 | final gcd | **total** | vs 1.8 |
|---|---|---|--:|--:|--:|--:|--:|--:|
| M25964951 | **Mp_p-1_gpu** | 1.8 | 19s | 40s | 47s | 3s | **1m50s** | **1.00x** |
| M25964951 | Mp_p-1_gpu | 1.7 | 19s | 41s | 46s | 41s | 2m28s | 1.35x |
| M25964951 | [PrMers](https://github.com/cherubrock-seb/PrMers) | v99.95 | — | 1m03s | 1m38s | *(in stage)* | 2m41s | 1.46x |
| M43112609 | **Mp_p-1_gpu** | 1.8 | 18s | 1m09s | 1m22s | 6s | **2m55s** | **1.00x** |
| M43112609 | Mp_p-1_gpu | 1.7 | 17s | 1m10s | 1m25s | 1m19s | 4m11s | 1.43x |
| M43112609 | PrMers | v99.95 | — | 1m38s | 2m26s | *(in stage)* | 4m05s | 1.40x |
| M57885161 | **Mp_p-1_gpu** | 1.8 | 20s | 1m34s | 2m14s | 8s | **4m18s** | **1.00x** |
| M57885161 | Mp_p-1_gpu | 1.7 | 20s | 1m34s | 2m11s | 2m09s | 6m14s | 1.45x |
| M57885161 | PrMers | v99.95 | — | 3m15s | 4m57s | *(in stage)* | 8m12s | 1.91x |
| M82589933 | **Mp_p-1_gpu** | 1.8 | 22s | 1m35s | 2m09s | 12s | **4m19s** | **1.00x** |
| M82589933 | Mp_p-1_gpu | 1.7 | 26s | 1m36s | 2m11s | 3m03s | 7m16s | 1.68x |
| M82589933 | PrMers | v99.95 | 1m06s | 3m22s | 4m57s | *(in stage)* | 9m25s | 2.18x |
| M82589933 | [gpuowl](https://github.com/preda/gpuowl) | v7.5 | — | — | — | — | *no P-1 support* | n/a |

**The gap is the gcd, and it grows with the exponent.** 1.8's stage 1/stage 2
match 1.7's to within a couple of seconds at every size — both run the same
GPU kernels, picking the same transform shape. The only thing that changed is
`gcd CPU`: this project's own half-GCD (1.7) vs. GMP's `mpz_gcd` (1.8), which
finishes in single-digit-to-teens seconds regardless of exponent while 1.7's
own gcd cost scales with it — 41s at 26M digits, 3m03s at 82.6M.
That is the entire 1.35x-1.68x gap against 1.7. Stage-2 accumulator res64
matched at every exponent both versions were run at, confirming the gcd swap
changed nothing about correctness.

**PrMers's gap widens for the same reason, on its side.** Its own gcd is
folded into the stage figures above (marked *in stage*) rather than broken
out, but the pattern is the same: 1.46x at 26M digits, 2.18x at 82.6M.
PrMers's `setup` column is only populated at M82589933, carried over from an
earlier, more tightly controlled measurement (two runs per tool, agreed
within 10%, PrMers's own startup isolated from its stage-1 gcd); the other
three exponents were each measured once, so PrMers's setup there is folded
into its `stage 1` figure instead of broken out — a `—` means "not isolated
this time", not "zero". gpuowl has no P-1 task type in this build (checked
against its own `Worktodo.cpp` parser), so it appears only at M82589933, as a
non-comparison rather than a fourth contender.

Full methodology, PrMers's V-trace stage 2, and older-version numbers are in
[MANUAL.md](MANUAL.md) and [CHANGELOG.md](CHANGELOG.md).

## Quick start

Build from source — Visual Studio 2019/2022 with "Desktop development with
C++" (or take the prebuilt binary from Download above). `build.bat` also
needs GMP installed once via the vcpkg that ships with Visual Studio; it
prints the exact command if it's missing. See [MANUAL.md](MANUAL.md) for
that one-time step:

```bash
build.bat
```

```bash
Mp_p-1_gpu.exe --selftest
Mp_p-1_gpu.exe
```

Add an exponent to `worktodo.txt` (one per line) and edit `config.txt` to say
how to work it; every setting is documented inline.

[MANUAL.md](MANUAL.md) is the manual: `config.txt` keys, the
worktodo queue, the command line, reading the output, results format,
resuming and raising bounds, and self-tests. Read it before running a real
job.

## Requirements

| | |
|---|---|
| GPU | any OpenCL 1.2 GPU — NVIDIA, AMD or Intel. Developed on an RTX 3070. |
| Runtime | `OpenCL.dll`, which **ships with your GPU driver**. Nothing to install. |
| OS | the prebuilt binaries are 64-bit Windows. The MSVC port took real work — see [ATTRIBUTION.md](ATTRIBUTION.md) — so a Linux build is not a given; no Makefile or CI exists yet. |

The `.exe` is statically linked and resolves OpenCL from the driver at run
time, so it plus `config.txt` are all you need to copy to another machine —
including one with a different GPU vendor.

## Scope and limitations

- **P-1 only implements factoring**, not PRP or Lucas–Lehmer — for that, see
  [PRPLL / gpuowl](https://github.com/preda/gpuowl), which this is built on.
- **P+1 is a secondary mode**, lower yield per unit of GPU time than P-1. It
  exists because no GIMPS *GPU* tool offers P+1 at all; run it in earnest only
  once P-1 has been tried. It has its own B1 and B2 model now, but still
  shares P-1's stage-2 pairing shape (a GPU-memory budget decision, not a cost
  gap), and (unlike P-1) does not yet reuse a completed stage-2 walk when you
  raise B2.
- **Windows / MSVC only.** No Makefile, no Linux build, no CI.
- **No PrimeNet automation by this program itself** — no networking of any
  kind. As of 1.8 it understands the worktodo/results shape
  [AutoPrimeNet](https://github.com/tdulcet/AutoPrimeNet) reads and writes
  (assignment IDs, known factors), so running that alongside this program
  covers the fetch/submit half; without it, results are written for you to
  upload manually. See `MANUAL.md` for the details.

The full, per-version list is in each version's own README.

## Name

"M_p, P-1, GPU". It implements P-1 and P+1 — it is **not** an ECM
implementation.
