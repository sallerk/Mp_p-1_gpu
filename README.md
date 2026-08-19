# Mp_p-1_gpu

GPU **P-1** and **P+1** factoring of Mersenne numbers `M_p = 2^p - 1`, both with
two stages. Windows / MSVC, no external dependencies — no GMP, no CUDA or
OpenCL SDK.

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
their binaries are no longer distributed — but the source of every version is
in this repository, one directory per release, kept exactly as it shipped, and
each builds the same way. [CHANGELOG.md](CHANGELOG.md) is the same history in
one file.

GitHub's green *Code → Download ZIP* button gives you the **whole repository**,
every version at once.

## Speed

A full P-1 run — stage 1 to B1, then stage 2 to B2 — on the same GPU (NVIDIA
RTX 3070), `M82589933` (8 digits, current GIMPS wavefront territory), and the
same bounds (B1=100,000, B2=2,000,000) for every tool. Every row was measured
back to back in one sitting under one FFT tuning, from no checkpoint, so they
are comparable to each other:

| tool | version | stage 1 | stage 2 | final gcd | **total** | relative |
|---|---|--:|--:|--:|--:|--:|
| **Mp_p-1_gpu** | 1.7 | 1m34s | 2m11s | 3m05s | **6m51s** | **1.00x** |
| Mp_p-1_gpu | 1.6 | 1m36s | 2m18s | 3m13s | 7m07s | 1.04x slower |
| [PrMers](https://github.com/cherubrock-seb/PrMers) | v99.95 | 3m14s | 4m56s | *(included)* | 8m10s | 1.19x slower |
| [gpuowl](https://github.com/preda/gpuowl) | v7.5 | — | — | — | *(no P-1 support)* | n/a |

**Reading the three columns.** The final gcd needs the stage-2 accumulator, so
stage 1, stage 2 and the final gcd run in series and simply add up: 1m34 +
2m11 + 3m05 = 6m50s against 6m51s observed for 1.7. There is a *second* gcd —
stage 1's — which since 1.6 runs alongside stage 2; it took 3m13s here and
finished well before the final gcd did, so it costs nothing and is not a
column. That overlap is what 1.6 bought. It also means a second saved in
stage 2 is a second saved overall, which is where 1.7's stage-2 work shows up.

Mp_p-1_gpu is ahead on both GPU stages — about 2.1x on stage 1 and 2.3x on
stage 2 — but gives most of it back in the `gcd CPU` phase, which PrMers runs
through GMP and this runs through its own big-integer code. The final gcd is
the single largest line in the table at 3m05s, longer than either GPU stage.
**The gcd, not the pairing, is where the remaining headroom is.**

**gpuowl has no P-1 factoring in this build.** Its `-B1`/`-B2` flags are
documented but, checked directly against its `Worktodo.cpp`, its worktodo
parser only recognises `PRP=`, `Test=`/`DoubleCheck=` (Lucas-Lehmer), and
`Cert=` assignments — no P-1 task type exists to invoke. Passing `-B1 -B2`
alongside `-prp` silently runs a plain PRP test instead; confirmed by running
it and watching it sail past the P-1 iteration count into millions of PRP
iterations. gpuowl is a primality tool (PRP/LL) here, not a factoring tool —
it isn't in this table because there's nothing P-1-shaped to time.

PrMers's stage 2 ("V-trace"/"Pair95") is a genuinely different construction —
a Lucas trace `V_n = H^n + H^-n` instead of the `x^(j^2)` table shared by
Mp_p-1_gpu and its gpuowl/PRPLL ancestor. It was scoped for 1.7 and not
adopted: its pairing buys candidate partners at exactly the same price per
unit of GPU memory as the existing `stage2_w` window, and porting it to P-1
would need a modular inverse costing more than the stage 2 it accelerates.
The one part that did transfer — building the table by a recurrence rather
than an exponentiation per entry — is 1.7's speedup. There is a fuller
write-up in [1.7/README.md](1.7/README.md).

Take this as one data point, not a definitive verdict — a different exponent,
GPU, bounds, or either tool's own tuning could shift these numbers. The
exponent matters more than it looks: at `M5378909`, a quarter of this size,
the ranking **reverses** and PrMers is about 2.2x faster (25s against 54s,
both at stock settings). This engine descends from gpuowl and is at its best
on wavefront-sized transforms; at 262144 words it does not fill the GPU.

## Quick start

Build from source — Visual Studio 2019/2022 with "Desktop development with
C++", nothing else (or take the prebuilt binary from Download above):

```bash
cd 1.7
build.bat
```

```bash
Mp_p-1_gpu.exe --selftest
Mp_p-1_gpu.exe
```

Add an exponent to `worktodo.txt` (one per line) and edit `config.txt` to say
how to work it; every setting is documented inline.

The [1.7/README.md](1.7/README.md) is the manual: `config.txt` keys, the
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
- **No PrimeNet automation.** `worktodo.txt` is a local queue, not a synced
  one — results are written for you to upload manually.

The full, per-version list is in each version's own README.

## Name

"M_p, P-1, GPU". It implements P-1 and P+1 — it is **not** an ECM
implementation.
