# Mp_p-1_gpu

GPU **P-1** and **P+1** factoring of Mersenne numbers `M_p = 2^p - 1`, both with
two stages. Windows / MSVC, no external dependencies — no GMP, no CUDA or
OpenCL SDK.

> **Licence: GPLv3.** Derived from [gpuowl / PRPLL](https://github.com/preda/gpuowl)
> by Mihai Preda and George Woltman — see [ATTRIBUTION.md](ATTRIBUTION.md).

---

## Download

The current version has its own release, with a **prebuilt 64-bit Windows
binary** inside — you do not need Visual Studio, a CUDA Toolkit or an OpenCL
SDK to run it. The only requirement is `OpenCL.dll`, which ships with your GPU
driver. Older versions' prebuilt binaries are no longer distributed; their
source is still in this repository (see below) and builds the same way.

| version | download | |
|---|---|---|
| **1.7** | [Mp_p-1_gpu-1.7-win64.zip](https://github.com/sallerk/Mp_p-1_gpu/releases/download/v1.7/Mp_p-1_gpu-1.7-win64.zip) | current. Stage 2's `T_j` table is built by a recurrence instead of an exponentiation per entry — one multiply per unit of `j` rather than `2*log2(j)` per entry. Results are bit-for-bit unchanged; it is setup cost that stopped being paid, so the win scales with table size: 27% off stage 2 at a 7-digit exponent, a few percent at a wavefront one where GPU memory caps the table. Note that stage 2 is only on a run's critical path when it outlasts the stage-1 gcd it now runs alongside — see Speed below. Includes everything in 1.6. |
| 1.6 | [release notes](https://github.com/sallerk/Mp_p-1_gpu/releases/tag/v1.6) | **Fixed a correctness bug**: resuming an interrupted stage 1 (Ctrl-C, a crash, a reboot) reprocessed one exponent bit, silently corrupting the residue — present in every version before it, 1.5 included. If you've ever resumed a stage-1 run under 1.5 or earlier, treat its result as unverified. Also: stage 1 checkpoints immediately on Ctrl-C now instead of only periodically, Ctrl-C is honored during the gcd phases (it previously was not), the stage-1 gcd runs alongside stage 2 instead of after it, and the schoolbook multiply base case is faster. |
| 1.5 | [release notes](https://github.com/sallerk/Mp_p-1_gpu/releases/tag/v1.5) | **binary withdrawn — contains the checkpoint-resume bug fixed in 1.6**, see that row. Otherwise: the `gcd CPU` phase that ends each stage is about 2.4x faster (parallel Toom-Cook-3, a thread-local allocator after profiling showed the gcd was allocator-bound, and Toom-Cook-4). Nothing about running it changed. |
| 1.4 | *(source only)* | exponents come from `worktodo.txt`, a queue processed in order, instead of a single `exponent =` in `config.txt`. Kept as-is; see [1.4/README.md](1.4/README.md). |
| 1.3 | [release notes](https://github.com/sallerk/Mp_p-1_gpu/releases/tag/v1.3) | P+1 now picks its own B1 — it was silently borrowing P-1's, optimised for the wrong smoothness target, and the auto-chosen value is typically much smaller now. Plus display fixes. |
| 1.2 | [release notes](https://github.com/sallerk/Mp_p-1_gpu/releases/tag/v1.2) | P+1 gained stage 2 — catches a factor whose `q+1` is B1-smooth apart from one prime in `(B1,B2]`, the way P-1's stage 2 already did for `q-1`. Kept as-is. |
| 1.1 | [release notes](https://github.com/sallerk/Mp_p-1_gpu/releases/tag/v1.1) | raising B2 on a finished P-1 stage 2 reuses the completed accumulator and walks only the new range, instead of starting over. Kept as-is. |
| 1.0 | [release notes](https://github.com/sallerk/Mp_p-1_gpu/releases/tag/v1.0) | first release. GPU P-1 with both stages, P+1 stage 1, automatic bound selection, checkpoint/resume, PrimeNet-format JSON results. Kept as-is. |

The [releases page](https://github.com/sallerk/Mp_p-1_gpu/releases) has the
full history; 1.0–1.2 there have their changelog notes but no attached binary.
GitHub's green *Code → Download ZIP* button gives you the **whole
repository**, every version at once.

Each release is also its own directory in this repository, kept exactly as it
shipped — clone it and run `build.bat` (Visual Studio 2019/2022, "Desktop
development with C++") to get a binary for any version. What changed between
them is in [CHANGELOG.md](CHANGELOG.md).

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
C++", nothing else (or take the prebuilt 1.7 binary from Download above):

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
