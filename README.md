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
| **1.5** | [Mp_p-1_gpu-1.5-win64.zip](https://github.com/sallerk/Mp_p-1_gpu/releases/download/v1.5/Mp_p-1_gpu-1.5-win64.zip) | current. Speed release — the `gcd CPU` phase that ends each stage is about **2.4x faster** (parallel Toom-Cook-3, a thread-local allocator after profiling showed the gcd was allocator-bound, and Toom-Cook-4). Nothing about running it changed. |
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

## Quick start

Build from source — Visual Studio 2019/2022 with "Desktop development with
C++", nothing else (1.4 has no packaged release yet; older versions can be
unzipped and run directly instead):

```bash
cd 1.4
build.bat
```

```bash
Mp_p-1_gpu.exe --selftest
Mp_p-1_gpu.exe
```

Add an exponent to `worktodo.txt` (one per line) and edit `config.txt` to say
how to work it; every setting is documented inline.

The [1.4/README.md](1.4/README.md) is the manual: `config.txt` keys, the
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
