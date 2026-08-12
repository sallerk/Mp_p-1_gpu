# Mp_p-1_gpu

GPU **P-1** and **P+1** factoring of Mersenne numbers `M_p = 2^p - 1`, both with
two stages. Windows / MSVC, no external dependencies — no GMP, no CUDA or
OpenCL SDK.

> **Licence: GPLv3.** Derived from [gpuowl / PRPLL](https://github.com/preda/gpuowl)
> by Mihai Preda and George Woltman — see [ATTRIBUTION.md](ATTRIBUTION.md).

---

## Download

Each version has its own release, with a **prebuilt 64-bit Windows binary**
inside — you do not need Visual Studio, a CUDA Toolkit or an OpenCL SDK to run
it. The only requirement is `OpenCL.dll`, which ships with your GPU driver.

| version | download | |
|---|---|---|
| **1.2** | [**Mp_p-1_gpu-1.2-win64.zip**](https://github.com/sallerk/Mp_p-1_gpu/releases/download/v1.2/Mp_p-1_gpu-1.2-win64.zip) | current. P+1 now has stage 2 — catches a factor whose `q+1` is B1-smooth apart from one prime in `(B1,B2]`, the way P-1's stage 2 already did for `q-1`. Start here — see [1.2/README.md](1.2/README.md). |
| 1.1 | [Mp_p-1_gpu-1.1-win64.zip](https://github.com/sallerk/Mp_p-1_gpu/releases/download/v1.1/Mp_p-1_gpu-1.1-win64.zip) | raising B2 on a finished P-1 stage 2 reuses the completed accumulator and walks only the new range, instead of starting over. Kept as-is. |
| 1.0 | [Mp_p-1_gpu-1.0-win64.zip](https://github.com/sallerk/Mp_p-1_gpu/releases/download/v1.0/Mp_p-1_gpu-1.0-win64.zip) | first release. GPU P-1 with both stages, P+1 stage 1, automatic bound selection, checkpoint/resume, PrimeNet-format JSON results. Kept as-is. |

All of them are on the [releases page](https://github.com/sallerk/Mp_p-1_gpu/releases).
GitHub's green *Code → Download ZIP* button gives you the **whole repository**,
every version at once — the links above are the way to get one version on its
own.

Each release is also its own directory in this repository, kept exactly as it
shipped. What changed between them is in [CHANGELOG.md](CHANGELOG.md).

## Quick start

Unzip a release and run it:

```bash
Mp_p-1_gpu.exe --selftest
Mp_p-1_gpu.exe
```

Edit `config.txt` to say what to factor and how; every setting is documented
inline. Or build from source — Visual Studio 2019/2022 with "Desktop
development with C++", nothing else:

```bash
cd 1.2
build.bat
```

The [1.2/README.md](1.2/README.md) is the manual: `config.txt` keys, the
command line, reading the output, results format, resuming and raising
bounds, and self-tests. Read it before running a real job.

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
  once P-1 has been tried. It borrows P-1's bounds and pairing shape, and
  (unlike P-1) does not yet reuse a completed stage-2 walk when you raise B2.
- **Windows / MSVC only.** No Makefile, no Linux build, no CI.
- **One exponent per run**, from `config.txt`. No worktodo queue, no PrimeNet
  automation — results are written for you to upload manually.

The full, per-version list is in each version's own README.

## Name

"M_p, P-1, GPU". It implements P-1 and P+1 — it is **not** an ECM
implementation.
