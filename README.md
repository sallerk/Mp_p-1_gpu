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

A full P-1 run — stage 1 to B1, then stage 2 to B2 — on the same GPU (NVIDIA
RTX 3070), `M82589933` (8 digits, current GIMPS wavefront territory), and the
same bounds (B1=100,000, B2=2,000,000) for every tool. Every row was measured
from no checkpoint, and **with no tuning of any kind** — no `tune.txt`, no
tuned kernel options, no remembered FFT verdicts. 1.8's row: two runs, best of
the two shown, the pair agreed to within 2% (4m19s and 4m23s). The 1.7,
PrMers and gpuowl rows are carried over unchanged from the prior measurement
(same GPU, same bounds, same no-tuning discipline):

| tool | version | setup | stage 1 | stage 2 | final gcd | **total** | relative |
|---|---|--:|--:|--:|--:|--:|--:|
| **Mp_p-1_gpu** | 1.8 | 22s | 1m35s | 2m09s | 12s | **4m19s** | **1.00x** |
| Mp_p-1_gpu | 1.7 | 26s | 1m36s | 2m11s | 3m03s | 7m16s | 1.68x slower |
| [PrMers](https://github.com/cherubrock-seb/PrMers) | v99.95 | 1m06s | 3m22s | 4m57s | *(in stage 2)* | 9m25s | 2.18x slower |
| [gpuowl](https://github.com/preda/gpuowl) | v7.5 | — | — | — | — | *(no P-1 support)* | n/a |

**1.8 against 1.7 is entirely the gcd.** Both land on the same transform
shape here (`1:256:16:256:101` / `1:512:8:256:101`, 620-650 us/it either way
— which exact candidate wins that pair is a coin flip run to run, a couple of
percent, not the story), and stage 1/stage 2 come in within a second of each
other. What changed is `gcd CPU`: 1.7's own hand-rolled half-GCD takes 3m03s
on this exponent; 1.8 hands the identical computation to GMP's `mpz_gcd` and
gets 12s back — about 15x faster, and the whole reason the run drops from
7m16s to 4m19s. Both versions produced the identical stage-2 accumulator
(`acc res64 5f12b67e85c663c3`), the check that matters when the gcd backend
underneath changes.

**Reading the columns.** Every row adds across to its own total, which is why
`setup` is a column: without it the rows look like they are missing a minute.
The final gcd needs the stage-2 accumulator, so stage 1, stage 2 and the final
gcd run in series — 22 + 1m35 + 2m09 + 12s = 4m19s for 1.8, matching the wall
clock.

What `setup` holds differs by tool. For 1.7 and 1.8 alike it is almost
entirely **choosing the transform** — about 22s timing candidates, plus
device init and process start; neither version's GPU kernels changed here, so
neither version's setup did either. For PrMers it is startup and kernel build
plus its stage-1 gcd, which runs *after* the stage-1 timer stops and before
stage 2 begins; its stage-2 gcd, by contrast, falls inside the stage-2
figure, hence *(in stage 2)*. So the two tools' gcd work is not in comparable
columns and only the totals compare cleanly.

There is also a *second* gcd on this side — stage 1's — which since 1.6 runs
alongside stage 2. The program prints its duration (12s here) but it is not a
column and does not add: it overlaps stage 2 and the final gcd, both of which
outlast it. In 1.7 this second gcd competed with the final gcd for the same
worker threads — both were this project's own multi-threaded half-GCD, and
they would land within seconds of each other, handing back some of what the
GPU stages had won. GMP's gcd is single-threaded, so 1.8 has nothing left to
contend over: both gcds finish in 12s each, a fifth of PrMers's own setup
time alone, and the second gcd is now genuinely free, not just cheaper.

Mp_p-1_gpu is well ahead of PrMers on both GPU stages — 2.1x on stage 1, 2.3x
on stage 2 — and, as of 1.8, no longer gives any of it back in the `gcd CPU`
phase either: at 12s it is now the *smallest* line in the table, not the
largest. **The headroom 1.7's own writeup called out here — "the gcd, not the
pairing, is where the remaining headroom is" — is what 1.8 closes.**

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
write-up in [MANUAL.md](MANUAL.md).

Take this as one data point, not a definitive verdict — a different exponent,
GPU, bounds, or either tool's own tuning could shift these numbers. The
exponent still matters more than it looks: at `M5378909`, a quarter of this
size and at the same bounds, PrMers is **still ahead**, 31s against 42s. But
that gap used to be 2.2x and is now 1.35x, and the reason is the same FFT
change — 1.7 picks `3:256:2:256:101` at 62 us/it where 1.6 took
`2:256:2:256:101` at 118, cutting the run from 65s to 42s. At that size 1.7's
stage 2 is now *faster* than PrMers's (15s against 16s) and its stage 1 close
behind (11s against 10s); what is left of PrMers's lead there is startup and
gcd, and about 9s of 1.7's 42s is the transform search itself, which a repeat
run of the same exponent does not pay again. This engine descends from gpuowl
and is still at its best on wavefront-sized transforms; at 262144 words it
does not fill the GPU.

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
