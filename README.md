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
tuned kernel options, no remembered FFT verdicts. Each tool chooses its own
transform, because on this side that choice is most of what changed in 1.7.
Two runs each; the best of the two is shown, and the pair agreed to within
10%:

| tool | version | setup | stage 1 | stage 2 | final gcd | **total** | relative |
|---|---|--:|--:|--:|--:|--:|--:|
| **Mp_p-1_gpu** | 1.7 | 26s | 1m36s | 2m11s | 3m08s | **7m21s** | **1.00x** |
| [PrMers](https://github.com/cherubrock-seb/PrMers) | v99.95 | 1m06s | 3m22s | 4m57s | *(in stage 2)* | 9m25s | 1.28x slower |
| Mp_p-1_gpu | 1.6 | 10s | 2m02s | 2m47s | 2m56s | 7m55s | *(different sitting — see note)* |
| [gpuowl](https://github.com/preda/gpuowl) | v7.5 | — | — | — | — | *(no P-1 support)* | n/a |

1.7 and PrMers were re-measured together in a second sitting, machine idle
and cooler beforehand; 1.6's row is carried over from the original (warmer)
sitting and is not read against these two as a precise ratio — only the
1.7-vs-PrMers comparison above is same-conditions. Between the two sittings,
PrMers came out **faster** (596s → 565s, −5.6%, confirmed across two rounds
each 3s apart) while 1.7 came out slightly **slower** (427s → 441s, +3.3%)
despite the cooler machine — the FFT this run actually measured *faster*
per-iteration (634 vs. 675 us/it), so the extra time is not the transform;
it landed in the CPU-only final gcd (183s → 188s) and stage-1 gcd (3m18s →
3m09s, still overlapping stage 2), which is exactly where background load
from other work running on the same machine during this second sitting would
show up. Read the 1.7 numbers below with that caveat; PrMers's improvement,
confirmed twice, is the more trustworthy of the two deltas.

**1.6 against 1.7 is almost entirely which transform gets picked** (from the
original, same-sitting measurement of both). 1.6 takes the first candidate
that passes its correctness check, in catalog order, and lands on
`4:256:16:256:101`. 1.7 times candidates and keeps the fastest, which there
was `1:512:8:256:101` — measured 675 against 887 us/it in that run. That one
decision is worth ~25% on both GPU stages. Both versions produced the
identical stage-2 accumulator (`acc res64 5f12b67e85c663c3`) — which the new
1.7 re-run reproduced exactly as well — the check that matters when the
transform underneath changes.

**Reading the columns.** Every row adds across to its own total, which is why
`setup` is a column: without it the rows look like they are missing a minute.
The final gcd needs the stage-2 accumulator, so stage 1, stage 2 and the final
gcd run in series — 26 + 1m36 + 2m11 + 3m08 = 7m21s for 1.7, matching the wall
clock.

What `setup` holds differs by tool. For 1.7 it is almost entirely **choosing
the transform**: 23.2s timing candidates (22.8s in the other round), plus
device init and process start. For PrMers it is startup and kernel build plus
its stage-1 gcd, which runs *after* the stage-1 timer stops and before stage 2
begins — unchanged between sittings, since it is not compute-bound the same
way the GPU stages are. Its stage-2 gcd, by contrast, falls inside the
stage-2 figure, hence *(in stage 2)*. So the two tools' gcd work is not in
comparable columns and only the totals compare cleanly.

There is also a *second* gcd on this side — stage 1's — which since 1.6 runs
alongside stage 2. The program prints its duration (3m09s in this 1.7 run)
but it is not a column and does not add: it overlaps stage 2 and the final
gcd, both of which outlast it.

**That second gcd is no longer free, and 1.7 is why.** In 1.6 it takes about as
long as stage 2 (2m57s against 2m47s), so the two end together and the final
gcd has the CPU to itself. 1.7's stage 2 finishes earlier while the stage-1
gcd still needs about as long, so the final gcd now starts while the stage-1
gcd is still running and the two compete for the same worker threads — in
this run the two land within a second of each other (3m09s against 3m08s).
Some of what 1.7 wins on the GPU is handed back here.

Mp_p-1_gpu is well ahead on both GPU stages — 2.1x on stage 1 and 2.3x on
stage 2 — but gives much of it back in the `gcd CPU` phase, which PrMers runs
through GMP and this runs through its own big-integer code. The final gcd is
the largest of its own four phases at 3m08s, longer than either of its own
GPU stages. **The gcd, not the pairing, is where the remaining headroom is.**

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
C++", nothing else (or take the prebuilt binary from Download above):

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
