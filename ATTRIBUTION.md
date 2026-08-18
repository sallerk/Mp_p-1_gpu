# Attribution

Mp_p-1_gpu is licensed under the **GNU General Public License, version 3**
(see [LICENSE](LICENSE)), because it is a derivative work of a GPLv3 project.

This applies to every release: the repository is organized as one subfolder per
version (`1.0/`, `1.1/`, `1.2/`, `1.3/`, ...), and the file paths below are
relative to each version's own `src/` — e.g. `Gpu.{h,cpp}` means
`1.6/src/Gpu.{h,cpp}` in the current release. The provenance split (what is
upstream vs. original) does not change between versions; only the code within
each file does.

## Upstream

**gpuowl / PRPLL** — Mihai Preda and George Woltman
<https://github.com/preda/gpuowl>

Everything that performs modular arithmetic in this repository comes from there:

| Area | Files |
|---|---|
| OpenCL FFT/NTT kernels | all of `src/cl/*.cl`, and `src/bundle.cpp` (generated from them) |
| Transform pipeline, carries, trig | `Gpu.{h,cpp}`, `Trig.cpp`, `TrigBufCache.cpp`, `fftbpw.h` |
| FFT selection and tuning | `FFTConfig.{h,cpp}`, `tune.{h,cpp}`, `TuneEntry.{h,cpp}` |
| OpenCL plumbing | `clwrap.{h,cpp}`, `Queue.{h,cpp}`, `Kernel.{h,cpp}`, `KernelCompiler.{h,cpp}`, `Buffer.h`, `Context.h`, `Event.{h,cpp}`, `AllocTrac.{h,cpp}`, `gpuid.{h,cpp}`, `tinycl.h` |
| Support | `Args.{h,cpp}`, `File.{h,cpp}`, `fs.{h,cpp}`, `log.{h,cpp}`, `timeutil.{h,cpp}`, `common.{h,cpp}`, `state.{h,cpp}`, `Primes.{h,cpp}`, `Profile.{h,cpp}`, `TimeInfo.{h,cpp}`, `md5.cpp`, `sha3.cpp`, `Sha3Hash.h`, `MD5.h`, `Hash.h`, `CycleFile.{h,cpp}`, `Saver.{h,cpp}`, `Signal.{h,cpp}`, `Proof.{h,cpp}`, `Background.h`, `Task.h`, `shared.h`, `typeName.h`, `GpuCommon.h` |
| Probability model | the Dickman-rho table and `pm1Prob` in `Bounds.cpp`, ported from upstream `pm1/pm1.cpp` |

Those files carry `// Copyright (C) Mihai Preda` or
`// Copyright (C) Mihai Preda and George Woltman`. **GPLv3 §5 requires those
notices to be preserved.** They must not be removed or rewritten — including
during any future renaming of this project. Removing them would breach the
licence and make the result undistributable.

Renaming the *project* is unaffected: the name, executable, banner and
configuration files are ours to change. Branding and authorship are separate
things.

## Original to this project

Written for this program, and GPLv3 along with the rest:

| Area | Files |
|---|---|
| Multiprecision integers | `BigInt.{h,cpp}` — Karatsuba, Knuth-D division, Miller-Rabin |
| Subquadratic GCD | `Gcd.{h,cpp}` — recursive half-GCD with parallel cofactor apply |
| P-1 / P+1 driver | `PM1.{h,cpp}` — stage 1, P+1 Lucas ladder, stage-2 driver |
| Stage-2 pairing | `Stage2Plan.{h,cpp}` — Montgomery pairing with a matching window |
| Checkpointing | `Save.{h,cpp}`, `Stage2Save.{h,cpp}`, `Pp1Stage2Save.{h,cpp}` |
| Bounds | the cost model and `chooseBounds` in `Bounds.{h,cpp}`; `pp1Prob` and `choosePP1Bounds` (P+1's own B1/B2 model — no upstream P+1 probability tool exists to port from, unlike `pm1Prob`) |
| Config and tests | `Config.{h,cpp}`, `Selftest.cpp`, `testBigInt.cpp` |
| Exponent queue | `Worktodo.{h,cpp}` (1.4) — `worktodo.txt` parsing and atomic per-entry removal |
| Build | `build.bat`, `tools/` |
| Stage-2 kernels | `subWords` in `src/cl/etc.cl` |
| P+1 stage 2 | `lucasVResidue`, `pp1Stage2`, `renormalize` in `Gpu.{h,cpp}` (no new kernel — built entirely from the existing `subWords`/`modMul`/`squareLL`) |

`main.cpp` is original; `clshim.cpp` is original (it resolves OpenCL entry
points from the driver at run time so no SDK is needed).

## Modifications to upstream files

The upstream sources were ported to MSVC. Changes include: replacing
`__int128`/`__float128` with two-limb emulation, removing GCC-only constructs,
splitting the kernel bundle into <64 KB chunks for MSVC's string-literal limit,
fixing a deadlock in `Queue::finish()`, fixing a missing `else` in
`FFTConfig.cpp`, adding a lower bound to FFT candidate selection, verifying the
tuner's baseline shape before trusting measurements against it, and adding the
`subWords` kernel plus stage-2 entry points to `Gpu`.
