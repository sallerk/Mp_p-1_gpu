// Copyright (C) Mp_p-1_gpu
//
// config.txt, in the same shape as mersenne_tf's: one key = value per line,
// '#' starts a comment. Numbers may be plain (4000000), separated
// (4_000_000), suffixed (4M, 500K) or "auto" where noted.

#pragma once

#include "common.h"

#include <string>
#include <vector>

// Reported in the banner and in every results line. Bump the version for any
// change that alters what a result MEANS, so a submitted result can be traced
// to the code that produced it.
inline constexpr const char* PROGRAM_NAME = "Mp_p-1_gpu";
inline constexpr const char* PROGRAM_VERSION = "1.8";

// What to do with the console window when the program exits.
//   AUTO   hold only when launched by double-click, i.e. when closing would
//          otherwise wipe the output before it can be read. Running from an
//          existing shell does not pause.
//   ALWAYS always wait for Enter.
//   NEVER  never wait -- the right choice for scripts and schedulers.
enum PauseMode { PAUSE_AUTO = 0, PAUSE_ALWAYS, PAUSE_NEVER };

// Which stages to run.
//   STAGE2_AUTO  let the bounds model decide. It may still decline stage 2: on
//                a small exponent at a low bias, P-1 barely pays for itself and
//                a second stage does not earn its place.
//   STAGE2_ON    always run stage 2, even where the model would not choose it.
//   STAGE2_OFF   stage 1 only. B1 is then optimised on its own, which makes it
//                LARGER than it would be with a stage 2 following -- with no
//                second stage there is nothing else to catch a near-miss, so it
//                pays to push B1 harder.
enum Stage2Mode { STAGE2_AUTO = 0, STAGE2_ON, STAGE2_OFF };

struct Config {
  // p in M_p = 2^p - 1; must be prime. Populated per worktodo.txt entry by
  // the driver, not parsed from config.txt -- see Worktodo.h.
  u32 exponent = 0;

  bool doPM1 = true;            // method = pm1 | pp1 | both
  bool doPP1 = false;
  // P+1 seeds. Each has ~50% chance the Lucas sequence lands in the group that
  // makes it a genuine P+1 rather than a disguised P-1, so several are tried.
  // 2 is excluded: V_n(2,1) == 2 for every n.
  std::vector<u32> pp1Seeds{3, 5, 7};

  u64 b1 = 0;                   // 0 == auto
  u64 b2 = 0;                   // 0 == auto; ignored when stage 2 does not run
  u32 factoredTo = 0;           // TF level in bits; 0 == pick from exponent
  double bias = 1.0;            // value of a factor vs a PRP result

  // How much worse than the cheapest expected cost is still acceptable, as a
  // fraction; the highest-yield bounds inside that band are chosen. The cost
  // minimum is very flat -- a few percent of extra expected cost can nearly
  // double P(factor) -- so the exact argmin is false precision. 0 = strict
  // argmin.
  double boundsTolerance = 0.05;

  int stage2Mode = STAGE2_AUTO; // stages = auto | both | 1

  // Pairing shape for stage 2. 0 == auto, sized to free GPU memory. A larger d
  // or window pairs more primes per multiply but needs more T_j buffers.
  u32 stage2D = 0;
  u32 stage2W = 0;

  std::string fftSpec;          // "" == auto
  int device = -1;              // -1 == auto

  // Exponents to work on, one per line, processed in order -- see Worktodo.h.
  std::string worktodoFile = "worktodo.txt";

  std::string resultsFile = "results.txt";
  // Optional, written into each result line when set. PrimeNet matches results
  // to an account by these; leave blank and they are simply omitted.
  std::string username;
  std::string computerName;
  bool checkpoint = true;
  std::string checkpointFile;   // "" == checkpoint_<exponent>.txt
  u32 checkpointSeconds = 300;

  u32 reportEvery = 1000;       // squarings between progress updates
  int pauseMode = PAUSE_AUTO;   // hold the window open at exit

  // Worker threads for the stage-1 gcd (phase 3). 0 == every hardware thread.
  // The gcd is CPU-only; the GPU is idle throughout, so this is the one knob
  // that speeds that phase up.
  u32 gcdThreads = 0;

  // Check the chosen FFT actually computes correctly before committing hours
  // to it. Costs a couple of seconds. Some configs really do produce wrong
  // results -- see the notes in Selftest.cpp.
  bool verifyFft = true;

  // Reuse completed work instead of redoing it:
  //   stage 1  a smaller completed B1 is raised to this one, saving exactly
  //            oldB1/newB1 of the squarings;
  //   stage 2  a smaller completed B2's accumulator seeds this run, which then
  //            walks only (oldB2, newB2].
  bool extend = true;

  // Instead of exiting when worktodo.txt is empty, wait and re-check it --
  // for running unattended alongside AutoPrimeNet (github.com/tdulcet/
  // AutoPrimeNet), which appends new Pfactor=/Pminus1= lines as it fetches
  // assignments. Off by default: an empty queue exiting cleanly is what a
  // script or a one-shot manual run expects.
  bool waitForWork = false;         // wait_for_work = yes | no
  u32 waitPollSeconds = 5;          // wait_poll_seconds

  // ---- Populated per worktodo.txt entry by the driver, same precedent as
  // exponent/factoredTo above -- NOT read from config.txt itself. Carries a
  // Pfactor=/Pminus1= assignment's identity through to the job-start
  // printout and writeResultJson(), so a submitted result can be matched
  // back to the assignment that produced it and AutoPrimeNet's own upload
  // step sees the known factors it expects. Empty/zero for a bare-exponent
  // entry, or when the assignment simply omitted them.
  std::string aid;
  std::vector<std::string> knownFactors;
  // Pfactor='s tests_saved: a GIMPS PRP-credit concept (PRP tests skipped if
  // a factor turns up) with no equivalent in this program's own cost model --
  // bias means "value of a factor vs. a composite result", not "tests
  // saved" -- so this is carried through only to print as FYI, never fed
  // into chooseBounds/choosePP1Bounds.
  double testsSaved = 0;
  // Pminus1='s optional B2_start (stage 2 already partly covered elsewhere)
  // cannot be honored -- see Worktodo.h's WorktodoEntry::b2Start comment --
  // so runOneJob prints one warning instead of silently ignoring it.
  bool b2StartIgnored = false;
  u64 ignoredB2Start = 0;
};

// Returns false and sets `err` on a malformed file.
bool loadConfig(const std::string& path, Config& cfg, std::string& err);

// Parse a number with optional _ separators and a K/M/G suffix.
bool parseNumber(const std::string& s, u64& out);

// A sensible trial-factoring depth when the user does not supply one. These
// track what GIMPS has actually done at each size, and are only a fallback --
// a wrong value skews the automatic B1.
u32 defaultFactoredTo(u32 exponent);
