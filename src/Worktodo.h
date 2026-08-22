// Copyright (C) Mp_p-1_gpu
//
// worktodo.txt: the exponent queue. One entry per line, top to bottom,
// processed in order; '#' starts a comment, blank lines are ignored -- same
// convention as config.txt. Every entry runs under the SAME config.txt
// settings (method, bias, stages, ...); a bare integer is this program's own
// shorthand for "factor M_p, everything else from config.txt".
//
// A line can also be a Pfactor= or Pminus1= assignment, the shape AutoPrimeNet
// (https://github.com/tdulcet/AutoPrimeNet) and Prime95 itself write. Such a
// line carries its own k,b,n,c (validated as a Mersenne number: k=1,b=2,c=-1)
// and, for Pminus1=, its own B1/B2, which are used whenever present -- a
// Pfactor= line or a bare exponent has none, and falls back to config.txt's
// own b1/b2 (auto by default). Both may also carry an assignment ID and
// known factors, echoed back into results.txt so AutoPrimeNet's own upload
// step stays consistent with what PrimeNet expects.
// This program never talks to PrimeNet itself; AutoPrimeNet is a separate
// process that reads/writes these same two files.

#pragma once

#include "common.h"

#include <string>
#include <vector>

struct WorktodoEntry {
  u32 exponent = 0;
  u32 lineNo = 0;   // physical line in the file, used to remove exactly this entry

  // Set only by a Pfactor=/Pminus1= line; a bare-integer line leaves every
  // field below at its default, identical to today's behavior.
  std::string aid;                        // "" == none (32 hex chars in the file)

  bool hasAssignedBounds = false;         // true only for Pminus1=
  u64 assignedB1 = 0;
  u64 assignedB2 = 0;

  // Pminus1='s optional B2_start: stage 2 was already walked up to this point
  // elsewhere. Not honored (see Stage2Save.h's extend path -- it needs a
  // local checkpoint to seed from, and there isn't one for externally
  // declared progress); kept only so runOneJob can warn instead of silently
  // redoing work AutoPrimeNet/PrimeNet believe is already covered.
  u64 b2Start = 0;

  bool hasFactoredTo = false;             // how_far_factored (Pfactor) or
                                          // sieve_depth (Pminus1) was present
  u32 factoredTo = 0;                     // bits, floored; same range as Config::factoredTo

  double testsSaved = 0;                  // Pfactor= only; FYI, see Config.h's note
  std::vector<std::string> knownFactors;  // decimal strings, echoed back verbatim
};

// Every valid entry, in file order. A missing file is an EMPTY queue, not an
// error -- there is nothing malformed about not having started a queue yet.
// A malformed non-comment, non-blank line is a hard error (out is left in
// whatever state it reached, but should not be used) -- same fail-fast
// philosophy as loadConfig: a typo (or an unsupported worktodo keyword, e.g.
// PRP= -- this program only factors) should stop the run, not silently
// vanish from the queue.
bool loadWorktodo(const std::string& path, std::vector<WorktodoEntry>& out,
                  std::string& err);

// Removes exactly this entry's line, atomically (temp file + rename). Before
// removing, re-reads that physical line and confirms it still parses to the
// same exponent (and, if the original entry had one, the same aid) --
// guards against the file having been hand-edited, or reissued by
// AutoPrimeNet under a new assignment, between load and consume. On any
// mismatch or I/O failure, returns false with err set and the file is left
// untouched.
bool consumeWorktodoEntry(const std::string& path, const WorktodoEntry& entry,
                          std::string& err);

// Which B1/B2 a queued entry actually runs at: the entry's own assigned
// bounds (Pminus1=) whenever it has them, else whatever config.txt says
// (configuredB1/B2 -- 0 means "auto", same meaning it always has). Pulled
// out of the queue loop in main.cpp so this precedence decision -- changed
// twice already in one session -- has one place to test instead of none.
struct ResolvedBounds { u64 b1, b2; };
ResolvedBounds resolveBounds(const WorktodoEntry& entry, u64 configuredB1, u64 configuredB2);

// Self-test: round-trip parsing (bare entries, Pfactor=/Pminus1= assignments,
// comments, blank lines, order), consumeWorktodoEntry correctness, malformed-
// line and missing-file handling, and the no-trailing-newline-on-the-last-line
// edge case. No GPU.
int runWorktodoTests();
