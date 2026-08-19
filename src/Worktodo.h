// Copyright (C) Mp_p-1_gpu
//
// worktodo.txt: the exponent queue. One exponent per line, top to bottom,
// processed in order; '#' starts a comment, blank lines are ignored -- same
// convention as config.txt. Every entry runs under the SAME config.txt
// settings (method, b1/b2, bias, stages, ...); only the exponent (and its
// factoredTo default, when auto) varies per entry -- there is no per-line
// override syntax.

#pragma once

#include "common.h"

#include <string>
#include <vector>

struct WorktodoEntry {
  u32 exponent = 0;
  u32 lineNo = 0;   // physical line in the file, used to remove exactly this entry
};

// Every valid entry, in file order. A missing file is an EMPTY queue, not an
// error -- there is nothing malformed about not having started a queue yet.
// A malformed non-comment, non-blank line is a hard error (out is left in
// whatever state it reached, but should not be used) -- same fail-fast
// philosophy as loadConfig: a typo should stop the run, not silently vanish
// from the queue.
bool loadWorktodo(const std::string& path, std::vector<WorktodoEntry>& out,
                  std::string& err);

// Removes exactly this entry's line, atomically (temp file + rename). Before
// removing, re-reads that physical line and confirms it still parses to the
// same exponent -- guards against the file having been hand-edited between
// load and consume. On any mismatch or I/O failure, returns false with err
// set and the file is left untouched.
bool consumeWorktodoEntry(const std::string& path, const WorktodoEntry& entry,
                          std::string& err);

// Self-test: round-trip parsing (entries, comments, blank lines, order),
// consumeWorktodoEntry correctness, malformed-line and missing-file handling,
// and the no-trailing-newline-on-the-last-line edge case. No GPU.
int runWorktodoTests();
