// Copyright (C) Mp_p-1_gpu
//
// Stage-1 checkpointing.
//
// A stage-1 run is hours long, so losing it to a crash, a reboot or a Ctrl-C is
// the single worst failure mode. The state needed to resume is just the current
// residue plus how far through the exponent's bits we are.
//
// The metadata exists to make a WRONG resume impossible. A residue is only
// meaningful for one (exponent, B1, base, E-construction) combination, and
// reusing one across a change is not a crash -- it silently produces a
// plausible answer that misses factors. That is not hypothetical: the stage-1
// exponent gained its 2*p factor partway through this project's development,
// and any residue saved before that change is worthless after it. So every
// field is compared on load and a mismatch is rejected rather than adapted.

#pragma once

#include "common.h"

#include <string>

// Bump when the meaning of the residue changes -- a different base, a different
// E construction, anything that makes an old file wrong rather than merely old.
//   1: E = prod q^floor(log_q B1)                 (WRONG: no 2p; never shipped)
//   2: E = 2*p * prod q^floor(log_q B1)           (current)
// The on-disk header also gained a second residue for P+1; the magic was bumped
// with it so older files are rejected cleanly rather than misparsed.
inline constexpr u32 E_FORMAT_VERSION = 2;

struct SaveState {
  u32 exponent = 0;
  u64 b1 = 0;
  u64 eBits = 0;        // bit length of E; cross-checks that E rebuilt identically
  u64 nextBit = 0;      // next bit index of E to process
  // 0 means the ladder started from base 3 (a fresh stage 1). Non-zero means
  // this is an EXTENSION whose base is the completed residue for that B1, so
  // that file must still exist and validate before this one can be resumed.
  u64 baseB1 = 0;
  bool complete = false;  // stage 1 finished for this (exponent, b1)
  u32 base = 3;
  u32 eVersion = E_FORMAT_VERSION;
  Words residue;
  // P+1's Lucas ladder carries the PAIR (V_k, V_{k+1}); both are needed to
  // resume. Empty for P-1.
  Words residue2;
  u32 seed = 0;         // P+1 seed; 0 for P-1
};

// Default path when none is configured.
std::string defaultSavePath(u32 exponent, u64 b1);

// Written via a temp file and rename, so an interrupted write cannot leave a
// corrupt checkpoint in place of a good one.
bool saveState(const std::string& path, const SaveState& s, std::string& err);

// Returns false (with a reason) if the file is missing, corrupt, or describes a
// different job. `want` supplies the values that must match.
bool loadState(const std::string& path, const SaveState& want, SaveState& out,
               std::string& err);
