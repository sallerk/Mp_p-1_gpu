// Copyright (C) Mp_p-1_gpu
//
// Stage-2 checkpointing.
//
// A separate file and a separate format from the stage-1 checkpoint in Save.h,
// on purpose. Stage-1 resume is load-bearing and already proven; folding new
// fields into its header would force a magic bump and invalidate every existing
// stage-1 checkpoint on disk for no benefit.
//
// Resuming needs the position in the plan plus the three live residues. The plan
// itself is not stored: it is a pure function of (b1, b2, d, w), and every one
// of those is in the header and compared on load. A mismatch means the residues
// describe a DIFFERENT walk, which would not crash -- it would quietly compute
// the wrong product and lose factors -- so each is a hard reject, never adapted.

#pragma once

#include "common.h"

#include <string>

// Bump when the meaning of the accumulator changes: a different pairing scheme,
// a different A/S recurrence, anything that makes an old file wrong rather than
// merely old.
inline constexpr u32 STAGE2_FORMAT_VERSION = 1;

struct Stage2State {
  u32 exponent = 0;
  u64 b1 = 0;
  u64 b2 = 0;
  u32 d = 0;
  u32 w = 0;
  u32 version = STAGE2_FORMAT_VERSION;

  // Position: the next slot to process is (m, jIdx). Everything strictly before
  // it in plan order is already folded into acc.
  u64 m = 0;
  u64 jIdx = 0;
  u64 done = 0;             // multiplies completed, for progress continuity

  Words acc;                // the accumulated product
  Words a;                  // A_m = x^((m*D)^2) at this m
  Words s;                  // S_m = x^(D^2 * (2m+1)) at this m
};

std::string defaultStage2Path(u32 exponent, u64 b1, u64 b2);

// Written via a temp file and rename, so an interrupted write cannot leave a
// corrupt checkpoint in place of a good one.
bool saveStage2(const std::string& path, const Stage2State& s, std::string& err);

// Returns false (with a reason) if the file is missing, corrupt, or describes a
// different walk. `want` supplies the values that must match.
bool loadStage2(const std::string& path, const Stage2State& want, Stage2State& out,
                std::string& err);
