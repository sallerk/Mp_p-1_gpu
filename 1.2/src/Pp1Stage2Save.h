// Copyright (C) Mp_p-1_gpu
//
// P+1 stage-2 checkpointing.
//
// A separate file and format from P-1's Stage2Save.h, on the SAME rule that
// keeps Stage2Save.h separate from stage-1's Save.h: P-1's stage-2 format is
// load-bearing and already proven, and folding a P+1-only field (`seed`) into
// it would force a magic bump that invalidates every existing P-1 stage-2
// checkpoint on disk for a change that has nothing to do with P-1.
//
// Resuming needs the position in the plan plus the three live residues. The
// plan itself is not stored: it is a pure function of (b1, b2, d, w), and
// every one of those -- plus `seed`, since P+1 tries several seeds and each
// has its own V_E -- is in the header and compared on load. A mismatch means
// the residues describe a DIFFERENT walk, which would not crash -- it would
// quietly compute the wrong product and lose factors -- so each is a hard
// reject, never adapted.
//
// `fromB2` exists for schema parity with Stage2State but is never set in this
// version: B2 extension for P+1 (seeding a larger-B2 walk from a completed
// smaller one, the way P-1's does) is out of scope here -- P+1 doesn't even
// have its own B1 yet (it borrows P-1's), and that gap is the one worth
// closing first. A future version can wire it up without another format bump.

#pragma once

#include "common.h"

#include <string>
#include <vector>

inline constexpr u32 PP1_STAGE2_FORMAT_VERSION = 1;

struct Pp1Stage2State {
  u32 exponent = 0;
  u64 b1 = 0;
  u64 b2 = 0;
  u32 d = 0;
  u32 w = 0;

  // Which P+1 seed this walk belongs to. A hard-reject compared field: a
  // seed-5 walk resumed against seed-3's checkpoint would not crash, it would
  // silently compute a meaningless accumulator, since the two seeds' V_E
  // residues are unrelated values.
  u32 seed = 0;

  u32 version = PP1_STAGE2_FORMAT_VERSION;

  // Position: the next slot to process is (m, jIdx). Everything strictly
  // before it in plan order is already folded into acc.
  u64 m = 0;
  u64 jIdx = 0;
  u64 done = 0;

  // The walk for this (b1, b2] finished and `acc` is final.
  bool complete = false;

  // Present for schema parity with P-1's Stage2State; always 0 in this
  // version. See the file header for why B2 extension is out of scope here.
  u64 fromB2 = 0;

  // res64 of y1 = V_E(seed,1), the stage-1 residue this accumulator was built
  // from. P+1's analogue of Stage2State's xRes64, checked on every resume.
  u64 yRes64 = 0;

  Words acc;                // the accumulated product
  Words a;                  // A_curr = V_(m*D)(y1,1) at this m
  Words s;                  // A_prev = V_((m-1)*D)(y1,1) at this m
};

std::string defaultPp1Stage2Path(u32 exponent, u64 b1, u64 b2, u32 seed);

// Written via a temp file and rename, so an interrupted write cannot leave a
// corrupt checkpoint in place of a good one.
bool savePp1Stage2(const std::string& path, const Pp1Stage2State& s, std::string& err);

// Strict resume: returns false (with a reason) if the file is missing,
// corrupt, or describes a different walk. `want` supplies the values that
// must match -- every one of exponent, b1, b2, d, w, seed, fromB2 and (when
// non-zero) yRes64.
bool loadPp1Stage2(const std::string& path, const Pp1Stage2State& want,
                   Pp1Stage2State& out, std::string& err);

// The B2 values with a stage-2 save file on disk for this (exponent, b1,
// seed), largest first. Mirrors findStage2Saves; unused until P+1 gets B2
// extension, kept for symmetry and because the directory scan itself is
// useful for diagnostics even without it.
std::vector<u64> findPp1Stage2Saves(u32 exponent, u64 b1, u32 seed);
