// Copyright (C) Mp_p-1_gpu
//
// Stage-2 checkpointing, and the completed-run record a B2 extension reuses.
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
//
// The same file, with `complete` set, outlives the run: the accumulator over
// (b1, b2] is exactly the seed a later run needs to extend to a larger B2,
// since the accumulator is a plain product and
//
//   acc(b1, b2new]  ==  acc(b1, b2old] * acc(b2old, b2new]
//
// The pairing shape does NOT have to match across the two: each range gets its
// own plan and the product does not care. What must match is the exponent, B1,
// and the stage-1 residue the accumulator was built from -- hence `xRes64`.

#pragma once

#include "common.h"

#include <string>
#include <vector>

// Bump when the meaning of the accumulator changes: a different pairing scheme,
// a different A/S recurrence, anything that makes an old file wrong rather than
// merely old.
//   1: position + acc/A/S                          (1.0)
//   2: adds complete, fromB2 and xRes64            (1.1, B2 extension)
inline constexpr u32 STAGE2_FORMAT_VERSION = 2;

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

  // The walk for this (b1, b2] finished and `acc` is final. A partial record is
  // only ever a resume point; only a complete one may seed an extension.
  bool complete = false;

  // 0 when the walk covered all of (b1, b2]. Otherwise the B2 whose completed
  // accumulator seeded it, i.e. this walk only covered (fromB2, b2].
  //
  // This is a COMPARED field on resume, not decoration. A checkpoint written by
  // a from-scratch walk to b2 has the same (exponent, b1, b2, d, w) and the same
  // position as one written by an extension, but its accumulator is missing
  // every term below fromB2. Resuming one against the other loses those primes
  // silently.
  u64 fromB2 = 0;

  // res64 of the stage-1 residue this accumulator was built from. Binds it to
  // the actual x, so a changed E construction can never be extended across.
  u64 xRes64 = 0;

  Words acc;                // the accumulated product
  Words a;                  // A_m = x^((m*D)^2) at this m
  Words s;                  // S_m = x^(D^2 * (2m+1)) at this m
};

std::string defaultStage2Path(u32 exponent, u64 b1, u64 b2);

// Written via a temp file and rename, so an interrupted write cannot leave a
// corrupt checkpoint in place of a good one.
bool saveStage2(const std::string& path, const Stage2State& s, std::string& err);

// Strict resume: returns false (with a reason) if the file is missing, corrupt,
// or describes a different walk. `want` supplies the values that must match --
// every one of exponent, b1, b2, d, w, fromB2 and (when non-zero) xRes64.
bool loadStage2(const std::string& path, const Stage2State& want, Stage2State& out,
                std::string& err);

// Load a COMPLETED record to seed an extension. Deliberately looser than
// loadStage2: B2 and the pairing shape are whatever the earlier run used and are
// read out rather than checked, because neither affects what the accumulator
// means. exponent, b1, xRes64 and the format version are still hard requirements,
// and `complete` must be set -- a partial accumulator is a resume point, never a
// seed.
bool loadCompletedStage2(const std::string& path, u32 exponent, u64 b1, u64 xRes64,
                         Stage2State& out, std::string& err);

// The B2 values with a stage-2 save file on disk for this (exponent, b1),
// largest first. Mirrors findCompletedB1 for stage 1: the naming convention
// lives here, the policy lives in PM1.cpp.
std::vector<u64> findStage2Saves(u32 exponent, u64 b1);
