// Copyright (C) Mp_p-1_gpu
//
// results.txt: one JSON object per line, the shape PrimeNet accepts and
// Prime95 writes, so the file can be uploaded to mersenne.org's "Manual
// Results" page as-is and read side by side with Prime95's own.
//
// This lives apart from main.cpp for one reason: the exact set of fields, and
// their order, is a compatibility contract with another program, and it has
// been got wrong twice -- once by emitting "seed" where Prime95 writes
// "start", and once by emitting "stage2-fft-length", which turns out to be a
// polymult-only field. Both were caught by reading Prime95's source, not by
// this program noticing. In its own translation unit the writer can be called
// by a self-test, which is what runResultsTests does.

#pragma once

#include "common.h"

#include <vector>

struct Config;
struct FoundFactor;
struct Pp1Start;

// Appends one line to cfg.resultsFile describing the outcome of one job.
//
//   worktype  "P-1" or "P+1", verbatim into the line
//   b1, b2    b2 is omitted when it does not exceed b1, which is how a
//             stage-1-only result is spelled
//   factors   what the gcd yielded. Primes that divide M_p are reported as
//             status "F"; anything else is status "C" and is not submittable,
//             recorded so the run is not silently lost. Empty means "NF".
//   start     P+1's rational starting point, null for P-1
//   stage2D   the stage-2 pairing modulus, 0 when stage 2 did not produce
//             this result (which is how Prime95 decides to omit "d" too)
void writeResultJson(const Config& cfg, const char* worktype, u64 b1, u64 b2,
                     const std::vector<FoundFactor>& factors, const Pp1Start* start,
                     u32 stage2D);

// Self-test: the field set and field ORDER of every line shape this program
// can emit, checked against Prime95's. No GPU.
int runResultsTests();
