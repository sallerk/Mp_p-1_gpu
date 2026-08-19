// Copyright (C) Mp_p-1_gpu
//
// Shared thread-budget dispatcher for independent large-operand work: Gcd.cpp's
// matCompose/matApply (8/4 independent matrix sub-multiplies) and BigInt.cpp's
// mulToom3/mulToom4 (5/7 independent evaluation-point sub-products). One global
// budget, not one per file, because these call sites nest -- a mulToom3 call
// happening inside an already-parallel matCompose batch must see the same
// budget matCompose already spent, or the machine gets oversubscribed.

#pragma once

#include "common.h"

#include <functional>
#include <vector>

// Cap the worker threads runParallel will spawn. 0 restores the default (all
// hardware threads).
void setParallelThreads(unsigned n);

// Operands smaller than this are not worth a thread hand-off.
inline constexpr size_t PARALLEL_MIN_LIMBS = 2000;

// Runs `tasks` across worker threads if `worthIt` and the global budget has
// room, otherwise runs them inline on the calling thread. The caller's own
// thread always does its share of the work rather than idling.
void runParallel(std::vector<std::function<void()>> tasks, bool worthIt);
