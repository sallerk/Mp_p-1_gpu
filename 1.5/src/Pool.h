// Copyright (C) Mp_p-1_gpu
//
// A pooled allocator for Nat::w (see BigInt.h). Every Nat-producing operation
// -- add, sub, every mul tier, shl, shr, divrem -- allocates a fresh backing
// array; at production scale (gcd operands ~1.29M limbs) the recursive
// Toom-3/4/half-GCD call tree does this an enormous number of times. PoolAlloc
// hands out and reclaims those arrays from a size-bucketed freelist instead of
// going through the general-purpose heap allocator each time.
//
// Global, not thread_local: worker threads here come from Parallel.h's
// runParallel, which dispatches via std::async -- on MSVC that typically
// spawns a fresh OS thread per call rather than reusing a warm pool, so a
// thread_local freelist would go cold on every dispatch and be discarded,
// defeating pooling on exactly the hot parallel path (matCompose/matApply/
// mulToom3/mulToom4). A global structure sidesteps that: a block allocated on
// one thread and freed on another just goes back into the same shared bucket.
//
// Buckets are sized to the next power of two, so a block freed by one call is
// immediately usable by any later call asking for up to that many bytes --
// recursive Toom splits and interpolation temporaries produce many distinct
// sizes, not just a couple of threshold values, so a general bucketing scheme
// is needed rather than a couple of fixed sizes. Allocations above
// POOL_MAX_BYTES bypass the pool entirely (plain new/delete) so worst-case
// retained memory stays bounded.

#pragma once

#include "common.h"

#include <atomic>
#include <cstddef>

namespace poolimpl {

inline constexpr size_t POOL_MAX_BYTES = 32u * 1024 * 1024;   // 32 MB cap per block

void* poolAllocate(size_t bytes);
void poolDeallocate(void* p, size_t bytes) noexcept;

} // namespace poolimpl

template<typename T>
struct PoolAlloc {
  using value_type = T;
  using propagate_on_container_move_assignment = std::true_type;
  using is_always_equal = std::true_type;

  PoolAlloc() noexcept = default;
  template<typename U> PoolAlloc(const PoolAlloc<U>&) noexcept {}

  template<typename U> struct rebind { using other = PoolAlloc<U>; };

  T* allocate(size_t n) {
    return static_cast<T*>(poolimpl::poolAllocate(n * sizeof(T)));
  }
  void deallocate(T* p, size_t n) noexcept {
    poolimpl::poolDeallocate(p, n * sizeof(T));
  }
};

template<typename T, typename U>
bool operator==(const PoolAlloc<T>&, const PoolAlloc<U>&) noexcept { return true; }
template<typename T, typename U>
bool operator!=(const PoolAlloc<T>&, const PoolAlloc<U>&) noexcept { return false; }
