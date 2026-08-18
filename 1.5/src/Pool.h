// Copyright (C) Mp_p-1_gpu
//
// A pooled allocator for Nat::w (see BigInt.h). Every Nat-producing operation
// -- add, sub, every mul tier, shl, shr, divrem -- allocates a fresh backing
// array, and at production scale the gcd does that hundreds of millions of
// times.
//
// This is not a guess about where the time goes. Sampling-profiling the real
// gcd workload (VS DiagnosticsHub CPU sampling, ~1.8M samples per run,
// aggregated by module) put ~45% of ALL process CPU time inside ntdll -- the
// heap -- against ~53% in the program's own arithmetic. The allocator, not
// the multiply, is the single largest cost centre in this workload.
//
// THREAD-LOCAL, and deliberately so. An earlier attempt used one global
// free list with a mutex per size bucket; it was correct but 3.4x SLOWER,
// because a contended global lock is strictly worse than the per-thread
// caching the Windows heap already does, and this workload runs the multiply
// recursion across ~19 threads at once. A thread_local free list has no lock
// at all: a block is reused only by the thread that freed it, so there is
// nothing to contend on. Blocks freed on a different thread than they were
// allocated on simply join that thread's list, which is harmless.
//
// This only pays off because Parallel.cpp's workers are long-lived. With the
// previous std::async-per-task dispatch the pools would have been born cold
// and thrown away on every batch.
//
// Buckets are powers of two, so a block freed by one call is usable by any
// later call asking for up to that many bytes -- the recursive Toom splits
// and interpolation temporaries produce many distinct sizes, not a couple of
// fixed ones. Requests above POOL_MAX_BYTES bypass the pool entirely so the
// retained memory per thread stays bounded.

#pragma once

#include "common.h"

#include <cstddef>

namespace poolimpl {

inline constexpr size_t POOL_MAX_BYTES = 8u * 1024 * 1024;   // 8 MB cap per block

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
