// Copyright (C) Mp_p-1_gpu
//
// See Pool.h. Each bucket is a thread_local intrusive free list: a freed
// block's own first bytes hold the "next" pointer, so no separate node
// allocation is needed, and there is no lock because no other thread can
// reach the list. Blocks in bucket i are always exactly 2^i bytes, so any
// block popped from a bucket is large enough for the request that mapped to
// it.

#include "Pool.h"

#include <new>

namespace poolimpl {

namespace {

// Buckets for 2^0 .. 2^23 bytes (1 B .. 8 MB); POOL_MAX_BYTES enforces the top.
constexpr int NUM_BUCKETS = 24;

// Per-thread free lists. A block is only ever reused by the thread that freed
// it, which is what removes the locking entirely.
thread_local void* tFree[NUM_BUCKETS] = {};

int bucketIndex(size_t bytes) {
  size_t cap = 1;
  int idx = 0;
  while (cap < bytes) { cap <<= 1; ++idx; }
  return idx;
}

} // namespace

void* poolAllocate(size_t bytes) {
  if (bytes == 0) { bytes = 1; }
  if (bytes > POOL_MAX_BYTES) { return ::operator new(bytes); }

  const int idx = bucketIndex(bytes);
  void* head = tFree[idx];
  if (head) {
    tFree[idx] = *reinterpret_cast<void**>(head);
    return head;
  }
  return ::operator new(size_t(1) << idx);
}

void poolDeallocate(void* p, size_t bytes) noexcept {
  if (!p) { return; }
  if (bytes == 0) { bytes = 1; }
  if (bytes > POOL_MAX_BYTES) { ::operator delete(p); return; }

  const int idx = bucketIndex(bytes);
  *reinterpret_cast<void**>(p) = tFree[idx];
  tFree[idx] = p;
}

} // namespace poolimpl
