// Copyright (C) Mp_p-1_gpu
//
// See Pool.h. Each bucket is a Treiber stack (lock-free singly linked list,
// pushed/popped via CAS on an atomic head pointer): a freed block's own first
// bytes are reinterpreted as the "next" pointer, so no separate node
// allocation is needed. Blocks in bucket i are always exactly 2^i bytes, so
// any block popped from a bucket is guaranteed large enough for the request
// that mapped to it.

#include "Pool.h"

#include <new>

namespace poolimpl {

namespace {

// Buckets for 2^0 .. 2^25 bytes (1 .. 32 MB); POOL_MAX_BYTES enforces the top.
constexpr int NUM_BUCKETS = 26;

struct Bucket {
  std::atomic<void*> head{nullptr};
};

Bucket gBuckets[NUM_BUCKETS];

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
  Bucket& bk = gBuckets[idx];
  void* head = bk.head.load(std::memory_order_acquire);
  while (head) {
    void* next = *reinterpret_cast<void* const*>(head);
    if (bk.head.compare_exchange_weak(head, next, std::memory_order_acq_rel,
                                       std::memory_order_acquire)) {
      return head;
    }
    // head was refreshed by the failed compare_exchange_weak; retry.
  }
  return ::operator new(size_t(1) << idx);
}

void poolDeallocate(void* p, size_t bytes) noexcept {
  if (!p) { return; }
  if (bytes == 0) { bytes = 1; }
  if (bytes > POOL_MAX_BYTES) { ::operator delete(p); return; }

  const int idx = bucketIndex(bytes);
  Bucket& bk = gBuckets[idx];
  void* head = bk.head.load(std::memory_order_relaxed);
  do {
    *reinterpret_cast<void**>(p) = head;
  } while (!bk.head.compare_exchange_weak(head, p, std::memory_order_release,
                                           std::memory_order_relaxed));
}

} // namespace poolimpl
