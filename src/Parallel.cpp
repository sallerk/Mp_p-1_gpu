// Copyright (C) Mp_p-1_gpu
//
// See Parallel.h.
//
// The workers are long-lived on purpose. This started as std::async per task,
// which is simpler, but sampling-profiling the real gcd workload showed ~45%
// of all CPU time inside ntdll's heap (every Nat operation allocates), so
// BigInt.cpp's Nat allocator became a thread-local free list -- and a
// thread_local pool is only worth anything if the threads themselves outlive
// a single task. std::async gives no such guarantee (MSVC spawns a fresh OS
// thread per call), so the pool below keeps a fixed set of workers alive for
// the process lifetime and hands them work through a queue.

#include "Parallel.h"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

namespace {

unsigned gThreadLimit = 0;   // 0 == all hardware threads
std::atomic<int> gActiveTasks{0};

int maxTasks() {
  if (gThreadLimit) { return int(gThreadLimit); }
  static const int n = std::max(1u, std::thread::hardware_concurrency());
  return n;
}

// A task plus the counter the submitter waits on.
struct Job {
  std::function<void()>* fn;
  std::atomic<int>* remaining;
};

class Pool {
public:
  // Deliberately leaked, not a plain function-local static: a `static Pool p`
  // IS destroyed at normal process exit (by the CRT's atexit machinery --
  // "lazily constructed" is not "never destroyed"), and workers_'s std::thread
  // destructors would then run on threads that are still joinable (they sit
  // forever in workerLoop()'s condvar wait, by design -- see ~Pool() below),
  // which std::terminate()s the process. Allocating with `new` and never
  // deleting is what actually achieves "never destroyed": the OS reclaims the
  // parked threads when the process exits, no destructor ever runs, nothing
  // to terminate() on.
  static Pool& instance() {
    static Pool* p = new Pool();
    return *p;
  }

  void submit(Job j) {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      queue_.push(j);
    }
    cv_.notify_one();
  }

  // Run queued work on the calling thread while waiting, so the submitter is
  // never merely blocked -- same "the caller does its share" property the
  // std::async version had.
  void helpUntilZero(std::atomic<int>& remaining) {
    while (remaining.load(std::memory_order_acquire) > 0) {
      Job j{};
      {
        std::lock_guard<std::mutex> lock(mutex_);
        if (queue_.empty()) { continue; }
        j = queue_.front();
        queue_.pop();
      }
      runJob(j);
    }
  }

private:
  Pool() {
    const int n = std::max(1, maxTasks());
    workers_.reserve(size_t(n));
    for (int i = 0; i < n; ++i) { workers_.emplace_back([this] { workerLoop(); }); }
  }

  // Never runs -- see instance()'s comment: the Pool itself is deliberately
  // leaked, specifically so this destructor is never reached. Joining at
  // process exit would buy nothing while risking a shutdown deadlock against
  // a worker mid-task, and destroying the still-joinable worker threads
  // without joining them is worse than not destroying them at all (that is
  // the bug this leak avoids).
  ~Pool() = default;

  static void runJob(Job& j) {
    (*j.fn)();
    j.remaining->fetch_sub(1, std::memory_order_release);
  }

  void workerLoop() {
    for (;;) {
      Job j{};
      {
        std::unique_lock<std::mutex> lock(mutex_);
        cv_.wait(lock, [this] { return !queue_.empty(); });
        j = queue_.front();
        queue_.pop();
      }
      runJob(j);
    }
  }

  std::mutex mutex_;
  std::condition_variable cv_;
  std::queue<Job> queue_;
  std::vector<std::thread> workers_;
};

} // namespace

void setParallelThreads(unsigned n) { gThreadLimit = n; }

void runParallel(std::vector<std::function<void()>> tasks, bool worthIt) {
  const int budget = maxTasks() - gActiveTasks.load() - 1;
  if (!worthIt || tasks.size() < 2 || budget < 1) {
    for (auto& t : tasks) { t(); }
    return;
  }
  const size_t offload = std::min(tasks.size() - 1, size_t(budget));
  gActiveTasks += int(offload);

  std::atomic<int> remaining{int(offload)};
  Pool& pool = Pool::instance();
  for (size_t i = 0; i < offload; ++i) { pool.submit(Job{&tasks[i], &remaining}); }
  // The caller's own thread does the remainder rather than idling.
  for (size_t i = offload; i < tasks.size(); ++i) { tasks[i](); }
  pool.helpUntilZero(remaining);

  gActiveTasks -= int(offload);
}
