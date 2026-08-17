// Copyright (C) Mp_p-1_gpu
//
// See Parallel.h.

#include "Parallel.h"

#include <algorithm>
#include <atomic>
#include <future>
#include <thread>

namespace {
unsigned gThreadLimit = 0;   // 0 == all hardware threads
std::atomic<int> gActiveTasks{0};

int maxTasks() {
  if (gThreadLimit) { return int(gThreadLimit); }
  static const int n = std::max(1u, std::thread::hardware_concurrency());
  return n;
}
} // namespace

void setParallelThreads(unsigned n) { gThreadLimit = n; }

void runParallel(std::vector<std::function<void()>> tasks, bool worthIt) {
  const int budget = maxTasks() - gActiveTasks.load() - 1;
  if (!worthIt || tasks.size() < 2 || budget < 1) {
    for (auto& t : tasks) { t(); }
    return;
  }
  const size_t spawn = std::min(tasks.size() - 1, size_t(budget));
  gActiveTasks += int(spawn);
  std::vector<std::future<void>> futures;
  futures.reserve(spawn);
  for (size_t i = 0; i < spawn; ++i) {
    futures.push_back(std::async(std::launch::async, tasks[i]));
  }
  // The caller's own thread does the remainder rather than idling.
  for (size_t i = spawn; i < tasks.size(); ++i) { tasks[i](); }
  for (auto& f : futures) { f.get(); }
  gActiveTasks -= int(spawn);
}
