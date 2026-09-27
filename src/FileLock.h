// Copyright (C) Mp_p-1_gpu
//
// The lock-file convention GIMPS programs share for worktodo.txt and
// results.txt (mfaktc, mfakto, CUDALucas, AutoPrimeNet): before touching
// FILE, create FILE.lck exclusively; while it exists, someone else is editing
// FILE; delete it when done. There is no stale-lock detection anywhere, so a
// holder must not die holding one.
//
// This program holds its lock as an open handle created with delete-on-close,
// so however the process ends -- Ctrl-C, closing the console window, a kill,
// a crash -- Windows removes the file. To everyone else it is an ordinary
// FILE.lck: their exclusive create fails because the name exists.
//
// Waiting follows the convention -- forever, while the file exists -- with a
// message after a couple of seconds and a stale-lock hint after a minute.
// Ctrl-C ends the wait after a short grace period, one of two ways (Options).

#pragma once

#include <atomic>
#include <functional>
#include <string>

class FileLock {
public:
  // What Ctrl-C while waiting means. GiveUp: leave FILE alone (state
  // Interrupted). Proceed: go ahead without the lock (state Unlocked) -- for
  // writing a finished result, which must never be dropped.
  enum class OnInterrupt { GiveUp, Proceed };

  // Unlocked: not holding it -- Proceed after Ctrl-C, a lock file that cannot
  // be created at all (a read-only folder, say) where waiting would never
  // end, or after release(). TimedOut exists for the self-tests only.
  enum class State { Held, Unlocked, Interrupted, TimedOut };

  enum class Msg { Waiting, StaleHint, Acquired, Unlocked };

  struct Options {
    OnInterrupt onInterrupt = OnInterrupt::GiveUp;
    // Self-test hooks. The defaults are gInterrupted, log(), and no limit.
    const std::atomic<bool>* stop = nullptr;
    std::function<void(Msg, const std::string&)> say;
    double maxWaitSec = -1;
    double quietSec = 2.0;       // silent this long: holders normally take milliseconds
    double staleHintSec = 60.0;
    double graceSec = 1.0;       // keep trying this long after Ctrl-C
    double transientSec = 5.0;   // access denied this long means "cannot lock here"
  };

  explicit FileLock(const std::string& file) : FileLock(file, Options{}) {}
  FileLock(const std::string& file, Options opt);
  ~FileLock() { release(); }
  FileLock(const FileLock&) = delete;
  FileLock& operator=(const FileLock&) = delete;

  State state() const { return state_; }
  bool held() const { return state_ == State::Held; }
  const std::string& lockPath() const { return lockPath_; }

  // Removes the lock file. Idempotent; the destructor calls it.
  void release();

  // Self-tests only. While one exists, every FileLock takes its stop flag,
  // messages, limit and thresholds from `o` (OnInterrupt still comes from the
  // call site), so code that locks internally -- the worktodo and results
  // functions -- can be tested without a real Ctrl-C or an endless wait.
  class TestDefaults {
  public:
    explicit TestDefaults(const Options& o);
    ~TestDefaults();
    TestDefaults(const TestDefaults&) = delete;
    TestDefaults& operator=(const TestDefaults&) = delete;
  };

private:
  std::string lockPath_;
  void* handle_ = nullptr;       // a HANDLE; windows.h stays out of this header
  State state_ = State::Unlocked;
};
