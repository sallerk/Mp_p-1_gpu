// Copyright (C) Mp_p-1_gpu

#include "FileLock.h"

#include "log.h"
#include "timeutil.h"

// common.h (via timeutil.h) says `using namespace std`, so the full
// windows.h would see std::byte clash with its own byte.
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#include <algorithm>
#include <cstdio>
#include <mutex>
#include <set>

using namespace std;

extern std::atomic<bool> gInterrupted;

namespace {

std::mutex gMutex;
// Lock files this process holds. A second FileLock on one of them would wait
// for itself forever.
set<string> gHeld;
// Lock files that could not be created at all. Checked once more per use,
// but without the transientSec retry, which would otherwise be paid on every
// wait_for_work poll.
set<string> gCannotLock;

const FileLock::Options* gTestDefaults = nullptr;

string winError(DWORD e) {
  char buf[256] = "";
  FormatMessageA(FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS, nullptr, e, 0,
                 buf, sizeof(buf), nullptr);
  string s = buf;
  while (!s.empty() && (s.back() == '\n' || s.back() == '\r' || s.back() == ' ' || s.back() == '.')) {
    s.pop_back();
  }
  return s.empty() ? "error " + to_string(e) : s;
}

string secs(double s) {
  char buf[32];
  snprintf(buf, sizeof(buf), "%.1f s", s);
  return buf;
}

} // namespace

FileLock::TestDefaults::TestDefaults(const Options& o) { gTestDefaults = &o; }
FileLock::TestDefaults::~TestDefaults() { gTestDefaults = nullptr; }

FileLock::FileLock(const string& file, Options opt) : lockPath_(file + ".lck") {
  if (gTestDefaults) {
    const OnInterrupt mode = opt.onInterrupt;
    opt = *gTestDefaults;
    opt.onInterrupt = mode;
  }
  std::function<void(Msg, const string&)> say = opt.say;
  if (!say) { say = [](Msg, const string& s) { log("%s\n", s.c_str()); }; }
  const std::atomic<bool>& stop = opt.stop ? *opt.stop : gInterrupted;

  bool knownCannotLock;
  {
    lock_guard<mutex> g(gMutex);
    if (gHeld.count(lockPath_)) {
      say(Msg::Unlocked, "  BUG: " + lockPath_ + " is already held by this program -- carrying on without it");
      state_ = State::Unlocked;
      return;
    }
    knownCannotLock = gCannotLock.count(lockPath_) != 0;
  }

  auto cannotLock = [&](DWORD e) {
    lock_guard<mutex> g(gMutex);
    if (gCannotLock.insert(lockPath_).second) {
      say(Msg::Unlocked, "  WARNING: cannot create " + lockPath_ + " (" + winError(e)
                         + ") -- using " + file + " without the lock");
    }
    state_ = State::Unlocked;
  };

  Timer waited;
  double stopSeenAt = -1, deniedSince = -1;
  bool saidWaiting = false, saidStale = false;
  u32 sleepMs = 0;
  for (;;) {
    // No FILE_SHARE_DELETE: a `del` of a live lock fails instead of quietly
    // breaking the exclusion. DELETE access is for the POSIX delete below.
    HANDLE h = CreateFileA(lockPath_.c_str(), GENERIC_WRITE | DELETE,
                           FILE_SHARE_READ | FILE_SHARE_WRITE, nullptr, CREATE_NEW,
                           FILE_ATTRIBUTE_NORMAL | FILE_FLAG_DELETE_ON_CLOSE, nullptr);
    if (h != INVALID_HANDLE_VALUE) {
      handle_ = h;
      state_ = State::Held;
      {
        lock_guard<mutex> g(gMutex);
        gHeld.insert(lockPath_);
        gCannotLock.erase(lockPath_);
      }
      if (saidWaiting) { say(Msg::Acquired, "  got " + lockPath_ + " after " + secs(waited.at())); }
      return;
    }
    const DWORD e = GetLastError();
    const double t = waited.at();

    if (e == ERROR_FILE_EXISTS || e == ERROR_ALREADY_EXISTS) {
      deniedSince = -1;                                  // held by someone: wait
    } else if (!knownCannotLock &&
               (e == ERROR_ACCESS_DENIED || e == ERROR_SHARING_VIOLATION || e == ERROR_LOCK_VIOLATION)) {
      // Also what a lock file that is being deleted returns, so retry -- but
      // a read-only folder says the same forever.
      if (deniedSince < 0) { deniedSince = t; }
      if (t - deniedSince >= opt.transientSec) { cannotLock(e); return; }
    } else {
      cannotLock(e);
      return;
    }

    if (stop.load()) {
      if (stopSeenAt < 0) { stopSeenAt = t; }
      if (t - stopSeenAt >= opt.graceSec) {
        if (opt.onInterrupt == OnInterrupt::GiveUp) {
          state_ = State::Interrupted;
        } else {
          say(Msg::Unlocked, "  interrupted while waiting for " + lockPath_ + " -- writing "
                             + file + " without it");
          state_ = State::Unlocked;
        }
        return;
      }
    }
    if (opt.maxWaitSec >= 0 && t >= opt.maxWaitSec) {
      state_ = State::TimedOut;
      return;
    }
    if (!saidWaiting && t >= opt.quietSec) {
      say(Msg::Waiting, "  waiting for " + lockPath_ + " -- another program (AutoPrimeNet?) is"
                        " using " + file);
      saidWaiting = true;
    }
    if (!saidStale && t >= opt.staleHintSec) {
      say(Msg::StaleHint, "  still waiting for " + lockPath_ + ". If no other program is"
                          " running here, the lock is stale: delete " + lockPath_);
      saidStale = true;
    }

    // mfaktc's backoff, 1 ms longer per try up to 1 s, slept in slices so a
    // Ctrl-C is noticed promptly. After one, retry often through the grace.
    sleepMs = stopSeenAt >= 0 ? 20 : std::min<u32>(sleepMs + 1, 1000);
    for (u32 slept = 0; slept < sleepMs; ) {
      const u32 slice = std::min<u32>(100, sleepMs - slept);
      Timer::usleep(slice * 1000);
      slept += slice;
      if (stopSeenAt < 0 && stop.load()) { break; }
    }
  }
}

void FileLock::release() {
  if (!handle_) { return; }
  // Delete with POSIX semantics first: the name goes at once even if another
  // process -- an antivirus scanner, say -- has the file open. Delete-on-close
  // alone would leave it "delete pending" until they close it, and an
  // exclusive create meanwhile fails with access denied, which AutoPrimeNet
  // treats as an error rather than as "busy". Best effort: delete-on-close
  // still removes it when the handle closes.
  struct { DWORD flags; } disposition = {0x1 | 0x2};   // FILE_DISPOSITION_FLAG_DELETE | _POSIX_SEMANTICS
  SetFileInformationByHandle(handle_, FILE_INFO_BY_HANDLE_CLASS(21),   // FileDispositionInfoEx
                             &disposition, sizeof(disposition));
  CloseHandle(handle_);
  handle_ = nullptr;
  state_ = State::Unlocked;
  lock_guard<mutex> g(gMutex);
  gHeld.erase(lockPath_);
}
