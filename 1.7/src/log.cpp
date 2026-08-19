// Copyright (C) Mihai Preda

#include "log.h"
#include "File.h"
#include "timeutil.h"

#include <mutex>

thread_local string context;
thread_local vector<string> contextParts;

thread_local File logFile;

File stdoutFile{stdout, "stdout"};

string logContext() { return context; }

void initLog(const char *logName) {
  assert(!logFile);
  logFile = File::openAppend(logName);
}

string longTimeStr()  { return timeStr("%Y-%m-%d %H:%M:%S %Z"); }
string shortTimeStr() { return timeStr("%Y%m%d %H:%M:%S"); }

static char logBuf[32 * 1024];

void log(const char *fmt, ...) {
  static std::mutex logMutex;

  string prefix = shortTimeStr() + ' ' + context;

  std::unique_lock lock(logMutex);

  // Upstream appended "%n" to the prefix format to learn its length. MSVC's CRT
  // disables %n by default: it invokes the invalid-parameter handler, which
  // terminates the process (0xC0000409) before anything is printed. snprintf's
  // return value is that same count, and is portable.
  int pos = snprintf(logBuf, sizeof(logBuf), "%s ", prefix.c_str());
  if (pos < 0) { pos = 0; }
  if (pos > int(sizeof(logBuf)) - 1) { pos = int(sizeof(logBuf)) - 1; }  // prefix truncated

  va_list va;
  va_start(va, fmt);
  vsnprintf(logBuf + pos, sizeof(logBuf) - pos, fmt, va);
  va_end(va);
  string_view s{logBuf};

  // stdout is already unbuffered at the CRT level (main.cpp's own
  // setvbuf(..., _IONBF, ...), for the identical reason), but a file opened
  // via File::openAppend is fully buffered by default -- File::write() is a
  // raw fwrite() with no flush, built for high-throughput data (checkpoints)
  // where flushing every call would be a real cost. Log lines are the
  // opposite: low-frequency, and the entire point of writing them to a file
  // is to survive a crash or a forced kill, which a full stdio buffer does
  // not -- it only reaches disk on a clean exit. Flush explicitly here,
  // where the cost is negligible.
  if (logFile) { logFile.write(s); logFile.flush(); }
  stdoutFile.write(s);
}

LogContext::LogContext(const string& s) : part{s} {
  contextParts.push_back(s);
  context += s;
}

LogContext::~LogContext() {
  assert(!contextParts.empty());
  assert(context.size() >= contextParts.back().size());
  context = context.substr(0, context.size() - contextParts.back().size());
  contextParts.pop_back();
}
