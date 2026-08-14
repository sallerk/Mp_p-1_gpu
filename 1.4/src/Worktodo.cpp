// Copyright (C) Mp_p-1_gpu

#include "Worktodo.h"
#include "Config.h"

#include <cctype>
#include <cstdio>
#include <filesystem>
#include <system_error>

using namespace std;

namespace {

// A private copy of Config.cpp's own trim() -- not worth sharing a header
// over five lines, and this file has no other reason to depend on Config.cpp.
string trim(const string& s) {
  size_t a = 0, b = s.size();
  while (a < b && isspace((unsigned char) s[a])) { ++a; }
  while (b > a && isspace((unsigned char) s[b - 1])) { --b; }
  return s.substr(a, b - a);
}

// Deliberately NOT File's range-for line iterator (used by
// TuneEntry::readTuneFile and the vendored, currently-unused
// deleteLine()/copyWithout() in fs.cpp): File::readLine() throws if a line
// does not end in '\n', which is fine for a machine-written file like
// tune.txt but wrong here -- worktodo.txt is hand-edited, and a text editor
// not adding a trailing newline to the last line is routine, not an error.
string stripEol(const char* buf) {
  string raw = buf;
  while (!raw.empty() && (raw.back() == '\n' || raw.back() == '\r')) { raw.pop_back(); }
  return raw;
}

} // namespace

bool loadWorktodo(const string& path, vector<WorktodoEntry>& out, string& err) {
  out.clear();
  FILE* f = fopen(path.c_str(), "r");
  if (!f) { return true; }   // no file yet == empty queue, not an error

  char buf[1024];
  u32 lineNo = 0;
  while (fgets(buf, sizeof(buf), f)) {
    ++lineNo;
    const string t = trim(stripEol(buf));
    if (t.empty() || t[0] == '#') { continue; }

    u64 n = 0;
    if (!parseNumber(t, n) || n < 3) {
      fclose(f);
      err = "worktodo.txt line " + to_string(lineNo) + ": '" + t
          + "' is not a valid exponent (>= 3)";
      return false;
    }
    out.push_back({u32(n), lineNo});
  }
  fclose(f);
  return true;
}

bool consumeWorktodoEntry(const string& path, const WorktodoEntry& entry, string& err) {
  FILE* fi = fopen(path.c_str(), "r");
  if (!fi) { err = "worktodo.txt disappeared"; return false; }

  const string tmp = path + ".tmp";
  FILE* fo = fopen(tmp.c_str(), "w");
  if (!fo) { fclose(fi); err = "cannot open " + tmp; return false; }

  char buf[1024];
  u32 lineNo = 0;
  bool removed = false;
  bool ok = true;
  while (fgets(buf, sizeof(buf), fi)) {
    ++lineNo;
    if (lineNo == entry.lineNo) {
      // Confirm this is still the entry we loaded before dropping it -- the
      // file may have been hand-edited since. On a mismatch, fall through
      // and copy the line through unchanged instead of removing it.
      u64 n = 0;
      if (parseNumber(trim(stripEol(buf)), n) && u32(n) == entry.exponent) {
        removed = true;
        continue;
      }
    }
    if (fputs(buf, fo) == EOF) { ok = false; }
  }
  fclose(fi);
  fclose(fo);

  if (!ok || !removed) {
    remove(tmp.c_str());
    err = !ok ? "write failed"
              : "worktodo.txt changed since it was read -- exponent "
                + to_string(entry.exponent) + " was not removed";
    return false;
  }

  std::error_code ec;
  filesystem::rename(tmp, path, ec);
  if (ec) {
    // Windows will not rename onto an existing file.
    filesystem::remove(path, ec);
    filesystem::rename(tmp, path, ec);
  }
  if (ec) { err = "rename failed: " + ec.message(); return false; }
  return true;
}

namespace {

// "cannot collide with a real worktodo.txt" -- same reasoning as the
// checkpoint round-trip self-test's choice of exponent 1.
const char* TEST_FILE = "selftest-worktodo.txt";

bool writeRaw(const string& path, const string& content) {
  FILE* f = fopen(path.c_str(), "wb");
  if (!f) { return false; }
  const bool ok = fwrite(content.data(), 1, content.size(), f) == content.size();
  fclose(f);
  return ok;
}

} // namespace

int runWorktodoTests() {
  printf("worktodo self-test\n");
  int fails = 0;
  filesystem::remove(TEST_FILE);
  filesystem::remove(string(TEST_FILE) + ".tmp");

  printf("\n  A. round-trip parse: entries, comments, blank lines, order\n");
  {
    writeRaw(TEST_FILE,
             "# leading comment\n"
             "786613\n"
             "\n"
             "   # indented comment\n"
             "7873417\n"
             "28733813\n");
    vector<WorktodoEntry> out;
    string err;
    const bool ok = loadWorktodo(TEST_FILE, out, err);
    const bool shapeOk = ok && out.size() == 3
        && out[0].exponent == 786613  && out[0].lineNo == 2
        && out[1].exponent == 7873417 && out[1].lineNo == 5
        && out[2].exponent == 28733813 && out[2].lineNo == 6;
    if (!shapeOk) { ++fails; printf("  FAIL round-trip: got ok=%d size=%zu\n", ok, out.size()); }
    printf("     %s  %zu entries, comments/blanks skipped, line numbers preserved\n",
           shapeOk ? "PASS" : "FAIL", out.size());
  }

  printf("\n  B. consumeWorktodoEntry removes exactly one line\n");
  {
    vector<WorktodoEntry> out;
    string err;
    loadWorktodo(TEST_FILE, out, err);
    const bool consumed = consumeWorktodoEntry(TEST_FILE, out[1], err);   // the 7873417 entry

    vector<WorktodoEntry> after;
    loadWorktodo(TEST_FILE, after, err);
    const bool ok = consumed && after.size() == 2
        && after[0].exponent == 786613 && after[1].exponent == 28733813;
    if (!ok) { ++fails; printf("  FAIL consume: %s\n", err.c_str()); }
    printf("     %s  removed the targeted entry, left the other two in order\n", ok ? "PASS" : "FAIL");
  }

  printf("\n  C. malformed line is a hard error, not a silent skip\n");
  {
    writeRaw(TEST_FILE, "786613\nnot-a-number\n28733813\n");
    vector<WorktodoEntry> out;
    string err;
    const bool ok = loadWorktodo(TEST_FILE, out, err);
    const bool rejected = !ok && err.find("line 2") != string::npos;
    if (!rejected) { ++fails; printf("  FAIL malformed: ok=%d err='%s'\n", ok, err.c_str()); }
    printf("     %s  rejected with '%s'\n", rejected ? "PASS" : "FAIL", err.c_str());
  }

  printf("\n  D. missing file is an empty queue, not an error\n");
  {
    filesystem::remove(TEST_FILE);
    vector<WorktodoEntry> out;
    string err;
    const bool ok = loadWorktodo(TEST_FILE, out, err);
    const bool empty = ok && out.empty();
    if (!empty) { ++fails; printf("  FAIL missing file: ok=%d size=%zu\n", ok, out.size()); }
    printf("     %s  %s\n", empty ? "PASS" : "FAIL", ok ? "empty queue" : err.c_str());
  }

  printf("\n  E. last line without a trailing newline still parses and consumes\n");
  {
    // Deliberately no '\n' after the last entry -- the routine case a text
    // editor produces, and the reason File's throwing line iterator was not
    // reused here (see the comment on stripEol above).
    writeRaw(TEST_FILE, "786613\n28733813");
    vector<WorktodoEntry> out;
    string err;
    const bool loaded = loadWorktodo(TEST_FILE, out, err)
                      && out.size() == 2 && out[1].exponent == 28733813;
    const bool consumed = loaded && consumeWorktodoEntry(TEST_FILE, out[1], err);
    vector<WorktodoEntry> after;
    const bool afterOk = consumed && loadWorktodo(TEST_FILE, after, err)
                       && after.size() == 1 && after[0].exponent == 786613;
    if (!afterOk) { ++fails; printf("  FAIL no-trailing-newline: %s\n", err.c_str()); }
    printf("     %s  parsed and consumed without a trailing newline\n", afterOk ? "PASS" : "FAIL");
  }

  filesystem::remove(TEST_FILE);
  filesystem::remove(string(TEST_FILE) + ".tmp");

  printf("\n%s\n", fails ? "worktodo: FAILED" : "worktodo: all tests passed");
  return fails ? 1 : 0;
}
