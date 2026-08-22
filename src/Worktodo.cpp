// Copyright (C) Mp_p-1_gpu

#include "Worktodo.h"
#include "Config.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
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

string lower(string s) {
  for (char& c : s) { c = (char) tolower((unsigned char) c); }
  return s;
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

// ---- Pfactor=/Pminus1= assignment-line parsing -----------------------------
//
// Ground truth verified against Prime95's own C parser (primesearch/Prime95,
// commonc.c, parseWorkToDoLine) rather than assumed -- both keywords, and the
// optional assignment-ID prefix, follow that file's exact field order and
// positional-optional-field quirks. AutoPrimeNet writes this same shape.

// A field is one comma-separated token, EXCEPT that a comma inside a "..."
// quoted span (the known-factors list) does not split -- that quoted field
// is one token with its quotes still attached, unwrapped separately below.
vector<string> splitTopLevel(const string& s) {
  vector<string> out;
  string cur;
  bool inQuotes = false;
  for (char c : s) {
    if (c == '"') { inQuotes = !inQuotes; cur += c; }
    else if (c == ',' && !inQuotes) { out.push_back(cur); cur.clear(); }
    else { cur += c; }
  }
  out.push_back(cur);
  return out;
}

bool parseDoubleField(const string& s, double& out) {
  if (s.empty()) { return false; }
  char* end = nullptr;
  out = strtod(s.c_str(), &end);
  return end == s.c_str() + s.size();
}

bool isQuoted(const string& s) { return !s.empty() && s[0] == '"'; }

// "12345,67890" (quotes included, as splitTopLevel leaves them) -> the two
// factor strings, trimmed, validated as non-empty decimal integers.
bool parseKnownFactors(const string& field, vector<string>& out, string& what) {
  if (field.size() < 2 || field.back() != '"') { what = "unterminated known-factors quote"; return false; }
  const string inner = field.substr(1, field.size() - 2);
  size_t pos = 0;
  while (pos <= inner.size()) {
    const size_t comma = inner.find(',', pos);
    const string one = trim(inner.substr(pos, comma == string::npos ? string::npos : comma - pos));
    if (one.empty() || one.find_first_not_of("0123456789") != string::npos) {
      what = "known-factors entry '" + one + "' is not a decimal integer";
      return false;
    }
    out.push_back(one);
    if (comma == string::npos) { break; }
    pos = comma + 1;
  }
  return !out.empty();
}

// Strips a leading 32-hex-char assignment ID (immediately followed by ',')
// from `value` if present -- generic to ANY worktodo keyword in Prime95, not
// specific to Pfactor=/Pminus1=. Ported verbatim from commonc.c's own loop.
string stripAid(string& value) {
  size_t i = 0;
  for (; i < value.size(); ++i) {
    const char c = value[i];
    const bool hex = (c >= '0' && c <= '9')
                   || (toupper((unsigned char) c) >= 'A' && toupper((unsigned char) c) <= 'F');
    if (!hex) { break; }
    if (i == 31) {
      if (value.size() > 32 && value[32] == ',') {
        const string aid = value.substr(0, 32);
        value = value.substr(33);
        return aid;
      }
      break;
    }
  }
  return "";
}

// value is everything after "PFACTOR=" or "PMINUS1=", with any AID already
// stripped. Parses the shared k,b,n,c prefix both keywords require, validates
// it describes a Mersenne number (this program factors nothing else), and
// returns the remaining fields for the caller's keyword-specific tail.
bool parseKbnc(const vector<string>& f, u32 lineNo, u32& exponent, string& err) {
  if (f.size() < 4) {
    err = "worktodo.txt line " + to_string(lineNo) + ": expected k,b,n,c fields";
    return false;
  }
  // k, b, n, c all parsed via double then range/value-checked, same approach
  // Prime95's own atof-then-cast takes -- b/n are unsigned in practice, c is
  // signed (-1 for a Mersenne number), k is always exactly 1.0 here.
  double k = 0, b = 0, n = 0, c = 0;
  if (!parseDoubleField(f[0], k) || !parseDoubleField(f[1], b)
      || !parseDoubleField(f[2], n) || !parseDoubleField(f[3], c)) {
    err = "worktodo.txt line " + to_string(lineNo) + ": malformed k,b,n,c";
    return false;
  }
  if (k != 1.0 || u64(b + 0.5) != 2 || c != -1.0) {
    err = "worktodo.txt line " + to_string(lineNo) + ": k=" + f[0] + " b=" + f[1]
        + " c=" + f[3] + " is not a Mersenne number (this program only factors M_p = 2^p - 1)";
    return false;
  }
  if (n < 3) {
    err = "worktodo.txt line " + to_string(lineNo) + ": exponent " + f[2] + " is too small";
    return false;
  }
  exponent = u32(n + 0.5);
  return true;
}

bool parsePfactor(string value, u32 lineNo, WorktodoEntry& out, string& err) {
  out.aid = stripAid(value);
  const vector<string> f = splitTopLevel(value);
  if (!parseKbnc(f, lineNo, out.exponent, err)) { return false; }
  // New style only (k,b,n,c,how_far_factored,tests_saved[,"known_factors"]) --
  // the pre-k,b,n,c old style predates AutoPrimeNet and PrimeNet's own API by
  // years and is not something it emits.
  if (f.size() < 6) {
    err = "worktodo.txt line " + to_string(lineNo)
        + ": Pfactor= needs k,b,n,c,how_far_factored,tests_saved";
    return false;
  }
  double howFar = 0, testsSaved = 0;
  if (!parseDoubleField(f[4], howFar) || !parseDoubleField(f[5], testsSaved)) {
    err = "worktodo.txt line " + to_string(lineNo) + ": malformed how_far_factored/tests_saved";
    return false;
  }
  if (howFar < 1 || howFar > 127) {
    err = "worktodo.txt line " + to_string(lineNo)
        + ": how_far_factored " + f[4] + " must be 1..127 bits";
    return false;
  }
  out.hasFactoredTo = true;
  out.factoredTo = u32(floor(howFar));
  out.testsSaved = testsSaved;
  if (f.size() > 6) {
    if (!isQuoted(f[6])) {
      err = "worktodo.txt line " + to_string(lineNo) + ": expected \"known_factors\" after tests_saved";
      return false;
    }
    string what;
    if (!parseKnownFactors(f[6], out.knownFactors, what)) {
      err = "worktodo.txt line " + to_string(lineNo) + ": " + what;
      return false;
    }
  }
  return true;
}

bool parsePminus1(string value, u32 lineNo, WorktodoEntry& out, string& err) {
  out.aid = stripAid(value);
  const vector<string> f = splitTopLevel(value);
  if (!parseKbnc(f, lineNo, out.exponent, err)) { return false; }
  if (f.size() < 6) {
    err = "worktodo.txt line " + to_string(lineNo) + ": Pminus1= needs k,b,n,c,B1,B2";
    return false;
  }
  double b1 = 0, b2 = 0;
  if (!parseDoubleField(f[4], b1) || !parseDoubleField(f[5], b2) || b1 < 2 || b2 < 2) {
    err = "worktodo.txt line " + to_string(lineNo) + ": malformed B1/B2";
    return false;
  }
  out.assignedB1 = u64(b1 + 0.5);
  out.assignedB2 = u64(b2 + 0.5);
  if (out.assignedB2 <= out.assignedB1) {
    err = "worktodo.txt line " + to_string(lineNo) + ": B2 must be greater than B1";
    return false;
  }
  out.hasAssignedBounds = true;

  // The optional tail (sieve_depth, B2_start, known_factors) is POSITIONAL,
  // not tagged -- each of the first two slots is consumed whether or not its
  // value passes its own sanity check, which is what lets Prime95 (and this
  // parser) tell "this is a number" from "this is the known-factors quote"
  // apart. Ported verbatim, including sieve_depth's own >300-is-discarded
  // rule -- that is not a bug, it is how a real line disambiguates a TF-depth
  // field from a B2_start field without either being tagged.
  size_t cursor = 6;
  if (cursor < f.size() && !isQuoted(f[cursor])) {
    double sieve = 0;
    if (parseDoubleField(f[cursor], sieve) && sieve >= 1 && sieve <= 127) {
      out.hasFactoredTo = true;
      out.factoredTo = u32(floor(sieve));
    }
    ++cursor;
  }
  if (cursor < f.size() && !isQuoted(f[cursor])) {
    double start = 0;
    if (parseDoubleField(f[cursor], start) && start > double(out.assignedB1)) {
      out.b2Start = u64(start + 0.5);
    }
    ++cursor;
  }
  if (cursor < f.size()) {
    if (!isQuoted(f[cursor])) {
      err = "worktodo.txt line " + to_string(lineNo) + ": unexpected extra field '" + f[cursor] + "'";
      return false;
    }
    string what;
    if (!parseKnownFactors(f[cursor], out.knownFactors, what)) {
      err = "worktodo.txt line " + to_string(lineNo) + ": " + what;
      return false;
    }
  }
  return true;
}

// Parses one non-empty, non-comment worktodo.txt line: either a bare
// exponent (this program's own shorthand) or a Pfactor=/Pminus1= assignment
// line (AutoPrimeNet's/Prime95's own format). Shared by loadWorktodo and
// consumeWorktodoEntry's re-validation guard, so the two can never drift.
bool parseWorktodoLine(const string& t, u32 lineNo, WorktodoEntry& out, string& err) {
  out = WorktodoEntry{};
  out.lineNo = lineNo;

  const size_t eq = t.find('=');
  if (eq == string::npos) {
    u64 n = 0;
    if (!parseNumber(t, n) || n < 3) {
      err = "worktodo.txt line " + to_string(lineNo) + ": '" + t + "' is not a valid exponent (>= 3)";
      return false;
    }
    out.exponent = u32(n);
    return true;
  }

  const string key = lower(trim(t.substr(0, eq)));
  const string value = trim(t.substr(eq + 1));
  if (key == "pfactor") { return parsePfactor(value, lineNo, out, err); }
  if (key == "pminus1") { return parsePminus1(value, lineNo, out, err); }
  err = "worktodo.txt line " + to_string(lineNo) + ": keyword '" + key
      + "' is not supported -- this program only factors (use Pfactor= or Pminus1=)";
  return false;
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

    WorktodoEntry entry;
    if (!parseWorktodoLine(t, lineNo, entry, err)) {
      fclose(f);
      return false;
    }
    out.push_back(entry);
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
      // file may have been hand-edited, or reissued by AutoPrimeNet under a
      // new assignment, since. On a mismatch, fall through and copy the line
      // through unchanged instead of removing it.
      WorktodoEntry reread;
      string perr;
      const string t = trim(stripEol(buf));
      const bool sameEntry = !t.empty() && t[0] != '#'
          && parseWorktodoLine(t, lineNo, reread, perr)
          && reread.exponent == entry.exponent
          && (entry.aid.empty() || reread.aid == entry.aid);
      if (sameEntry) {
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

  printf("\n  F. Pfactor= round-trip (AID, known-factors, both, neither)\n");
  {
    struct Case { const char* line; bool wantAid; bool wantKF; };
    const Case cases[] = {
      {"Pfactor=1,2,86243,-1,74.5,1.30", false, false},
      {"Pfactor=AABBCCDD00112233445566778899AABB,1,2,86243,-1,74.5,1.30", true, false},
      {"Pfactor=1,2,86243,-1,74.5,1.30,\"170141183460469231731687303715884105727\"", false, true},
      {"Pfactor=AABBCCDD00112233445566778899AABB,1,2,86243,-1,74.5,1.30,\"123,456\"", true, true},
    };
    bool allOk = true;
    for (const Case& c : cases) {
      writeRaw(TEST_FILE, string(c.line) + "\n");
      vector<WorktodoEntry> out;
      string err;
      const bool ok = loadWorktodo(TEST_FILE, out, err) && out.size() == 1
          && out[0].exponent == 86243 && out[0].hasFactoredTo && out[0].factoredTo == 74
          && out[0].testsSaved > 1.29 && out[0].testsSaved < 1.31
          && (!c.wantAid || out[0].aid == "AABBCCDD00112233445566778899AABB")
          && (c.wantAid  || out[0].aid.empty())
          && (!c.wantKF  || !out[0].knownFactors.empty())
          && (c.wantKF   || out[0].knownFactors.empty());
      if (!ok) { allOk = false; printf("  FAIL Pfactor case '%s': %s\n", c.line, err.c_str()); }
    }
    if (!allOk) { ++fails; }
    printf("     %s  AID/known-factors present or absent, both parse correctly\n", allOk ? "PASS" : "FAIL");
  }

  printf("\n  G. Pminus1= round-trip (bare, sieve_depth, B2_start, both+known-factors, B2_start<=B1)\n");
  {
    bool allOk = true;

    writeRaw(TEST_FILE, "Pminus1=1,2,86243,-1,50000,3000000\n");
    { vector<WorktodoEntry> out; string err;
      const bool ok = loadWorktodo(TEST_FILE, out, err) && out.size() == 1
          && out[0].hasAssignedBounds && out[0].assignedB1 == 50000 && out[0].assignedB2 == 3000000
          && !out[0].hasFactoredTo && out[0].b2Start == 0;
      if (!ok) { allOk = false; printf("  FAIL Pminus1 bare: %s\n", err.c_str()); } }

    writeRaw(TEST_FILE, "Pminus1=1,2,86243,-1,50000,3000000,74.0\n");
    { vector<WorktodoEntry> out; string err;
      const bool ok = loadWorktodo(TEST_FILE, out, err) && out.size() == 1
          && out[0].hasFactoredTo && out[0].factoredTo == 74 && out[0].b2Start == 0;
      if (!ok) { allOk = false; printf("  FAIL Pminus1 sieve_depth only: %s\n", err.c_str()); } }

    // B2_start present with NO preceding sieve_depth -- the field-shift edge
    // case. Prime95's own positional walk means this value lands in the
    // sieve_depth SLOT first; it is > 300 so it is discarded there, then the
    // walk's second slot is never reached (only 7 fields total), so
    // b2Start stays 0. This is upstream's real, if surprising, behavior --
    // an AutoPrimeNet-written line never omits sieve_depth while giving
    // B2_start, so this shape should not occur in practice.
    writeRaw(TEST_FILE, "Pminus1=1,2,86243,-1,50000,3000000,3100000\n");
    { vector<WorktodoEntry> out; string err;
      const bool ok = loadWorktodo(TEST_FILE, out, err) && out.size() == 1
          && !out[0].hasFactoredTo && out[0].b2Start == 0;
      if (!ok) { allOk = false; printf("  FAIL Pminus1 B2_start-in-sieve-slot: %s\n", err.c_str()); } }

    writeRaw(TEST_FILE,
             "Pminus1=AABBCCDD00112233445566778899AABB,1,2,86243,-1,50000,3000000,74.0,3100000,\"12345\"\n");
    { vector<WorktodoEntry> out; string err;
      const bool ok = loadWorktodo(TEST_FILE, out, err) && out.size() == 1
          && out[0].aid == "AABBCCDD00112233445566778899AABB"
          && out[0].hasFactoredTo && out[0].factoredTo == 74
          && out[0].b2Start == 3100000 && out[0].knownFactors.size() == 1
          && out[0].knownFactors[0] == "12345";
      if (!ok) { allOk = false; printf("  FAIL Pminus1 full: %s\n", err.c_str()); } }

    writeRaw(TEST_FILE, "Pminus1=1,2,86243,-1,50000,3000000,74.0,40000\n");
    { vector<WorktodoEntry> out; string err;
      const bool ok = loadWorktodo(TEST_FILE, out, err) && out.size() == 1 && out[0].b2Start == 0;
      if (!ok) { allOk = false; printf("  FAIL Pminus1 B2_start<=B1 discarded: %s\n", err.c_str()); } }

    if (!allOk) { ++fails; }
    printf("     %s  every optional-field shape matches Prime95's own positional walk\n", allOk ? "PASS" : "FAIL");
  }

  printf("\n  H. AID-strip boundary: 32 hex chars strip, 31 do not\n");
  {
    writeRaw(TEST_FILE, "Pminus1=AABBCCDD00112233445566778899AABB,1,2,86243,-1,50000,3000000\n");
    vector<WorktodoEntry> out32; string err;
    const bool ok32 = loadWorktodo(TEST_FILE, out32, err) && out32.size() == 1
        && out32[0].aid == "AABBCCDD00112233445566778899AABB" && out32[0].exponent == 86243;

    // 31 hex chars: the AID-strip loop must NOT fire (needs exactly 32 then
    // a comma) -- this line is then a malformed k,b,n,c (k would swallow
    // the extra digit) and correctly rejected rather than silently
    // misparsed with a truncated AID.
    writeRaw(TEST_FILE, "Pminus1=AABBCCDD00112233445566778899AAB,1,2,86243,-1,50000,3000000\n");
    vector<WorktodoEntry> out31;
    const bool rejected31 = !loadWorktodo(TEST_FILE, out31, err);

    const bool ok = ok32 && rejected31;
    if (!ok) { ++fails; printf("  FAIL AID boundary: ok32=%d rejected31=%d err='%s'\n", ok32, rejected31, err.c_str()); }
    printf("     %s  32 hex chars + comma strips as AID, 31 does not\n", ok ? "PASS" : "FAIL");
  }

  printf("\n  I. rejected shapes: wrong k/b/c, unterminated known-factors quote\n");
  {
    bool allOk = true;
    const char* bad[] = {
      "Pminus1=2,2,86243,-1,50000,3000000",     // k != 1
      "Pminus1=1,3,86243,-1,50000,3000000",     // b != 2
      "Pminus1=1,2,86243,1,50000,3000000",      // c != -1
      "Pfactor=1,2,86243,-1,74.0,1.3,\"12345",  // unterminated quote
    };
    for (const char* line : bad) {
      writeRaw(TEST_FILE, string(line) + "\n");
      vector<WorktodoEntry> out; string err;
      if (loadWorktodo(TEST_FILE, out, err)) { allOk = false; printf("  FAIL accepted bad line '%s'\n", line); }
    }
    if (!allOk) { ++fails; }
    printf("     %s  non-Mersenne k/b/c and an unterminated quote are hard errors\n", allOk ? "PASS" : "FAIL");
  }

  printf("\n  J. unsupported worktodo keyword is a hard error naming it\n");
  {
    writeRaw(TEST_FILE, "PRP=1,2,86243,-1,74,1\n");
    vector<WorktodoEntry> out; string err;
    const bool ok = !loadWorktodo(TEST_FILE, out, err) && err.find("prp") != string::npos;
    if (!ok) { ++fails; printf("  FAIL unsupported keyword: err='%s'\n", err.c_str()); }
    printf("     %s  rejected with '%s'\n", ok ? "PASS" : "FAIL", err.c_str());
  }

  printf("\n  K. mixed file: bare + Pfactor= + Pminus1=, consume removes the right one\n");
  {
    writeRaw(TEST_FILE,
             "786613\n"
             "Pfactor=1,2,7873417,-1,72.0,1.0\n"
             "Pminus1=1,2,28733813,-1,50000,3000000\n");
    vector<WorktodoEntry> out; string err;
    const bool loaded = loadWorktodo(TEST_FILE, out, err) && out.size() == 3
        && out[0].exponent == 786613   && !out[0].hasAssignedBounds
        && out[1].exponent == 7873417  && out[1].hasFactoredTo
        && out[2].exponent == 28733813 && out[2].hasAssignedBounds;
    const bool consumed = loaded && consumeWorktodoEntry(TEST_FILE, out[1], err);
    vector<WorktodoEntry> after;
    const bool ok = consumed && loadWorktodo(TEST_FILE, after, err) && after.size() == 2
        && after[0].exponent == 786613 && after[1].exponent == 28733813;
    if (!ok) { ++fails; printf("  FAIL mixed file: %s\n", err.c_str()); }
    printf("     %s  three line shapes coexist, consume targets the right physical line\n", ok ? "PASS" : "FAIL");
  }

  filesystem::remove(TEST_FILE);
  filesystem::remove(string(TEST_FILE) + ".tmp");

  printf("\n%s\n", fails ? "worktodo: FAILED" : "worktodo: all tests passed");
  return fails ? 1 : 0;
}
