// Copyright (C) Mp_p-1_gpu

#include "Worktodo.h"
#include "BigInt.h"
#include "Config.h"
#include "Bounds.h"
#include "Stage2Plan.h"

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
//
// One WHOLE line, terminator and all, however long. A plain
// fgets(buf, 1024, f) stops at its buffer and hands the remainder back as
// though it were the next line, which both loops below then parsed as an
// entry in its own right: a comment longer than 1023 bytes whose tail
// happened to be digits queued an exponent that appears nowhere in the file
// and wrote its result to results.txt under that exponent. Consuming the
// phantom was worse -- the truncated head went back without the newline it
// never had, gluing the next real entry onto the end of a comment and
// dropping it from the queue with nothing said. Both loops share this, so
// their line numbering and what consume copies through cannot disagree.
bool readRawLine(FILE* f, string& raw) {
  raw.clear();
  char buf[1024];
  while (fgets(buf, sizeof(buf), f)) {
    raw += buf;
    if (!raw.empty() && raw.back() == '\n') { return true; }
  }
  return !raw.empty();   // last line, no trailing newline: still a line
}

string stripEol(const string& line) {
  string raw = line;
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
// factor strings, trimmed, validated as decimal integers that actually
// divide M_p.
bool parseKnownFactors(const string& field, u32 exponent, vector<string>& out, string& what) {
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
    // Being made of digits is not enough. This list is echoed verbatim into
    // results.txt as "known-factors", where it is a CLAIM about M_p that
    // PrimeNet reads -- and it is the ONLY thing the program does with it, so
    // a typo here had nothing downstream that could notice. 2^p == 1 (mod f)
    // is exactly "f divides 2^p - 1"; it costs microseconds, and it is the
    // same test PM1.cpp's describe() already applies to every factor this
    // program reports on its own account.
    Nat f;
    if (!fromDecimal(one, f) || f < Nat(3) || !powMod(Nat(2), Nat(exponent), f).isOne()) {
      what = "known factor " + one + " does not divide M" + to_string(exponent);
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
// A plain decimal integer, for the fields where strtod's flexibility is a
// hazard rather than a convenience. Before 1.9.1 every numeric field went
// through strtod, which accepts far more than the format specifies:
// "1000099.7" rounded to a DIFFERENT exponent (M1000100) and ran, "0x10"
// became M16, "1e9" became M1000000000, and "nan" defeated every range check
// (NaN compares false against everything) to arrive as M0. Underscores are
// tolerated as digit separators, matching config.txt's own convention;
// nothing else is. Prime95 and AutoPrimeNet always write plain integers here,
// so this is stricter without being incompatible.
bool parseExactU64(const string& s, u64& out) {
  string d;
  for (char c : s) { if (c != '_') { d += c; } }
  if (d.empty() || d.find_first_not_of("0123456789") != string::npos) { return false; }
  u64 v = 0;
  for (char c : d) {
    const u64 digit = u64(c - '0');
    if (v > (~0ull - digit) / 10) { return false; }   // would overflow u64
    v = v * 10 + digit;
  }
  out = v;
  return true;
}

// M_p is only a GIMPS candidate for prime p, and config.txt has always said
// so -- but nothing enforced it until 1.9.1, so a typo'd composite exponent
// ran to completion. That is not merely wasted time: for composite p every
// d | p contributes an algebraic factor 2^d - 1 of 2^p - 1, so such a run can
// "find" a factor and write it to results.txt as though it were a discovery.
// Trial division is instant at u32 range (at most ~65536 steps).
bool isPrimeExponent(u64 n) {
  if (n < 2) { return false; }
  if (n % 2 == 0) { return n == 2; }
  for (u64 d = 3; d * d <= n; d += 2) {
    if (n % d == 0) { return false; }
  }
  return true;
}

// Shared by both entry shapes: the exponent must be a plain integer, in range,
// and prime. `what` is the raw token, quoted back in every message so a bad
// line names itself.
bool validateExponent(const string& what, u32 lineNo, u64 n, u32& exponent, string& err) {
  const string where = "worktodo.txt line " + to_string(lineNo) + ": exponent " + what;
  if (n < 3) { err = where + " is too small"; return false; }
  if (n > 4294967295ull) {
    err = where + " is above this program's limit of 4294967295";
    return false;
  }
  if (!isPrimeExponent(n)) {
    err = where + " is not prime (M_p = 2^p - 1 is only a candidate for prime p)";
    return false;
  }
  exponent = u32(n);
  return true;
}

bool parseKbnc(const vector<string>& f, u32 lineNo, u32& exponent, string& err) {
  if (f.size() < 4) {
    err = "worktodo.txt line " + to_string(lineNo) + ": expected k,b,n,c fields";
    return false;
  }
  // k, b and c stay on the strtod path -- they are compared against exact
  // constants below, so a surprising parse is rejected rather than acted on.
  // n does NOT: it is the one field whose value is used rather than checked,
  // so it gets the strict integer parse. See parseExactU64.
  double k = 0, b = 0, c = 0;
  u64 n = 0;
  if (!parseDoubleField(f[0], k) || !parseDoubleField(f[1], b)
      || !parseDoubleField(f[3], c)) {
    err = "worktodo.txt line " + to_string(lineNo) + ": malformed k,b,n,c";
    return false;
  }
  if (!parseExactU64(f[2], n)) {
    err = "worktodo.txt line " + to_string(lineNo) + ": exponent " + f[2]
        + " is not a plain decimal integer";
    return false;
  }
  // Written as a range test, not "b + 0.5 != 2": u64(b + 0.5) accepted
  // anything in [1.5, 2.5) as base 2, and NaN reached the cast at all.
  if (k != 1.0 || !(b >= 2.0 && b <= 2.0) || c != -1.0) {
    err = "worktodo.txt line " + to_string(lineNo) + ": k=" + f[0] + " b=" + f[1]
        + " c=" + f[3] + " is not a Mersenne number (this program only factors M_p = 2^p - 1)";
    return false;
  }
  return validateExponent(f[2], lineNo, n, exponent, err);
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
    if (!parseKnownFactors(f[6], out.exponent, out.knownFactors, what)) {
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
  if (!parseDoubleField(f[4], b1) || !parseDoubleField(f[5], b2)
      || !(b1 >= 2) || !(b2 >= 2)) {
    err = "worktodo.txt line " + to_string(lineNo) + ": malformed B1/B2";
    return false;
  }
  // Positive-form range tests above, and an explicit ceiling here. "b1 < 2"
  // let NaN straight through (NaN < 2 is false), and an out-of-range B2 was
  // worse than useless: chooseBounds discards every candidate above 4e9, its
  // candidate list comes back empty, and the default-constructed Bounds it
  // then returns means B1 = 0 -- a run that does no stage 1 at all and
  // reports "no factor" for an assignment that asked for real work.
  if (!(b2 <= STAGE2_B2_MAX)) {
    err = "worktodo.txt line " + to_string(lineNo) + ": B2 " + f[5]
        + " is above this program's limit of " + to_string(u64(STAGE2_B2_MAX))
        + " (stage 2's setup exponents must fit a machine word)";
    return false;
  }
  out.assignedB1 = u64(b1 + 0.5);
  out.assignedB2 = u64(b2 + 0.5);
  if (out.assignedB2 <= out.assignedB1) {
    err = "worktodo.txt line " + to_string(lineNo) + ": B2 must be greater than B1";
    return false;
  }
  out.hasAssignedBounds = true;

  // A Pminus1= line always asks for a stage 2 (B2 > B1 is enforced just
  // above), and stage 2's pairing needs every prime it walks to exceed w*D/2
  // -- with the smallest shape this program has, D = 210 and w = 1, that
  // floor is 105. An assigned B1 under it is not a small job but an
  // impossible one: buildStage2Plan throws, and it throws only once stage 1
  // has already run to completion, leaving the entry queued to fail the same
  // way on every restart. Refuse it here, where the line can name itself.
  if (out.assignedB1 < STAGE2_MIN_B1) {
    err = "worktodo.txt line " + to_string(lineNo) + ": B1 " + f[4]
        + " is below " + to_string(STAGE2_MIN_B1) + ", the smallest B1 stage 2 can pair against";
    return false;
  }

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
    if (!parseKnownFactors(f[cursor], out.exponent, out.knownFactors, what)) {
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
    if (!parseExactU64(t, n)) {
      err = "worktodo.txt line " + to_string(lineNo) + ": '" + t
          + "' is not a plain decimal integer";
      return false;
    }
    return validateExponent(t, lineNo, n, out.exponent, err);
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

  string raw;
  u32 lineNo = 0;
  while (readRawLine(f, raw)) {
    ++lineNo;
    const string t = trim(stripEol(raw));
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

ResolvedBounds resolveBounds(const WorktodoEntry& entry, u64 configuredB1, u64 configuredB2) {
  if (entry.hasAssignedBounds) { return {entry.assignedB1, entry.assignedB2}; }
  return {configuredB1, configuredB2};
}

bool consumeWorktodoEntry(const string& path, const WorktodoEntry& entry, string& err) {
  FILE* fi = fopen(path.c_str(), "r");
  if (!fi) { err = "worktodo.txt disappeared"; return false; }

  const string tmp = path + ".tmp";
  FILE* fo = fopen(tmp.c_str(), "w");
  if (!fo) { fclose(fi); err = "cannot open " + tmp; return false; }

  string raw;
  u32 lineNo = 0;
  bool removed = false;
  bool ok = true;
  while (readRawLine(fi, raw)) {
    ++lineNo;
    if (lineNo == entry.lineNo) {
      // Confirm this is still the entry we loaded before dropping it -- the
      // file may have been hand-edited, or reissued by AutoPrimeNet under a
      // new assignment, since. On a mismatch, fall through and copy the line
      // through unchanged instead of removing it.
      WorktodoEntry reread;
      string perr;
      const string t = trim(stripEol(raw));
      const bool sameEntry = !t.empty() && t[0] != '#'
          && parseWorktodoLine(t, lineNo, reread, perr)
          && reread.exponent == entry.exponent
          && (entry.aid.empty() || reread.aid == entry.aid);
      if (sameEntry) {
        removed = true;
        continue;
      }
    }
    if (fwrite(raw.data(), 1, raw.size(), fo) != raw.size()) { ok = false; }
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
    // The exponent here used to be 86243 -- which is a Mersenne PRIME, so
    // every known-factors fixture written against it was necessarily made
    // up. That was harmless while nothing checked them and wrong the moment
    // something did. M4444091 is composite and these are real factors of it.
    struct Case { const char* line; bool wantAid; bool wantKF; };
    const Case cases[] = {
      {"Pfactor=1,2,4444091,-1,74.5,1.30", false, false},
      {"Pfactor=AABBCCDD00112233445566778899AABB,1,2,4444091,-1,74.5,1.30", true, false},
      {"Pfactor=1,2,4444091,-1,74.5,1.30,\"636358278473\"", false, true},
      {"Pfactor=AABBCCDD00112233445566778899AABB,1,2,4444091,-1,74.5,1.30,\"8888183,319974553\"", true, true},
    };
    bool allOk = true;
    for (const Case& c : cases) {
      writeRaw(TEST_FILE, string(c.line) + "\n");
      vector<WorktodoEntry> out;
      string err;
      const bool ok = loadWorktodo(TEST_FILE, out, err) && out.size() == 1
          && out[0].exponent == 4444091 && out[0].hasFactoredTo && out[0].factoredTo == 74
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

    writeRaw(TEST_FILE, "Pminus1=1,2,4444091,-1,50000,3000000\n");
    { vector<WorktodoEntry> out; string err;
      const bool ok = loadWorktodo(TEST_FILE, out, err) && out.size() == 1
          && out[0].hasAssignedBounds && out[0].assignedB1 == 50000 && out[0].assignedB2 == 3000000
          && !out[0].hasFactoredTo && out[0].b2Start == 0;
      if (!ok) { allOk = false; printf("  FAIL Pminus1 bare: %s\n", err.c_str()); } }

    writeRaw(TEST_FILE, "Pminus1=1,2,4444091,-1,50000,3000000,74.0\n");
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
    writeRaw(TEST_FILE, "Pminus1=1,2,4444091,-1,50000,3000000,3100000\n");
    { vector<WorktodoEntry> out; string err;
      const bool ok = loadWorktodo(TEST_FILE, out, err) && out.size() == 1
          && !out[0].hasFactoredTo && out[0].b2Start == 0;
      if (!ok) { allOk = false; printf("  FAIL Pminus1 B2_start-in-sieve-slot: %s\n", err.c_str()); } }

    // Same slot, but a value a real typo could plausibly produce (200, not
    // an astronomical B2_start-shaped number) -- still > 127, still fails
    // sieve_depth, still discarded rather than misread as something else.
    writeRaw(TEST_FILE, "Pminus1=1,2,4444091,-1,50000,3000000,200\n");
    { vector<WorktodoEntry> out; string err;
      const bool ok = loadWorktodo(TEST_FILE, out, err) && out.size() == 1
          && !out[0].hasFactoredTo && out[0].b2Start == 0;
      if (!ok) { allOk = false; printf("  FAIL Pminus1 mid-range value in sieve-slot: %s\n", err.c_str()); } }

    writeRaw(TEST_FILE,
             "Pminus1=AABBCCDD00112233445566778899AABB,1,2,4444091,-1,50000,3000000,74.0,3100000,\"636358278473\"\n");
    { vector<WorktodoEntry> out; string err;
      const bool ok = loadWorktodo(TEST_FILE, out, err) && out.size() == 1
          && out[0].aid == "AABBCCDD00112233445566778899AABB"
          && out[0].hasFactoredTo && out[0].factoredTo == 74
          && out[0].b2Start == 3100000 && out[0].knownFactors.size() == 1
          && out[0].knownFactors[0] == "636358278473";
      if (!ok) { allOk = false; printf("  FAIL Pminus1 full: %s\n", err.c_str()); } }

    writeRaw(TEST_FILE, "Pminus1=1,2,4444091,-1,50000,3000000,74.0,40000\n");
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

  printf("\n  I. rejected shapes: wrong k/b/c, malformed known-factors\n");
  {
    bool allOk = true;
    const char* bad[] = {
      "Pminus1=2,2,86243,-1,50000,3000000",              // k != 1
      "Pminus1=1,3,86243,-1,50000,3000000",               // b != 2
      "Pminus1=1,2,86243,1,50000,3000000",                // c != -1
      "Pfactor=1,2,4444091,-1,74.0,1.3,\"8888183",       // unterminated quote
      "Pfactor=1,2,4444091,-1,74.0,1.3,\"8888183,abc\"",  // non-decimal entry
      "Pfactor=1,2,4444091,-1,74.0,1.3,\"8888183,,319974553\"",  // empty entry between two real ones
    };
    for (const char* line : bad) {
      writeRaw(TEST_FILE, string(line) + "\n");
      vector<WorktodoEntry> out; string err;
      if (loadWorktodo(TEST_FILE, out, err)) { allOk = false; printf("  FAIL accepted bad line '%s'\n", line); }
    }
    if (!allOk) { ++fails; }
    printf("     %s  non-Mersenne k/b/c and malformed known-factors are hard errors\n", allOk ? "PASS" : "FAIL");
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

  // Real bug this guards: the queue loop's bounds precedence has already
  // changed direction twice in one session (gated behind a config switch
  // that defaulted to ignoring the assignment, then that switch removed in
  // favor of always honoring it) and, before this test existed, had zero
  // automated coverage either way -- a regression here would silently run
  // the wrong B1/B2 against a real assignment and ship undetected.
  printf("\n  L. resolveBounds: assignment wins when present, config.txt wins otherwise\n");
  {
    WorktodoEntry withBounds; withBounds.hasAssignedBounds = true;
    withBounds.assignedB1 = 700; withBounds.assignedB2 = 10500;
    WorktodoEntry bare;   // hasAssignedBounds == false, the bare-exponent/Pfactor= shape

    const ResolvedBounds a = resolveBounds(withBounds, /*configuredB1=*/999, /*configuredB2=*/999);
    const ResolvedBounds b = resolveBounds(bare, /*configuredB1=*/12345, /*configuredB2=*/67890);
    const ResolvedBounds c = resolveBounds(bare, /*configuredB1=*/0, /*configuredB2=*/0);

    const bool ok = a.b1 == 700 && a.b2 == 10500      // assignment wins even though config differs
        && b.b1 == 12345 && b.b2 == 67890             // no assignment -> config.txt's own value
        && c.b1 == 0 && c.b2 == 0;                    // no assignment, config is "auto" -> stays auto
    if (!ok) {
      ++fails;
      printf("  FAIL resolveBounds: a={%llu,%llu} b={%llu,%llu} c={%llu,%llu}\n",
             (unsigned long long) a.b1, (unsigned long long) a.b2,
             (unsigned long long) b.b1, (unsigned long long) b.b2,
             (unsigned long long) c.b1, (unsigned long long) c.b2);
    }
    printf("     %s  assignment bounds override config.txt; a bare/Pfactor= entry falls back to it\n",
           ok ? "PASS" : "FAIL");
  }

  filesystem::remove(TEST_FILE);
  filesystem::remove(string(TEST_FILE) + ".tmp");

  printf("\n  M. Pfactor= corpus: a real assignment sweep, and out-of-range refused\n");
  {
    // A PrimeNet-shaped sweep: 60 assignments spanning 1e6 to 1e10 and TF 67
    // to 96 bits. Section F covers the SHAPE of a Pfactor= line (AID, known
    // factors); this covers its RANGE, which is where the parser was actually
    // wrong -- before 1.9 every exponent above 2^32-1 wrapped into a
    // different, entirely plausible exponent and ran. The five over-range
    // entries are the point of this corpus, not an afterthought: they must be
    // REFUSED, with an error, never silently reinterpreted.
    struct C { u64 exponent; u32 tf; };
    static const C CORPUS[] = {
      { 1000099u,    67 },
      { 3000077u,    67 },
      { 4444091u,    70 },
      { 10000831u,   68 },
      { 24000577u,   70 },
      { 30000853u,   67 },
      { 50001781u,   74 },
      { 51558151u,   74 },
      { 54447193u,   74 },
      { 58610467u,   74 },
      { 61012769u,   74 },
      { 81229789u,   75 },
      { 100000081u,  76 },
      { 110505011u,  76 },
      { 120002191u,  77 },
      { 150000713u,  77 },
      { 200001187u,  79 },
      { 230086243u,  79 },
      { 249500501u,  79 },
      { 290001377u,  80 },
      { 300008497u,  80 },
      { 301000159u,  80 },
      { 332230189u,  81 },
      { 353466917u,  81 },
      { 407363239u,  81 },
      { 423000089u,  82 },
      { 464000021u,  82 },
      { 490000003u,  82 },
      { 502000027u,  82 },
      { 563021377u,  81 },
      { 640402457u,  84 },
      { 654036583u,  84 },
      { 745964951u,  84 },
      { 840859433u,  83 },
      { 901000031u,  85 },
      { 940216091u,  85 },
      { 980000521u,  78 },
      { 1100000081u, 68 },
      { 1138000001u, 80 },
      { 1199999579u, 68 },
      { 1299999919u, 68 },
      { 1414000001u, 87 },
      { 1552999957u, 87 },
      { 1553000003u, 87 },
      { 1553000143u, 87 },
      { 1708000339u, 88 },
      { 1862000159u, 88 },
      { 2000000089u, 89 },
      { 3321928601u, 91 },
      { 3321928619u, 91 },
      { 3321928703u, 91 },
      { 3321928787u, 91 },
      { 3330000293u, 91 },
      { 4000000229u, 92 },
      { 4294967111u, 92 },
      { 4294967357u, 92 },
      { 6330000829u, 94 },
      { 8230000949u, 96 },
      { 8883334777u, 96 },
      { 9993000001u, 96 },
    };
    constexpr u64 MAX_EXP = 4294967295ull;
    bool allOk = true;
    u32 accepted = 0, refused = 0;
    for (const C& c : CORPUS) {
      const string line = "Pfactor=1,2," + to_string(c.exponent) + ",-1,"
                        + to_string(c.tf) + ",2";
      writeRaw(TEST_FILE, line + "\n");
      vector<WorktodoEntry> out;
      string err;
      const bool parsed = loadWorktodo(TEST_FILE, out, err) && out.size() == 1;
      const bool inRange = c.exponent <= MAX_EXP;
      bool ok;
      if (inRange) {
        ok = parsed && out[0].exponent == u32(c.exponent)
             && out[0].hasFactoredTo && out[0].factoredTo == c.tf;
        if (ok) { ++accepted; }
      } else {
        // Refused is the ONLY acceptable outcome here. Parsing "successfully"
        // into some other exponent is precisely the bug this catches.
        ok = !parsed;
        if (ok) {
          ++refused;
        } else if (out.size() == 1) {
          printf("  FAIL exponent %llu silently became M%u\n",
                 (unsigned long long) c.exponent, out[0].exponent);
        }
      }
      if (!ok) {
        allOk = false;
        printf("  FAIL Pfactor corpus '%s': %s\n", line.c_str(),
               err.empty() ? "unexpected outcome" : err.c_str());
      }
    }
    if (!allOk) { ++fails; }
    printf("     %s  %u in-range parsed, %u over-range refused (of %zu)\n",
           allOk ? "PASS" : "FAIL", accepted, refused,
           sizeof(CORPUS) / sizeof(CORPUS[0]));
  }

  printf("\n  N. inputs that must be REFUSED, not reinterpreted\n");
  {
    // Every entry here was, at some point, accepted and silently turned into a
    // different job. The bare-integer shape is listed alongside the assignment
    // shape on purpose: 1.9 guarded parseKbnc and shipped with the bare path
    // still wrapping, because the corpus above only ever exercised Pfactor=.
    struct Bad { const char* line; const char* why; };
    static const Bad BAD[] = {
      // over-range exponent -- used to wrap modulo 2^32
      { "9993000001",                              "bare: 9993000001 -> M1403065409" },
      { "Pfactor=1,2,9993000001,-1,96,2",          "Pfactor: same wrap" },
      { "Pminus1=1,2,4294967357,-1,100000,2000000","Pminus1: 4294967357 -> M61" },
      // non-integer exponent -- used to round to a different exponent
      { "1000099.7",                               "bare: rounds to M1000100" },
      { "Pfactor=1,2,1000099.7,-1,74,2",           "Pfactor: rounds to M1000100" },
      { "1e9",                                     "bare: strtod exponent form" },
      { "0x10",                                    "bare: strtod hex form" },
      { "Pfactor=1,2,0x10,-1,74,2",                "Pfactor: strtod hex form" },
      // NaN -- compares false against every bound, so every range test missed it
      { "Pfactor=1,2,nan,-1,74,2",                 "NaN exponent -> M0" },
      { "Pminus1=1,2,1000099,-1,nan,2000000",      "NaN B1" },
      { "Pminus1=1,2,1000099,-1,100000,nan",       "NaN B2 -> B1 = 0 run" },
      // composite exponent -- 2^p-1 has algebraic factors, so a run can
      // "find" one and report it as a discovery
      { "1000100",                                 "bare: composite exponent" },
      { "Pfactor=1,2,1000100,-1,74,2",             "Pfactor: composite exponent" },
      { "Pfactor=1,2,86244,-1,74,2",               "Pfactor: composite exponent" },
      // B2 the stage-2 engine cannot represent -- used to yield B1 = 0 and a
      // completed "no factor" result for an assignment that asked for work
      { "Pminus1=1,2,1000099,-1,100000,1e18",      "B2 above the machine-word cap" },
      { "Pminus1=1,2,1000099,-1,100000,5000000000","B2 just above the cap" },
      // base that is not 2 -- u64(b + 0.5) accepted anything in [1.5, 2.5)
      { "Pfactor=1,2.4,1000099,-1,74,2",           "b=2.4 rounded to base 2" },
      { "Pfactor=1,1.6,1000099,-1,74,2",           "b=1.6 rounded to base 2" },
    };
    bool allOk = true;
    for (const Bad& b : BAD) {
      writeRaw(TEST_FILE, string(b.line) + "\n");
      vector<WorktodoEntry> out;
      string err;
      const bool parsed = loadWorktodo(TEST_FILE, out, err) && !out.empty();
      // Refusal must also SAY something -- a silently dropped line is how
      // several of these hid for as long as they did.
      const bool ok = !parsed && !err.empty();
      if (!ok) {
        allOk = false;
        if (parsed) {
          printf("  FAIL accepted '%s' as M%u (%s)\n", b.line, out[0].exponent, b.why);
        } else {
          printf("  FAIL refused '%s' with no diagnostic (%s)\n", b.line, b.why);
        }
      }
    }
    if (!allOk) { ++fails; }
    printf("     %s  %zu malformed/out-of-range inputs each refused with a reason\n",
           allOk ? "PASS" : "FAIL", sizeof(BAD) / sizeof(BAD[0]));
  }

  printf("\n  O. lines longer than the read buffer\n");
  {
    // Both loops used to read with fgets(buf, 1024, f), which hands back the
    // tail of a longer line as though it were the next line -- and the tail of
    // this comment is digits, so it queued M1000003, an exponent that appears
    // nowhere in the file, and wrote its result to results.txt under it.
    const string longComment = "# " + string(2000, 'x') + "1000003";
    writeRaw(TEST_FILE, longComment + "\n786613\n");
    vector<WorktodoEntry> out;
    string err;
    const bool loaded = loadWorktodo(TEST_FILE, out, err);
    const bool oneEntry = loaded && out.size() == 1 && out[0].exponent == 786613;
    if (!oneEntry) {
      ++fails;
      printf("  FAIL long comment: ok=%d size=%zu first=%u\n", loaded, out.size(),
             out.empty() ? 0 : out[0].exponent);
    }
    printf("     %s  a %zu-byte comment is one comment, not a comment plus M1000003\n",
           oneEntry ? "PASS" : "FAIL", longComment.size());

    // The consume half had the worse end of the same bug: the truncated head
    // went back to disk without the newline it never had, so the next real
    // entry was glued onto the end of a comment and left the queue silently.
    bool consumeOk = false;
    if (oneEntry) {
      string cerr;
      vector<WorktodoEntry> after;
      string aerr;
      consumeOk = consumeWorktodoEntry(TEST_FILE, out[0], cerr)
               && loadWorktodo(TEST_FILE, after, aerr) && after.empty();
      // ... and the comment it left behind is still exactly one line.
      size_t newlines = 0;
      if (FILE* fp = fopen(TEST_FILE, "rb")) {
        char c;
        while (fread(&c, 1, 1, fp) == 1) { if (c == '\n') { ++newlines; } }
        fclose(fp);
      }
      consumeOk = consumeOk && newlines == 1;
    }
    if (!consumeOk) { ++fails; }
    printf("     %s  consuming the entry after it leaves the long line whole\n",
           consumeOk ? "PASS" : "FAIL");

    // The shape that makes this reachable in a real file is a Pminus1= line
    // whose known-factors list runs past the buffer. These three are genuine
    // factors of M4444091; the padding between them stands in for the dozens
    // more a heavily-factored exponent carries, without inventing any.
    const string many = "Pminus1=1,2,4444091,-1,100000,3000000,70,0,"
                        "\"8888183," + string(400, ' ')
                        + "319974553," + string(400, ' ')
                        + "636358278473" + string(300, ' ') + "\"";
    writeRaw(TEST_FILE, many + "\n");
    vector<WorktodoEntry> big;
    string berr;
    const bool bigOk = loadWorktodo(TEST_FILE, big, berr) && big.size() == 1
                    && big[0].exponent == 4444091 && big[0].knownFactors.size() == 3;
    if (!bigOk) {
      ++fails;
      printf("  FAIL long Pminus1: size=%zu err=%s\n", big.size(), berr.c_str());
    }
    printf("     %s  a %zu-byte Pminus1= line is one entry with all its known factors\n",
           bigOk ? "PASS" : "FAIL", many.size());
  }

  printf("\n  P. known factors are checked against M_p, not just against the digits\n");
  {
    // This list is echoed verbatim into results.txt as "known-factors", where
    // PrimeNet reads it as a claim about M_p -- and it was the only thing the
    // program ever did with it, so a typo had nothing downstream to catch it.
    struct KF { const char* line; bool want; const char* why; };
    static const KF KFS[] = {
      { "Pfactor=1,2,11,-1,10,1.0,\"23,89\"",                              true,
        "2047 = 23 * 89" },
      { "Pfactor=1,2,11,-1,10,1.0,\"23,91\"",                              false,
        "91 does not divide 2047" },
      { "Pfactor=1,2,1000099,-1,67,1.6,\"155058493988826487335266033969\"", true,
        "a real 30-digit factor of M1000099" },
      { "Pfactor=1,2,1000099,-1,67,1.6,\"155058493988826487335266033967\"", false,
        "the same factor with its last digit typed wrong" },
      { "Pfactor=1,2,11,-1,10,1.0,\"1\"",                                  false,
        "1 divides everything and factors nothing" },
      { "Pfactor=1,2,11,-1,10,1.0,\"0\"",                                  false,
        "0 would have been a modulus of zero" },
      { "Pminus1=1,2,11,-1,100000,3000000,10,0,\"89\"",                    true,
        "Pminus1= reaches the same check through its positional tail" },
      { "Pminus1=1,2,11,-1,100000,3000000,10,0,\"87\"",                    false,
        "and refuses there too" },
    };
    bool allOk = true;
    for (const KF& k : KFS) {
      writeRaw(TEST_FILE, string(k.line) + "\n");
      vector<WorktodoEntry> out;
      string err;
      const bool got = loadWorktodo(TEST_FILE, out, err) && !out.empty();
      if (got != k.want || (!got && err.empty())) {
        allOk = false;
        printf("  FAIL '%s': %s (%s)\n", k.line,
               got ? "accepted, should be refused" : "refused, should be accepted", k.why);
      }
    }
    if (!allOk) { ++fails; }
    printf("     %s  %zu known-factor lists, each accepted or refused on divisibility\n",
           allOk ? "PASS" : "FAIL", sizeof(KFS) / sizeof(KFS[0]));
  }

  printf("\n  Q. an assigned B1 stage 2 cannot pair against\n");
  {
    // A Pminus1= line always asks for a stage 2, and the smallest shape
    // (D=210, w=1) cannot walk primes at or below 105. Below that,
    // buildStage2Plan throws -- but only once stage 1 has already run, and
    // again on every restart, because the entry is still in the queue.
    struct B1C { const char* line; bool want; };
    static const B1C B1S[] = {
      { "Pminus1=1,2,1000099,-1,10,1000",     false },
      { "Pminus1=1,2,1000099,-1,104,200000",  false },
      { "Pminus1=1,2,1000099,-1,105,200000",  true  },
      { "Pminus1=1,2,1000099,-1,100000,3000000", true },
    };
    bool allOk = true;
    for (const B1C& b : B1S) {
      writeRaw(TEST_FILE, string(b.line) + "\n");
      vector<WorktodoEntry> out;
      string err;
      const bool got = loadWorktodo(TEST_FILE, out, err) && !out.empty();
      if (got != b.want || (!got && err.empty())) {
        allOk = false;
        printf("  FAIL '%s': %s\n", b.line,
               got ? "accepted, should be refused" : "refused, should be accepted");
      }
    }
    if (!allOk) { ++fails; }
    printf("     %s  B1 below %llu refused at parse time, not an hour into stage 2\n",
           allOk ? "PASS" : "FAIL", (unsigned long long) STAGE2_MIN_B1);
  }

  filesystem::remove(TEST_FILE);
  filesystem::remove(string(TEST_FILE) + ".tmp");

  printf("\n%s\n", fails ? "worktodo: FAILED" : "worktodo: all tests passed");
  return fails ? 1 : 0;
}
