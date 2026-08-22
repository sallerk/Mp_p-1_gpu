// Copyright (C) Mp_p-1_gpu

#include "Config.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <string>

using namespace std;

namespace {

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

// True when `v` is a spelling of `want`. Used where the two truth values need
// separate tests rather than one out-parameter.
bool parseBoolValue(const string& v, bool want) {
  bool got = false;
  const string t = lower(v);
  if (t == "1" || t == "true"  || t == "yes" || t == "on")  { got = true; }
  else if (t == "0" || t == "false" || t == "no" || t == "off") { got = false; }
  else { return false; }
  return got == want;
}

bool parseBool(const string& v, bool& out) {
  const string t = lower(v);
  if (t == "1" || t == "true"  || t == "yes" || t == "on")  { out = true;  return true; }
  if (t == "0" || t == "false" || t == "no"  || t == "off") { out = false; return true; }
  return false;
}

} // namespace

bool parseNumber(const string& in, u64& out) {
  string s;
  for (char c : in) { if (c != '_' && c != ',') { s += c; } }
  if (s.empty()) { return false; }

  u64 mult = 1;
  const char last = (char) toupper((unsigned char) s.back());
  if (last == 'K') { mult = 1000; s.pop_back(); }
  else if (last == 'M') { mult = 1000000; s.pop_back(); }
  else if (last == 'G') { mult = 1000000000ull; s.pop_back(); }
  if (s.empty()) { return false; }

  char* end = nullptr;
  const double v = strtod(s.c_str(), &end);
  if (end != s.c_str() + s.size() || v < 0) { return false; }
  out = u64(v * double(mult) + 0.5);
  return true;
}

u32 defaultFactoredTo(u32 exponent) {
  // Roughly the depth GIMPS trial-factors to at each size. Only a fallback;
  // the real value for a given exponent comes from mersenne.org.
  if (exponent <  20000000u) { return 69; }
  if (exponent <  40000000u) { return 72; }
  if (exponent <  60000000u) { return 74; }
  if (exponent <  80000000u) { return 75; }
  if (exponent < 100000000u) { return 76; }
  if (exponent < 200000000u) { return 77; }
  return 78;
}

bool loadConfig(const string& path, Config& cfg, string& err) {
  FILE* f = fopen(path.c_str(), "r");
  if (!f) { err = "cannot open '" + path + "'"; return false; }

  char buf[1024];
  int lineNo = 0;
  while (fgets(buf, sizeof(buf), f)) {
    ++lineNo;
    string line = buf;
    const size_t hash = line.find('#');
    if (hash != string::npos) { line = line.substr(0, hash); }
    line = trim(line);
    if (line.empty()) { continue; }

    const size_t eq = line.find('=');
    if (eq == string::npos) {
      err = "line " + to_string(lineNo) + ": expected key = value";
      fclose(f);
      return false;
    }
    const string key = lower(trim(line.substr(0, eq)));
    const string val = trim(line.substr(eq + 1));
    const string lval = lower(val);

    auto bad = [&](const string& what) {
      err = "line " + to_string(lineNo) + ": " + key + " " + what;
      fclose(f);
      return false;
    };

    u64 n = 0;
    if (key == "exponent") {
      err = "line " + to_string(lineNo)
          + ": exponent moved to worktodo.txt in 1.4 -- put p there instead, one per line";
      fclose(f);
      return false;
    } else if (key == "method") {
      if      (lval == "pm1")  { cfg.doPM1 = true;  cfg.doPP1 = false; }
      else if (lval == "pp1" || lval == "p+1") { cfg.doPM1 = false; cfg.doPP1 = true; }
      else if (lval == "both") { cfg.doPM1 = true;  cfg.doPP1 = true; }
      else { return bad("must be pm1, pp1 or both"); }
    } else if (key == "pp1_seeds" || key == "seeds") {
      cfg.pp1Seeds.clear();
      size_t pos = 0;
      while (pos <= val.size()) {
        const size_t comma = val.find(',', pos);
        const string one = trim(val.substr(pos, comma == string::npos ? string::npos : comma - pos));
        if (!one.empty()) {
          u64 sv = 0;
          if (!parseNumber(one, sv) || sv < 3 || sv > 1000000) {
            return bad("seeds must be 3..1000000 (2 is degenerate: V_n(2,1) == 2)");
          }
          cfg.pp1Seeds.push_back(u32(sv));
        }
        if (comma == string::npos) { break; }
        pos = comma + 1;
      }
      if (cfg.pp1Seeds.empty()) { return bad("needs at least one seed"); }
    } else if (key == "b1") {
      if (lval == "auto") { cfg.b1 = 0; }
      else if (!parseNumber(val, n) || n < 2) { return bad("must be >= 2 or 'auto'"); }
      else { cfg.b1 = n; }
    } else if (key == "b2") {
      if (lval == "auto") { cfg.b2 = 0; }
      else if (!parseNumber(val, n) || n < 2) { return bad("must be >= 2 or 'auto'"); }
      else { cfg.b2 = n; }
    } else if (key == "stages") {
      // Counts a number of stages, so "1" means stage 1 only and "2"/"both"
      // means run both.
      if      (lval == "auto") { cfg.stage2Mode = STAGE2_AUTO; }
      else if (lval == "both" || lval == "2" || lval == "12" || lval == "1+2") {
        cfg.stage2Mode = STAGE2_ON;
      } else if (lval == "1" || lval == "one") {
        cfg.stage2Mode = STAGE2_OFF;
      } else { return bad("must be auto, both, or 1"); }
    } else if (key == "stage2") {
      // A yes/no switch on stage 2, so "1" means ON -- the opposite of what it
      // means to `stages`. Kept as separate keys precisely so neither reads
      // backwards.
      if      (lval == "auto") { cfg.stage2Mode = STAGE2_AUTO; }
      else if (parseBoolValue(lval, true))  { cfg.stage2Mode = STAGE2_ON; }
      else if (parseBoolValue(lval, false)) { cfg.stage2Mode = STAGE2_OFF; }
      else { return bad("must be auto, yes or no"); }
    } else if (key == "stage2_d") {
      if (lval == "auto") { cfg.stage2D = 0; }
      else if (!parseNumber(val, n) || n < 2 || n > 100000 || (n & 1)) {
        return bad("must be an even number 2..100000 or 'auto'");
      } else { cfg.stage2D = u32(n); }
    } else if (key == "stage2_w" || key == "stage2_window") {
      if (lval == "auto") { cfg.stage2W = 0; }
      else if (!parseNumber(val, n) || n < 1 || n > 99) { return bad("must be 1..99 or 'auto'"); }
      else { cfg.stage2W = u32(n); }
    } else if (key == "factored_to" || key == "factoredto") {
      if (lval == "auto") { cfg.factoredTo = 0; }
      else if (!parseNumber(val, n) || n < 1 || n > 127) { return bad("must be 1..127 or 'auto'"); }
      else { cfg.factoredTo = u32(n); }
    } else if (key == "bias") {
      cfg.bias = atof(val.c_str());
      if (!(cfg.bias > 0)) { return bad("must be > 0"); }
    } else if (key == "bounds_tolerance") {
      cfg.boundsTolerance = atof(val.c_str());
      if (!(cfg.boundsTolerance >= 0) || cfg.boundsTolerance > 1) {
        return bad("must be 0..1 (a fraction of expected cost)");
      }
    } else if (key == "fft") {
      cfg.fftSpec = (lval == "auto") ? "" : val;
    } else if (key == "device") {
      cfg.device = (lval == "auto") ? -1 : atoi(val.c_str());
    } else if (key == "worktodo_file") {
      cfg.worktodoFile = val;
    } else if (key == "results_file") {
      cfg.resultsFile = val;
    } else if (key == "username" || key == "user") {
      cfg.username = val;
    } else if (key == "computer_name" || key == "computer") {
      cfg.computerName = val;
    } else if (key == "checkpoint") {
      if (!parseBool(val, cfg.checkpoint)) { return bad("must be 0/1"); }
    } else if (key == "checkpoint_file") {
      cfg.checkpointFile = (lval == "auto") ? "" : val;
    } else if (key == "checkpoint_seconds") {
      if (!parseNumber(val, n)) { return bad("must be a number"); }
      cfg.checkpointSeconds = u32(n);
    } else if (key == "gcd_threads" || key == "cpus" || key == "threads") {
      if (lval == "auto") { cfg.gcdThreads = 0; }
      else if (!parseNumber(val, n) || n > 256) { return bad("must be 0..256 or 'auto'"); }
      else { cfg.gcdThreads = u32(n); }
    } else if (key == "extend") {
      if (!parseBool(val, cfg.extend)) { return bad("must be 0/1"); }
    } else if (key == "bounds_source") {
      if      (lval == "auto")       { cfg.boundsSource = BOUNDS_AUTO; }
      else if (lval == "assignment") { cfg.boundsSource = BOUNDS_ASSIGNMENT; }
      else { return bad("must be auto or assignment"); }
    } else if (key == "wait_for_work") {
      if (!parseBool(val, cfg.waitForWork)) { return bad("must be 0/1"); }
    } else if (key == "wait_poll_seconds") {
      if (!parseNumber(val, n) || n == 0) { return bad("must be > 0"); }
      cfg.waitPollSeconds = u32(n);
    } else if (key == "verify_fft") {
      if (!parseBool(val, cfg.verifyFft)) { return bad("must be 0/1"); }
    } else if (key == "pause_on_exit" || key == "pause") {
      if      (lval == "auto")   { cfg.pauseMode = PAUSE_AUTO; }
      else if (lval == "1" || lval == "true"  || lval == "always" || lval == "yes") { cfg.pauseMode = PAUSE_ALWAYS; }
      else if (lval == "0" || lval == "false" || lval == "never"  || lval == "no")  { cfg.pauseMode = PAUSE_NEVER; }
      else { return bad("must be auto, always or never"); }
    } else if (key == "report_every") {
      if (!parseNumber(val, n) || n == 0) { return bad("must be > 0"); }
      cfg.reportEvery = u32(n);
    } else {
      err = "line " + to_string(lineNo) + ": unknown key '" + key + "'";
      fclose(f);
      return false;
    }
  }
  fclose(f);

  // No exponent requirement here -- exponents come from worktodo.txt, read
  // separately by the driver once this config is loaded. factoredTo's
  // auto-default (0 == pick from exponent) is likewise resolved by the
  // driver per queue entry, not here: this file has no exponent to resolve
  // it against, and different entries can need different defaults.

  // Reject settings that contradict each other rather than silently picking one.
  // A b2 that can never be used almost always means the stages line is not what
  // the user thought it was.
  if (cfg.b2 && cfg.stage2Mode == STAGE2_OFF) {
    err = "b2 is set but stages = 1, so stage 2 would never run";
    return false;
  }
  if (cfg.b1 && cfg.b2 && cfg.b2 <= cfg.b1) {
    err = "b2 must be greater than b1";
    return false;
  }
  if (cfg.stage2W && !cfg.stage2D) {
    err = "stage2_w needs stage2_d as well: the window is measured in halves of D";
    return false;
  }
  if (cfg.stage2D && cfg.stage2W && u64(cfg.stage2D) * cfg.stage2W / 2 > (cfg.b1 ? cfg.b1 : ~0ull)) {
    err = "stage2_d * stage2_w / 2 must not exceed b1";
    return false;
  }
  return true;
}
