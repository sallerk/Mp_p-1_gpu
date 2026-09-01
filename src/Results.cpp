// Copyright (C) Mp_p-1_gpu

#include "Results.h"
#include "Config.h"
#include "PM1.h"

#include <cstdio>
#include <cstring>
#include <ctime>
#include "BigInt.h"

#include <algorithm>
#include <string>

// One line per PRIME. A gcd routinely carries several factors multiplied
// together -- every one whose k was smooth comes out of the same gcd -- so
// reporting it raw would be a composite masquerading as a factor.
// One JSON object per line, the shape PrimeNet accepts and the same convention
// Prime95 uses for its own results.txt -- so the file can be uploaded to
// mersenne.org's "Manual Results" page as-is.
//
//   {"status":"F","exponent":81679223,"worktype":"P-1","b1":2000000,
//    "b2":60000000,"factors":["..."],"program":{"name":"...","version":"1.0"},
//    "timestamp":"2026-07-28 17:29:54"}
//
// When the job came from a Pfactor=/Pminus1= worktodo entry (see Worktodo.h),
// "aid" and "known-factors" are added the same way "user"/"computer" already
// are -- present only when the assignment carried them, so AutoPrimeNet's own
// upload step can match this line back to the assignment it came from.
//
// Only factors that are PRIME and verified to divide M_p are reported as
// factors: a gcd routinely carries several multiplied together, and submitting
// that product would be a composite masquerading as a factor. Anything that
// could not be split is written with status "C" and is not a submittable
// result -- it is recorded so the run is not silently lost.
void writeResultJson(const Config& cfg, const char* worktype, u64 b1, u64 b2,
                     const std::vector<FoundFactor>& factors, const Pp1Start* start,
                     u32 stage2D) {
  FILE* f = fopen(cfg.resultsFile.c_str(), "a");
  if (!f) {
    printf("  WARNING: could not append to %s\n", cfg.resultsFile.c_str());
    return;
  }

  char stamp[32] = "";
  const time_t now = time(nullptr);
  struct tm utc;
  // MSVC and POSIX spell the reentrant gmtime differently, and swap the
  // argument order while they are at it.
#ifdef _WIN32
  const bool haveUtc = gmtime_s(&utc, &now) == 0;
#else
  const bool haveUtc = gmtime_r(&now, &utc) != nullptr;
#endif
  if (haveUtc) { strftime(stamp, sizeof(stamp), "%Y-%m-%d %H:%M:%S", &utc); }

  // The fields between "worktype" and the trailing metadata, in Prime95's own
  // order and spacing, so a line from this program and a line from Prime95
  // differ only where they must. "factors" is emitted by the caller between
  // head() and tail(), which is where Prime95 puts it.
  auto head = [&]() {
    fprintf(f, ", \"exponent\":%u, \"worktype\":\"%s\"", cfg.exponent, worktype);
  };
  auto tail = [&]() {
    fprintf(f, ", \"b1\":%llu", (unsigned long long) b1);
    // Omitted when stage 2 never ran, exactly as Prime95 omits it for a
    // stage-1-only result.
    const bool ranStage2 = b2 > b1;
    if (ranStage2) { fprintf(f, ", \"b2\":%llu", (unsigned long long) b2); }
    // P-1 stage 2 only. Prime95 gates "b2" and "d" together -- on C > B for a
    // no-factor result, and on "did stage 2 produce the factor" for a find --
    // and this program follows both rules: b2 > b1 here, and the P-1 callers
    // pass stage2D only when the result came out of stage 2.
    //
    // Everything else Prime95 puts on such a line lives inside its
    //     if (stage2_type == PM1_STAGE2_POLYMULT)
    // branch: "poly1-size", "poly2-size" AND "stage2-fft-length". Those three
    // describe the polynomial multiplication its stage 2 uses from 30.8 on.
    // This program pairs primes in Montgomery's scheme, so none of the three
    // applies -- there are no polynomials to size, and emitting the FFT length
    // under a name Prime95 only ever uses for polymult would imply an
    // algorithm this program does not run. A draft of this release emitted
    // stage2-fft-length on the reasoning that it was at least a true
    // statement; it is, but it is a true statement filed under a key that
    // means "polymult ran".
    //
    // Prime95 emits none of these for P+1 either, so the P+1 callers pass 0.
    // `start` non-null means P+1, and Prime95 writes no "d" for P+1. Every
    // P+1 caller already passes stage2D as 0, but the rule belongs here: a
    // contract with another program should not depend on all six call sites
    // remembering it.
    if (ranStage2 && stage2D && !start) { fprintf(f, ", \"d\":%u", stage2D); }
    // P+1 only: which rational the Lucas sequence started from. P-1 has no
    // equivalent, and Prime95 emits none for it either.
    if (start) { fprintf(f, ", \"start\":\"%s\"", start->label().c_str()); }
    if (cfg.fftLength) { fprintf(f, ", \"fft-length\":%llu", (unsigned long long) cfg.fftLength); }
    // No "security-code": that is Prime95's SEC5 checksum over its own state,
    // and a value computed here would not be one. PrimeNet treats it as
    // optional; a fabricated one would be worse than none. Same for "build"
    // and "port", which describe a Prime95 install.
    fprintf(f, ", \"program\":{\"name\":\"%s\", \"version\":\"%s\"}", PROGRAM_NAME, PROGRAM_VERSION);
    fprintf(f, ", \"timestamp\":\"%s\"", stamp);
    // Prime95's order: known-factors between the timestamp and the user, aid
    // last of all.
    if (!cfg.knownFactors.empty()) {
      fprintf(f, ", \"known-factors\":[");
      for (size_t i = 0; i < cfg.knownFactors.size(); ++i) {
        fprintf(f, "%s\"%s\"", i ? ", " : "", cfg.knownFactors[i].c_str());
      }
      fprintf(f, "]");
    }
    if (!cfg.username.empty()) { fprintf(f, ", \"user\":\"%s\"", cfg.username.c_str()); }
    if (!cfg.computerName.empty()) { fprintf(f, ", \"computer\":\"%s\"", cfg.computerName.c_str()); }
    if (!cfg.aid.empty()) { fprintf(f, ", \"aid\":\"%s\"", cfg.aid.c_str()); }
  };

  std::vector<const FoundFactor*> good, bad;
  for (const FoundFactor& ff : factors) {
    ((ff.prime && ff.dividesMp) ? good : bad).push_back(&ff);
  }

  if (!good.empty()) {
    fprintf(f, "{\"status\":\"F\"");
    head();
    fprintf(f, ", \"factors\":[");
    for (size_t i = 0; i < good.size(); ++i) {
      fprintf(f, "%s\"%s\"", i ? ", " : "", good[i]->value.dec().c_str());
    }
    fprintf(f, "]");
    tail();
    fprintf(f, "}\n");
  }
  for (const FoundFactor* ff : bad) {
    fprintf(f, "{\"status\":\"C\"");
    head();
    fprintf(f, ", \"composite\":\"%s\", \"note\":\"%s\"", ff->value.dec().c_str(),
            ff->dividesMp ? "could not be split into primes"
                          : "DOES NOT DIVIDE M_p -- please report");
    tail();
    fprintf(f, "}\n");
  }
  if (good.empty() && bad.empty()) {
    fprintf(f, "{\"status\":\"NF\"");
    head();
    tail();
    fprintf(f, "}\n");
  }
  fclose(f);
}

// ---------------------------------------------------------------------------
// Self-test: the field set and ORDER of every line shape, against Prime95's.
//
// This exists because the contract is with another program and was got wrong
// twice: "seed" where Prime95 writes "start", and "stage2-fft-length", which
// Prime95 emits only inside its `stage2_type == PM1_STAGE2_POLYMULT` branch
// and which therefore has no business on a line from this program at all.
// Neither was caught by running anything -- both were caught by reading
// ecm.cpp. Section E below is the direct guard against a third.
// ---------------------------------------------------------------------------

namespace {

// The top-level keys of one JSON object, in order. A small hand scanner
// rather than a parser on purpose: what is under test is the literal text
// this program emits, and anything that normalised it would test less.
std::vector<std::string> topLevelKeys(const std::string& line) {
  std::vector<std::string> keys;
  int depth = 0;
  bool inStr = false, esc = false;
  std::string cur;
  for (size_t i = 0; i < line.size(); ++i) {
    const char c = line[i];
    if (inStr) {
      if (esc) { esc = false; }
      else if (c == '\\') { esc = true; }
      else if (c == '"') {
        inStr = false;
        size_t j = i + 1;
        while (j < line.size() && (line[j] == ' ' || line[j] == '\t')) { ++j; }
        // A key is a string at depth 1 immediately followed by ':'. Array
        // elements and the keys inside "program" sit at depth 2 and are skipped.
        if (depth == 1 && j < line.size() && line[j] == ':') { keys.push_back(cur); }
      } else { cur += c; }
      continue;
    }
    if (c == '"') { inStr = true; cur.clear(); }
    else if (c == '{' || c == '[') { ++depth; }
    else if (c == '}' || c == ']') { --depth; }
  }
  return keys;
}

std::string join(const std::vector<std::string>& v) {
  std::string s;
  for (size_t i = 0; i < v.size(); ++i) { s += (i ? ", " : "") + v[i]; }
  return s;
}

} // namespace

int runResultsTests() {
  printf("G8: results.txt line shapes\n\n");
  int checks = 0, failures = 0;
  const char* PATH = "selftest-results.txt";

  Config base;
  base.exponent = 4444091;
  base.resultsFile = PATH;
  base.fftLength = 262144;

  auto factor = [](const char* dec, bool prime, bool divides) {
    FoundFactor ff;
    fromDecimal(dec, ff.value);
    ff.prime = prime;
    ff.dividesMp = divides;
    return ff;
  };
  const std::vector<FoundFactor> found{factor("636358278473", true, true)};
  const std::vector<FoundFactor> composite{factor("55277543419074012358186647", false, true)};

  auto emit = [&](const Config& cfg, const char* wt, u64 b1, u64 b2,
                  const std::vector<FoundFactor>& fs, const Pp1Start* st, u32 d) {
    remove(PATH);
    writeResultJson(cfg, wt, b1, b2, fs, st, d);
    std::string out;
    if (FILE* f = fopen(PATH, "rb")) {
      char buf[8192];
      size_t n;
      while ((n = fread(buf, 1, sizeof(buf), f)) > 0) { out.append(buf, n); }
      fclose(f);
    }
    while (!out.empty() && (out.back() == '\n' || out.back() == '\r')) { out.pop_back(); }
    return out;
  };

  auto shape = [&](const char* what, const std::string& line,
                   const std::vector<std::string>& want) {
    ++checks;
    const std::vector<std::string> got = topLevelKeys(line);
    const bool ok = got == want;
    if (!ok) {
      ++failures;
      printf("  FAIL %s\n       want: %s\n       got : %s\n",
             what, join(want).c_str(), join(got).c_str());
    }
    printf("     %s  %s\n", ok ? "PASS" : "FAIL", what);
  };

  const Pp1Start start27{2, 7};

  printf("  A. P-1, the three shapes Prime95 writes\n");
  shape("NF, stage 1 only          -- no b2, no d",
        emit(base, "P-1", 1000, 1000, {}, nullptr, 0),
        {"status", "exponent", "worktype", "b1", "fft-length", "program", "timestamp"});
  shape("NF, stage 2 ran           -- b2 and d together",
        emit(base, "P-1", 1000, 50000, {}, nullptr, 630),
        {"status", "exponent", "worktype", "b1", "b2", "d", "fft-length", "program", "timestamp"});
  shape("F, found in stage 1       -- factors early, still no b2",
        emit(base, "P-1", 1000, 1000, found, nullptr, 0),
        {"status", "exponent", "worktype", "factors", "b1", "fft-length", "program", "timestamp"});
  shape("F, found in stage 2",
        emit(base, "P-1", 1000, 50000, found, nullptr, 630),
        {"status", "exponent", "worktype", "factors", "b1", "b2", "d", "fft-length", "program", "timestamp"});

  printf("\n  B. P+1: start instead of d, on both statuses\n");
  shape("NF, stage 2 ran           -- start, and NO d",
        emit(base, "P+1", 1000, 50000, {}, &start27, 0),
        {"status", "exponent", "worktype", "b1", "b2", "start", "fft-length", "program", "timestamp"});
  shape("F, stage 1 only",
        emit(base, "P+1", 1000, 1000, found, &start27, 0),
        {"status", "exponent", "worktype", "factors", "b1", "start", "fft-length", "program", "timestamp"});

  printf("\n  C. a divisor that could not be split is status C, not F\n");
  shape("C, composite divisor",
        emit(base, "P-1", 1000, 50000, composite, nullptr, 630),
        {"status", "exponent", "worktype", "composite", "note", "b1", "b2", "d",
         "fft-length", "program", "timestamp"});

  printf("\n  D. the metadata tail, in Prime95's order\n");
  {
    Config cfg = base;
    cfg.username = "Kriesel";
    cfg.computerName = "roc";
    cfg.aid = "88D8BAFFFF12E5DDD3FB093FEFE04025";
    cfg.knownFactors = {"8888183", "319974553"};
    shape("known-factors before user, aid last",
          emit(cfg, "P-1", 1000, 50000, found, nullptr, 630),
          {"status", "exponent", "worktype", "factors", "b1", "b2", "d", "fft-length",
           "program", "timestamp", "known-factors", "user", "computer", "aid"});
  }
  {
    Config cfg = base;
    cfg.fftLength = 0;                 // not known: the field goes away entirely
    shape("fft-length omitted when unknown",
          emit(cfg, "P-1", 1000, 1000, {}, nullptr, 0),
          {"status", "exponent", "worktype", "b1", "program", "timestamp"});
  }

  printf("\n  E. fields Prime95 writes that this program must NOT\n");
  {
    // security-code is Prime95's SEC5 over its own state. poly1-size,
    // poly2-size and stage2-fft-length live inside its polymult branch and
    // describe an algorithm this program does not run. build and port
    // describe a Prime95 install. "seed" was this program's own old spelling
    // of "start". None of them may appear again.
    static const char* FORBIDDEN[] = {"security-code", "poly1-size", "poly2-size",
                                      "stage2-fft-length", "build", "port", "seed"};
    Config cfg = base;
    cfg.username = "u";
    cfg.computerName = "c";
    cfg.aid = "A";
    cfg.knownFactors = {"8888183"};
    const std::string lines[] = {
      emit(cfg, "P-1", 1000, 1000,  {},        nullptr,  0),
      emit(cfg, "P-1", 1000, 50000, {},        nullptr,  630),
      emit(cfg, "P-1", 1000, 50000, found,     nullptr,  630),
      emit(cfg, "P+1", 1000, 50000, {},        &start27, 0),
      emit(cfg, "P+1", 1000, 50000, found,     &start27, 0),
      emit(cfg, "P-1", 1000, 50000, composite, nullptr,  630),
    };
    bool clean = true;
    for (const std::string& line : lines) {
      for (const std::string& k : topLevelKeys(line)) {
        for (const char* bad : FORBIDDEN) {
          if (k == bad) { clean = false; printf("  FAIL emitted %s\n", bad); }
        }
      }
    }
    ++checks;
    if (!clean) { ++failures; }
    printf("     %s  none of security-code, poly1-size, poly2-size,\n"
           "           stage2-fft-length, build, port, seed appears on any line\n",
           clean ? "PASS" : "FAIL");
  }

  printf("\n  F. d belongs to P-1 stage 2, start to P+1, never the reverse\n");
  {
    // Every P+1 caller passes stage2D as 0; prove the writer would still not
    // emit "d" for P+1 if that ever changed, and that P-1 grows no "start".
    const std::vector<std::string> pp1 =
        topLevelKeys(emit(base, "P+1", 1000, 50000, {}, &start27, 630));
    const std::vector<std::string> pm1 =
        topLevelKeys(emit(base, "P-1", 1000, 50000, {}, nullptr, 630));
    auto has = [](const std::vector<std::string>& v, const char* k) {
      return std::find(v.begin(), v.end(), k) != v.end();
    };
    ++checks;
    const bool ok = has(pp1, "start") && !has(pp1, "d")
                 && has(pm1, "d") && !has(pm1, "start");
    if (!ok) { ++failures; }
    printf("     %s  P+1 carries start and no d; P-1 carries d and no start\n",
           ok ? "PASS" : "FAIL");
  }

  remove(PATH);
  printf("\nG8: %d/%d checks passed.\n", checks - failures, checks);
  return failures == 0 ? 0 : 1;
}
