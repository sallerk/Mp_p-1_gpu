// Copyright (C) Mp_p-1_gpu

#include "Stage2Save.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <string>

using namespace std;

namespace {

// Unchanged from v1. The layout did change, but `version` sits at a fixed offset
// straight after the magic and is checked before anything else is believed, so a
// v1 file is rejected with a message that says what it is rather than with a
// generic "not a stage-2 checkpoint".
const char MAGIC[8] = {'P','M','1','S','2','V','0','1'};

struct Header {
  char magic[8];
  u32 version;
  u32 exponent;
  u64 b1;
  u64 b2;
  u32 d;
  u32 w;
  u64 m;
  u64 jIdx;
  u64 done;
  u32 complete;
  u32 pad;                  // keep the u64s that follow 8-byte aligned
  u64 fromB2;
  u64 xRes64;
  u32 nAcc, crcAcc;
  u32 nA, crcA;
  u32 nS, crcS;
};

bool writeWordsTo(FILE* f, const Words& w) {
  return w.empty() || fwrite(w.data(), sizeof(u32), w.size(), f) == w.size();
}

bool readWordsFrom(FILE* f, u32 n, Words& out) {
  out.assign(n, 0);
  return n == 0 || fread(out.data(), sizeof(u32), n, f) == n;
}

// The half both loaders share: read and validate the file, check nothing about
// what the caller wanted. Returns false with a reason on a missing, truncated,
// foreign or corrupt file.
bool readStage2File(const string& path, Stage2State& out, string& err) {
  FILE* f = fopen(path.c_str(), "rb");
  if (!f) { err = "no stage-2 checkpoint"; return false; }

  Header h{};
  if (fread(&h, sizeof(h), 1, f) != 1) { fclose(f); err = "truncated header"; return false; }
  if (memcmp(h.magic, MAGIC, sizeof(MAGIC)) != 0) { fclose(f); err = "not a stage-2 checkpoint"; return false; }
  if (h.version != STAGE2_FORMAT_VERSION) {
    fclose(f);
    err = "written by stage-2 format v" + to_string(h.version) +
          ", this build uses v" + to_string(STAGE2_FORMAT_VERSION) + " -- discarding";
    return false;
  }

  Words acc, a, s;
  auto reject = [&](const string& why) { fclose(f); err = why; return false; };
  if (!readWordsFrom(f, h.nAcc, acc)) { return reject("truncated accumulator"); }
  if (!readWordsFrom(f, h.nA, a)) { return reject("truncated A"); }
  if (!readWordsFrom(f, h.nS, s)) { return reject("truncated S"); }
  fclose(f);

  if (crc32(acc) != h.crcAcc) { err = "CRC mismatch on the accumulator -- file is corrupt"; return false; }
  if (crc32(a) != h.crcA) { err = "CRC mismatch on A -- file is corrupt"; return false; }
  if (crc32(s) != h.crcS) { err = "CRC mismatch on S -- file is corrupt"; return false; }

  out.exponent = h.exponent;
  out.b1 = h.b1;
  out.b2 = h.b2;
  out.d = h.d;
  out.w = h.w;
  out.version = h.version;
  out.m = h.m;
  out.jIdx = h.jIdx;
  out.done = h.done;
  out.complete = h.complete != 0;
  out.fromB2 = h.fromB2;
  out.xRes64 = h.xRes64;
  out.acc = std::move(acc);
  out.a = std::move(a);
  out.s = std::move(s);
  return true;
}

} // namespace

string defaultStage2Path(u32 exponent, u64 b1, u64 b2) {
  char buf[160];
  snprintf(buf, sizeof(buf), "pm1_%u_s2_%llu_%llu.save", exponent,
           (unsigned long long) b1, (unsigned long long) b2);
  return buf;
}

bool saveStage2(const string& path, const Stage2State& s, string& err) {
  const string tmp = path + ".tmp";
  FILE* f = fopen(tmp.c_str(), "wb");
  if (!f) { err = "cannot open " + tmp; return false; }

  Header h{};
  memcpy(h.magic, MAGIC, sizeof(MAGIC));
  h.version = s.version;
  h.exponent = s.exponent;
  h.b1 = s.b1;
  h.b2 = s.b2;
  h.d = s.d;
  h.w = s.w;
  h.m = s.m;
  h.jIdx = s.jIdx;
  h.done = s.done;
  h.complete = s.complete ? 1u : 0u;
  h.fromB2 = s.fromB2;
  h.xRes64 = s.xRes64;
  h.nAcc = u32(s.acc.size()); h.crcAcc = crc32(s.acc);
  h.nA = u32(s.a.size());     h.crcA = crc32(s.a);
  h.nS = u32(s.s.size());     h.crcS = crc32(s.s);

  bool ok = fwrite(&h, sizeof(h), 1, f) == 1;
  if (ok) { ok = writeWordsTo(f, s.acc); }
  if (ok) { ok = writeWordsTo(f, s.a); }
  if (ok) { ok = writeWordsTo(f, s.s); }
  if (fclose(f) != 0) { ok = false; }

  if (!ok) { err = "write failed"; remove(tmp.c_str()); return false; }

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

bool loadStage2(const string& path, const Stage2State& want, Stage2State& out,
                string& err) {
  Stage2State got;
  if (!readStage2File(path, got, err)) { return false; }

  // Each of these means the residues came from a different walk of a different
  // plan. Adapting rather than rejecting would silently drop primes.
  auto reject = [&](const string& why) { err = why; return false; };
  if (got.exponent != want.exponent) { return reject("for a different exponent (M" + to_string(got.exponent) + ")"); }
  if (got.b1 != want.b1) { return reject("for B1=" + to_string(got.b1) + ", this run wants B1=" + to_string(want.b1)); }
  if (got.b2 != want.b2) { return reject("for B2=" + to_string(got.b2) + ", this run wants B2=" + to_string(want.b2)); }
  if (got.d != want.d || got.w != want.w) {
    return reject("for pairing shape D=" + to_string(got.d) + " w=" + to_string(got.w) +
                  ", this run uses D=" + to_string(want.d) + " w=" + to_string(want.w));
  }
  // Same bounds, same shape, same position -- and yet a different accumulator,
  // because one of the two walks started from a seeded product and the other
  // from 1.
  if (got.fromB2 != want.fromB2) {
    return reject(got.fromB2 ? "seeded from a completed B2=" + to_string(got.fromB2) +
                               ", this run walks the whole range from B1"
                             : "walked the whole range from B1, this run extends a completed B2=" +
                               to_string(want.fromB2));
  }
  if (want.xRes64 && got.xRes64 != want.xRes64) {
    return reject("built from a different stage-1 residue");
  }

  out = std::move(got);
  return true;
}

bool loadCompletedStage2(const string& path, u32 exponent, u64 b1, u64 xRes64,
                         Stage2State& out, string& err) {
  Stage2State got;
  if (!readStage2File(path, got, err)) { return false; }

  if (got.exponent != exponent) { err = "for a different exponent (M" + to_string(got.exponent) + ")"; return false; }
  if (got.b1 != b1) { err = "for B1=" + to_string(got.b1) + ", this run wants B1=" + to_string(b1); return false; }
  if (!got.complete) { err = "stage 2 to B2=" + to_string(got.b2) + " never finished"; return false; }
  // A pre-1.1 run left no xRes64 to compare against; there are no such files,
  // since the format version bump already rejected them, but the check is
  // written so that a zero on either side means "unknown", not "matches".
  if (!got.xRes64 || !xRes64 || got.xRes64 != xRes64) {
    err = "built from a different stage-1 residue";
    return false;
  }

  out = std::move(got);
  return true;
}

std::vector<u64> findStage2Saves(u32 exponent, u64 b1) {
  // Files are named pm1_<exponent>_s2_<b1>_<b2>.save. Largest B2 first, so the
  // most valuable seed is found before the ones it supersedes.
  std::vector<u64> out;
  const string prefix = "pm1_" + to_string(exponent) + "_s2_" + to_string(b1) + "_";
  std::error_code ec;
  for (const auto& e : filesystem::directory_iterator(".", ec)) {
    const string name = e.path().filename().string();
    if (name.rfind(prefix, 0) != 0 || name.size() <= prefix.size() + 5) { continue; }
    if (name.substr(name.size() - 5) != ".save") { continue; }
    const string mid = name.substr(prefix.size(), name.size() - prefix.size() - 5);
    char* end = nullptr;
    const u64 v = strtoull(mid.c_str(), &end, 10);
    if (end && *end == 0 && v) { out.push_back(v); }
  }
  std::sort(out.begin(), out.end(), std::greater<u64>());
  return out;
}
