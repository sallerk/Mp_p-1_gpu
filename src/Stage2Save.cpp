// Copyright (C) Mp_p-1_gpu

#include "Stage2Save.h"

#include <cstdio>
#include <cstring>
#include <filesystem>
#include <string>

using namespace std;

namespace {

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
  FILE* f = fopen(path.c_str(), "rb");
  if (!f) { err = "no stage-2 checkpoint"; return false; }

  Header h{};
  if (fread(&h, sizeof(h), 1, f) != 1) { fclose(f); err = "truncated header"; return false; }
  if (memcmp(h.magic, MAGIC, sizeof(MAGIC)) != 0) { fclose(f); err = "not a stage-2 checkpoint"; return false; }

  // Each of these means the residues came from a different walk of a different
  // plan. Adapting rather than rejecting would silently drop primes.
  auto reject = [&](const string& why) { fclose(f); err = why; return false; };
  if (h.version != want.version) {
    return reject("written by stage-2 format v" + to_string(h.version) +
                  ", this build uses v" + to_string(want.version) + " -- discarding");
  }
  if (h.exponent != want.exponent) { return reject("for a different exponent (M" + to_string(h.exponent) + ")"); }
  if (h.b1 != want.b1) { return reject("for B1=" + to_string(h.b1) + ", this run wants B1=" + to_string(want.b1)); }
  if (h.b2 != want.b2) { return reject("for B2=" + to_string(h.b2) + ", this run wants B2=" + to_string(want.b2)); }
  if (h.d != want.d || h.w != want.w) {
    return reject("for pairing shape D=" + to_string(h.d) + " w=" + to_string(h.w) +
                  ", this run uses D=" + to_string(want.d) + " w=" + to_string(want.w));
  }

  Words acc, a, s;
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
  out.acc = std::move(acc);
  out.a = std::move(a);
  out.s = std::move(s);
  return true;
}
