// Copyright (C) Mp_p-1_gpu

#include "Save.h"

#include <cstdio>
#include <cstring>
#include <filesystem>
#include <string>

using namespace std;

namespace {

// Bumped from ...E1 when the second residue was added: the header layout
// changed, so an old file must be rejected, not reinterpreted.
const char MAGIC[8] = {'P','M','1','S','A','V','E','2'};

struct Header {
  char magic[8];
  u32 eVersion;
  u32 exponent;
  u64 b1;
  u64 eBits;
  u64 nextBit;
  u64 baseB1;
  u32 base;
  u32 complete;
  u32 nWords;
  u32 crc;              // crc32 of the residue words
  u32 nWords2;          // P+1 second residue; 0 for P-1
  u32 crc2;
  u32 seed;
};

} // namespace

string defaultSavePath(u32 exponent, u64 b1) {
  char buf[128];
  snprintf(buf, sizeof(buf), "pm1_%u_b1_%llu.save", exponent, (unsigned long long) b1);
  return buf;
}

bool saveState(const string& path, const SaveState& s, string& err) {
  // Write to a temp file and rename. A checkpoint interrupted mid-write must
  // not replace a good one with a truncated file.
  const string tmp = path + ".tmp";

  FILE* f = fopen(tmp.c_str(), "wb");
  if (!f) { err = "cannot open " + tmp; return false; }

  Header h{};
  memcpy(h.magic, MAGIC, sizeof(MAGIC));
  h.eVersion = s.eVersion;
  h.exponent = s.exponent;
  h.b1 = s.b1;
  h.eBits = s.eBits;
  h.nextBit = s.nextBit;
  h.baseB1 = s.baseB1;
  h.base = s.base;
  h.complete = s.complete ? 1u : 0u;
  h.nWords = u32(s.residue.size());
  h.crc = crc32(s.residue);
  h.nWords2 = u32(s.residue2.size());
  h.crc2 = s.residue2.empty() ? 0 : crc32(s.residue2);
  h.seed = s.seed;

  bool ok = fwrite(&h, sizeof(h), 1, f) == 1;
  if (ok && h.nWords) {
    ok = fwrite(s.residue.data(), sizeof(u32), h.nWords, f) == h.nWords;
  }
  if (ok && h.nWords2) {
    ok = fwrite(s.residue2.data(), sizeof(u32), h.nWords2, f) == h.nWords2;
  }
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

bool loadState(const string& path, const SaveState& want, SaveState& out,
               string& err) {
  FILE* f = fopen(path.c_str(), "rb");
  if (!f) { err = "no checkpoint"; return false; }

  Header h{};
  if (fread(&h, sizeof(h), 1, f) != 1) { fclose(f); err = "truncated header"; return false; }
  if (memcmp(h.magic, MAGIC, sizeof(MAGIC)) != 0) { fclose(f); err = "not a checkpoint file"; return false; }

  // Every mismatch below means the residue describes a DIFFERENT computation.
  // Reusing it would not crash -- it would quietly compute the wrong thing --
  // so each is a hard reject.
  if (h.eVersion != want.eVersion) {
    fclose(f);
    err = "written by exponent format v" + to_string(h.eVersion) +
          ", this build uses v" + to_string(want.eVersion) + " -- discarding";
    return false;
  }
  if (h.exponent != want.exponent) { fclose(f); err = "for a different exponent (M" + to_string(h.exponent) + ")"; return false; }
  if (h.base != want.base) { fclose(f); err = "different base"; return false; }
  if (h.seed != want.seed) { fclose(f); err = "different P+1 seed"; return false; }
  if (h.b1 != want.b1) {
    fclose(f);
    err = "for B1=" + to_string(h.b1) + ", this run wants B1=" + to_string(want.b1);
    return false;
  }
  if (h.eBits != want.eBits) { fclose(f); err = "exponent bit length differs -- E rebuilt differently"; return false; }

  Words residue(h.nWords);
  if (h.nWords && fread(residue.data(), sizeof(u32), h.nWords, f) != h.nWords) {
    fclose(f);
    err = "truncated residue";
    return false;
  }
  Words residue2(h.nWords2);
  if (h.nWords2 && fread(residue2.data(), sizeof(u32), h.nWords2, f) != h.nWords2) {
    fclose(f);
    err = "truncated second residue";
    return false;
  }
  fclose(f);

  if (crc32(residue) != h.crc) { err = "CRC mismatch -- file is corrupt"; return false; }
  if (h.nWords2 && crc32(residue2) != h.crc2) { err = "CRC mismatch on second residue"; return false; }

  out.exponent = h.exponent;
  out.b1 = h.b1;
  out.eBits = h.eBits;
  out.nextBit = h.nextBit;
  out.baseB1 = h.baseB1;
  out.complete = h.complete != 0;
  out.base = h.base;
  out.eVersion = h.eVersion;
  out.residue = std::move(residue);
  out.residue2 = std::move(residue2);
  out.seed = h.seed;
  return true;
}
