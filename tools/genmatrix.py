#!/usr/bin/env python3
"""Generate a worktodo.txt acceptance matrix, and the table of what each line
should produce.

Usage, from a scratch directory holding the binary and a config.txt:

    python tools/genmatrix.py accept    > writes worktodo.txt + expected.json
    <run the program>
    python tools/checkmatrix.py accept  > compares results.txt and the log

    python tools/genmatrix.py resume    > the checkpoint-reuse scenario
    <run the program>
    python tools/checkmatrix.py resume

The "accept" scenario needs config.txt to pin:

    method = pm1   stages = both   b1 = 200   b2 = 2000
    factored_to = auto   bias = 2.0   checkpoint = 0

checkpoint = 0 matters: it makes each line independent, so a line cannot be
served by the residue another line left behind. The "resume" scenario needs
checkpoint = 1 for exactly the opposite reason.

WHY M1000273
------------
It is a low 7-digit prime with three factors whose smoothness makes the outcome
steerable rather than incidental:

    F1 = 16004369      k = 8     = 2^3       needs B1 >= 8
    F2 = 390106471     k = 195   = 3*5*13    needs B1 >= 13
    F3 = 13747752113   k = 6872  = 2^3*859   needs B1 >= 859, or 859 in reach
                                             of stage 2's grid

Two rules that are easy to get wrong, both measured rather than assumed:

  * Stage 1's exponent carries q^floor(log_q B1), so a factor whose k is 2^3
    needs B1 >= 8, NOT B1 >= 2. Measured: B1=7 misses F1, B1=8 finds it.
  * Stage 2 pairs over m*D +- j, so its grid reaches PAST B2, to
    mLast*D + jmax. A prime above B2 landing on a covered slot is found even
    though it was never targeted. Measured with D=210: B2=858 reaches 945 and
    does find F3 (missing prime 859); B2=525 reaches 735 and does not.

Declaring the easy factors as known_factors is what lets a line reach stage 2:
they still fall out of the stage-1 gcd, but are dropped from the report, so the
job does not end there.
"""
import json
import sys

E   = 1000273
F1, F2, F3 = "16004369", "390106471", "13747752113"
K12  = F1 + "," + F2
KALL = K12 + "," + F3
AID  = "88D8BAFFFF12E5DDD3FB093FEFE04025"
CFG_B1, CFG_B2 = 200, 2000          # what config.txt pins, for Pfactor/bare lines
CFG_FT   = 75                       # config.txt says factored_to = auto, which is 75
CFG_BIAS = 2.0                      # and bias = 2.0

KPOWERS = {F1: [8], F2: [3, 5, 13], F3: [8, 859]}

def grid_reach(b2, d=210, jmax=105):
    return (b2 + jmax) // d * d + jmax

def pm1_expect(b1, b2, known):
    stage1 = [f for f, powers in KPOWERS.items() if all(q <= b1 for q in powers)]
    stage2 = []
    ran2 = b2 > b1 and not [f for f in stage1 if f not in known]
    if ran2:
        for f, powers in KPOWERS.items():
            if f in stage1:
                continue
            missing = [q for q in powers if q > b1]
            if len(missing) == 1 and missing[0] <= grid_reach(b2):
                stage2.append(f)
    return [f for f in (stage1 + stage2) if f not in known], ran2

def split_fields(line):
    """The line's fields after the keyword and any AID, split at top level so a
    comma inside the known-factors quote does not split."""
    if "=" not in line:
        return []
    val, out, cur, inq = line.split("=", 1)[1], [], "", False
    for c in val:
        if c == '"':      inq = not inq; cur += c
        elif c == "," and not inq: out.append(cur); cur = ""
        else:             cur += c
    out.append(cur)
    # a 32-hex-char first field followed by more is an AID, not k
    if len(out) > 1 and len(out[0]) == 32 and all(c in "0123456789abcdefABCDEF" for c in out[0]):
        out = out[1:]
    return out

def _slot(fields, i):
    """A positional numeric slot: its value if it parses, else None."""
    if i >= len(fields):
        return None
    try:
        return float(fields[i])
    except ValueError:
        return None

def derive_ft(line):
    """What factored_to should resolve to. Pfactor='s how_far_factored is a
    required field validated 1..127; Pminus1=/Pplus1='s sieve_depth is a
    positional SLOT that silently drops anything outside 1..127, falling back
    to config.txt. See parsePfactor and parsePminus1 in Worktodo.cpp."""
    f = split_fields(line)
    low = line.split("=", 1)[0].strip().lower() if "=" in line else ""
    idx = {"pfactor": 4, "pminus1": 6, "pplus1": 7}.get(low)
    if idx is None:
        return CFG_FT                       # a bare exponent
    v = _slot(f, idx)
    if v is None or not (1 <= v <= 127):
        return CFG_FT if low != "pfactor" else CFG_FT
    return int(v)                           # floored, as u32(floor(x)) does

def derive_bias(line):
    """tests_saved, which only Pfactor= carries. 0 means "P-1 already done"
    and is treated as absent, so config.txt's bias stands."""
    f = split_fields(line)
    low = line.split("=", 1)[0].strip().lower() if "=" in line else ""
    if low != "pfactor":
        return CFG_BIAS
    v = _slot(f, 5)
    return CFG_BIAS if v is None or v <= 0 else v

CASES = []
def case(line, note, worktype="P-1", b1=None, b2=0, known=(), aid=False, predict=True,
         warn_b2start=0, ft=None, bias=None):
    """warn_b2start: the B2_start this line should be seen to carry, 0 for none.

    This program cannot honour a B2_start -- there is no local checkpoint to
    seed externally-declared progress from -- so it warns and walks the whole
    of (B1, B2] instead. Correct, just more work than PrimeNet expected. What
    is under test is that the warning fires exactly when the field was read.
    """
    known = list(known)
    found, ran2 = pm1_expect(b1, b2, known) if (predict and worktype == "P-1") \
                  else (None, b2 > b1)
    # ft / bias default to what the LINE carries, derived by mirroring the
    # program's own positional rules, so a case does not have to restate what
    # it already says. Pass them explicitly only to assert something the
    # derivation cannot know.
    CASES.append(dict(line=line, note=note, worktype=worktype, b1=b1, b2=b2,
                      known=known, aid=aid, factors=found, ran_stage2=ran2,
                      predict=predict, warn_b2start=warn_b2start,
                      ft=derive_ft(line) if ft is None else ft,
                      bias=derive_bias(line) if bias is None else bias))

def build_accept():
    # ---- Pfactor: no bounds of its own; tests_saved drives the bias ----------
    case(f"Pfactor=1,2,{E},-1,70,2",   "tests_saved 2 -> bias 2 (first-time test)", b1=CFG_B1, b2=CFG_B2, ft=70, bias=2.0)
    case(f"Pfactor=1,2,{E},-1,70,1",   "tests_saved 1 -> bias 1 (double-check)",    b1=CFG_B1, b2=CFG_B2, ft=70, bias=1.0)
    case(f"Pfactor=1,2,{E},-1,70,1.4", "fractional tests_saved",                    b1=CFG_B1, b2=CFG_B2, ft=70, bias=1.4)
    case(f"Pfactor=1,2,{E},-1,70,0.5", "tests_saved below 1",                       b1=CFG_B1, b2=CFG_B2, ft=70, bias=0.5)
    case(f"Pfactor=1,2,{E},-1,70,0",   "tests_saved 0 = P-1 already done -> config bias stands", b1=CFG_B1, b2=CFG_B2, ft=70, bias=CFG_BIAS)
    case(f"Pfactor=1,2,{E},-1,1,2",    "how_far_factored at its minimum, 1",        b1=CFG_B1, b2=CFG_B2, ft=1)
    case(f"Pfactor=1,2,{E},-1,127,2",  "how_far_factored at its maximum, 127",      b1=CFG_B1, b2=CFG_B2, ft=127)
    case(f"Pfactor=1,2,{E},-1,74.5,2", "fractional how_far_factored, floors to 74", b1=CFG_B1, b2=CFG_B2, ft=74)
    case(f"Pfactor=1,2,{E},-1,70,100", "tests_saved at its maximum, 100",           b1=CFG_B1, b2=CFG_B2, ft=70, bias=100.0)
    case(f"Pfactor={AID},1,2,{E},-1,70,2", "AID present",                           b1=CFG_B1, b2=CFG_B2, aid=True)
    case(f'Pfactor=1,2,{E},-1,70,2,"{KALL}"', "every reachable factor already known -> NF", b1=CFG_B1, b2=CFG_B2, known=KALL.split(","))
    case(f'Pfactor={AID},1,2,{E},-1,70,2,"{F1}"', "AID + one known factor",         b1=CFG_B1, b2=CFG_B2, known=[F1], aid=True)

    # ---- Pminus1: its own bounds, overriding config --------------------------
    case(f"Pminus1=1,2,{E},-1,200,2000",  "plain two-stage line", b1=200, b2=2000)
    case(f"Pminus1=1,2,{E},-1,200,0",     "B2 = 0 -> stage 1 alone", b1=200, b2=0)
    case(f"Pminus1=1,2,{E},-1,200,200",   "B2 == B1 -> stage 1 alone", b1=200, b2=200)
    case(f"Pminus1=1,2,{E},-1,200,100",   "B2 < B1 -> clamped, stage 1 alone", b1=200, b2=100)
    case(f"Pminus1=1,2,{E},-1,105,2000",  "B1 exactly at the stage-2 floor", b1=105, b2=2000)
    case(f"Pminus1=1,2,{E},-1,2,0",       "B1 = 2: below 2^3, so not even F1", b1=2, b2=0)
    case(f"Pminus1=1,2,{E},-1,7,0",       "B1 = 7: one short of 2^3", b1=7, b2=0)
    case(f"Pminus1=1,2,{E},-1,8,0",       "B1 = 8: exactly 2^3, finds F1", b1=8, b2=0)
    case(f"Pminus1=1,2,{E},-1,12,0",      "B1 = 12: one short of reaching F2", b1=12, b2=0)
    case(f"Pminus1=1,2,{E},-1,13,0",      "B1 = 13: exactly reaches F2", b1=13, b2=0)
    case(f'Pminus1=1,2,{E},-1,200,2000,70,0,"{K12}"', "easy factors known -> stage 2 finds F3", b1=200, b2=2000, known=K12.split(","))
    case(f'Pminus1=1,2,{E},-1,200,525,70,0,"{K12}"',  "B2 whose pairing grid cannot reach 859 -> NF", b1=200, b2=525, known=K12.split(","))
    case(f'Pminus1=1,2,{E},-1,200,858,70,0,"{K12}"',  "B2 = 858 < 859, but the grid reaches 945 -> F", b1=200, b2=858, known=K12.split(","))
    case(f'Pminus1=1,2,{E},-1,200,859,70,0,"{K12}"',  "B2 = 859 exactly", b1=200, b2=859, known=K12.split(","))
    case(f'Pminus1=1,2,{E},-1,859,0,70,0,"{K12}"',    "B1 = 859 -> stage 1 finds F3, no stage 2", b1=859, b2=0, known=K12.split(","))
    case(f'Pminus1=1,2,{E},-1,858,0,70,0,"{K12}"',    "B1 = 858 -> misses F3 entirely -> NF", b1=858, b2=0, known=K12.split(","))
    case(f'Pminus1={AID},1,2,{E},-1,200,2000,70,0,"{K12}"', "AID + a stage-2 find", b1=200, b2=2000, known=K12.split(","), aid=True)
    case(f"Pminus1=1,2,{E},-1,200,2000,70",     "sieve_depth only, no B2_start or factors", b1=200, b2=2000, ft=70)
    case(f"Pminus1=1,2,{E},-1,200,2000,127",    "sieve_depth at its maximum", b1=200, b2=2000, ft=127)
    case(f"Pminus1=1,2,{E},-1,200,2000,128",    "sieve_depth one over: discarded, falls back to config", b1=200, b2=2000, ft=CFG_FT)
    case(f"Pminus1=1,2,{E},-1,200,2000,400",    "400 in the sieve_depth slot: discarded", b1=200, b2=2000, ft=CFG_FT)
    # The sieve_depth SLOT takes 1..127 and silently drops anything else, which
    # is how a positional field distinguishes a TF depth from a B2_start. Each
    # of these must fall back to config.txt's value rather than land as one.
    case(f"Pminus1=1,2,{E},-1,200,2000,1",      "sieve_depth at its minimum, 1", b1=200, b2=2000, ft=1)
    case(f"Pminus1=1,2,{E},-1,200,2000,0",      "sieve_depth 0: below the slot range, discarded", b1=200, b2=2000, ft=CFG_FT)
    case(f"Pminus1=1,2,{E},-1,200,2000,0.5",    "sieve_depth 0.5: below 1, discarded", b1=200, b2=2000, ft=CFG_FT)
    case(f"Pminus1=1,2,{E},-1,200,2000,-1",     "sieve_depth negative, discarded", b1=200, b2=2000, ft=CFG_FT)
    case(f"Pminus1=1,2,{E},-1,200,2000,74.9",   "fractional sieve_depth floors to 74", b1=200, b2=2000, ft=74)
    # Prime95 keeps anything <= 300 in this slot; this program keeps only 1..127,
    # the range a bit depth can actually take. 300 is the case where the two
    # differ, and it is the safer direction: a 300-bit TF depth is not a depth.
    case(f"Pminus1=1,2,{E},-1,200,2000,300",    "sieve_depth 300: Prime95 keeps it, this does not", b1=200, b2=2000, ft=CFG_FT)
    # A trailing comma leaves an empty sieve_depth slot, which is consumed and
    # ignored. Prime95 does the same (atof("") is 0, which passes its <= 300
    # test), so this is accepted rather than refused.
    case(f"Pminus1=1,2,{E},-1,200,2000,",      "trailing comma: empty slot, consumed and ignored", b1=200, b2=2000)
    # B2_start: PrimeNet saying "(B1, B2_start] was already walked elsewhere".
    # Not the same thing as this program's own checkpoint-based B2 extension,
    # which needs a local completed record -- see the "resume" scenario.
    case(f"Pminus1=1,2,{E},-1,200,2000,70,1500","B2_start = 1500 > B1: read, warned about", b1=200, b2=2000, warn_b2start=1500)
    case(f"Pminus1=1,2,{E},-1,200,2000,70,2000","B2_start == B2: read, warned about", b1=200, b2=2000, warn_b2start=2000)
    case(f"Pminus1=1,2,{E},-1,200,2000,70,99999","B2_start beyond B2: read, warned about", b1=200, b2=2000, warn_b2start=99999)
    case(f"Pminus1=1,2,{E},-1,200,2000,70,150", "B2_start = 150 < B1: discarded, no warning", b1=200, b2=2000)
    case(f"Pminus1=1,2,{E},-1,200,2000,70,200", "B2_start == B1 exactly: discarded, needs strictly >", b1=200, b2=2000)
    # The round-trip hole, and it is upstream's: Prime95's writer omits
    # sieve_depth when it is 0, which puts B2_start in the sieve_depth SLOT,
    # where both readers discard it as out of range and never look further. The
    # value is lost by Prime95 too. Asserted here so that if it is ever fixed,
    # it is fixed deliberately.
    case(f"Pminus1=1,2,{E},-1,200,2000,8000",   "B2_start with sieve_depth omitted: silently lost", b1=200, b2=2000)
    case(f'Pminus1=1,2,{E},-1,1000,20000,70,0,"{KALL}"', "everything known, both stages run -> NF", b1=1000, b2=20000, known=KALL.split(","))
    case(f'Pminus1=1,2,{E},-1,200,2000,127,1500,"{K12}"', "max sieve_depth + B2_start + known + stage-2 find", b1=200, b2=2000, known=K12.split(","), warn_b2start=1500)
    case(f"Pminus1=1,2,{E},-1,1000,1001",  "B2 exactly one above B1", b1=1000, b2=1001)
    case(f"Pminus1=1,2,{E},-1,105,106",    "floor B1 with the narrowest possible stage 2", b1=105, b2=106)
    case(f'Pminus1=1,2,{E},-1,200,2000,70,0,"{F1},{F1}"', "the same known factor twice", b1=200, b2=2000, known=[F1, F1])
    case(f'Pminus1=1,2,{E},-1,200,2000,70,0,"00{F1}"', "known factor with leading zeros", b1=200, b2=2000, known=["00" + F1])
    case(f"PMINUS1=1,2,{E},-1,200,0",      "keyword in upper case", b1=200, b2=0)
    case(f"pminus1=1,2,{E},-1,200,0",      "keyword in lower case", b1=200, b2=0)

    # ---- Pplus1: a different method, so factors are not predicted ------------
    for nth, note in ((1, "start 2/7"), (2, "start 6/5"), (3, "a derived random start")):
        case(f"Pplus1=1,2,{E},-1,200,2000,{nth}", f"nth_run {nth} -> {note}",
             worktype="P+1", b1=200, b2=2000, predict=False)
    case(f"Pplus1=1,2,{E},-1,200,0,1",    "P+1 with B2 = 0 -> stage 1 alone", worktype="P+1", b1=200, b2=0, predict=False)
    case(f"Pplus1=1,2,{E},-1,200,200,2",  "P+1 with B2 == B1", worktype="P+1", b1=200, b2=200, predict=False)
    case(f"Pplus1={AID},1,2,{E},-1,200,2000,1,70", "P+1 with AID and sieve_depth", worktype="P+1", b1=200, b2=2000, aid=True, predict=False)
    case(f'Pplus1=1,2,{E},-1,200,2000,1,70,"{F1}"', "P+1 with a known factor", worktype="P+1", b1=200, b2=2000, known=[F1], predict=False)
    case(f"Pplus1=1,2,{E},-1,200,2000,2,1",   "P+1 sieve_depth at its minimum", worktype="P+1", b1=200, b2=2000, predict=False, ft=1)
    case(f"Pplus1=1,2,{E},-1,200,2000,3,127", "P+1 sieve_depth at its maximum", worktype="P+1", b1=200, b2=2000, predict=False, ft=127)
    case(f"Pplus1=1,2,{E},-1,200,2000,1,128", "P+1 sieve_depth over range, discarded", worktype="P+1", b1=200, b2=2000, predict=False, ft=CFG_FT)

    # ---- shapes that are not assignment lines at all -------------------------
    case(f"Pminus1=1,2,{E},-1,200,2000", "no depth on the line: config's value applies", b1=200, b2=2000, ft=CFG_FT)
    case(f"{E}", "bare exponent: config's method, bounds and everything else",
         worktype="both", b1=CFG_B1, b2=CFG_B2, predict=False)
    case(f"Pminus1=1,2,{E},-1,13,0", "the same exponent again: consume must take the right line",
         b1=13, b2=0)

def build_resume():
    """Three lines at one (exponent, B1) with rising B2, run WITH checkpointing.

    This is the program's OWN B2 extension -- reusing a completed stage-2
    record on disk -- and has nothing to do with a Pminus1= line's B2_start,
    which is PrimeNet declaring progress made somewhere else entirely and
    cannot be honoured. Nothing in the accept matrix can reach this path,
    because that scenario sets checkpoint = 0 to keep its lines independent,
    and this path is precisely one line being served by another's leftovers.

    Line 2 should walk only (525, 2000] rather than redoing (200, 2000], and
    line 3 should be recognised as already complete and walk nothing.
    """
    case(f'Pminus1=1,2,{E},-1,200,525,70,0,"{K12}"',  "first: completes stage 2 over (200, 525]", b1=200, b2=525, known=K12.split(","))
    case(f'Pminus1=1,2,{E},-1,200,2000,70,0,"{K12}"', "second: should EXTEND from 525 rather than redo (200, 2000]", b1=200, b2=2000, known=K12.split(","))
    case(f'Pminus1=1,2,{E},-1,200,2000,70,0,"{K12}"', "third: identical bounds, should be recognised as already complete", b1=200, b2=2000, known=K12.split(","))

# Lines that must be REFUSED. The queue is fail-fast -- one bad line stops the
# whole run -- so these cannot share a file with the accept cases and are fed
# in one at a time by checkmatrix.py.
REJECT = [
    ("exponent not prime",            f"Pminus1=1,2,1000272,-1,200,2000"),
    ("exponent 0",                    f"Pminus1=1,2,0,-1,200,2000"),
    ("exponent 1",                    f"Pminus1=1,2,1,-1,200,2000"),
    ("exponent 2 (prime but tiny)",   f"Pminus1=1,2,2,-1,200,2000"),
    ("exponent above 2^32-1",         f"Pminus1=1,2,4294967296,-1,200,2000"),
    ("exponent in e-notation",        f"Pminus1=1,2,1e6,-1,200,2000"),
    ("exponent in hex",               f"Pminus1=1,2,0x10,-1,200,2000"),
    ("exponent with a fraction",      f"Pminus1=1,2,1000273.5,-1,200,2000"),
    ("exponent nan",                  f"Pminus1=1,2,nan,-1,200,2000"),
    ("k is not 1",                    f"Pminus1=2,2,{E},-1,200,2000"),
    ("b is not 2",                    f"Pminus1=1,3,{E},-1,200,2000"),
    ("c is not -1",                   f"Pminus1=1,2,{E},1,200,2000"),
    ("b rounds toward 2",             f"Pminus1=1,2.4,{E},-1,200,2000"),
    ("B2 field omitted (Pminus1)",    f"Pminus1=1,2,{E},-1,200"),
    ("B2 field omitted (Pplus1)",     f"Pplus1=1,2,{E},-1,200"),
    ("nth_run missing",               f"Pplus1=1,2,{E},-1,200,2000"),
    ("nth_run zero",                  f"Pplus1=1,2,{E},-1,200,2000,0"),
    ("nth_run negative",              f"Pplus1=1,2,{E},-1,200,2000,-1"),
    ("nth_run not an integer",        f"Pplus1=1,2,{E},-1,200,2000,1.5"),
    ("nth_run 4: only 1, 2, 3 exist", f"Pplus1=1,2,{E},-1,200,2000,4"),
    ("nth_run 1000",                  f"Pplus1=1,2,{E},-1,200,2000,1000"),
    ("B1 below 2",                    f"Pminus1=1,2,{E},-1,1,0"),
    ("B1 zero",                       f"Pminus1=1,2,{E},-1,0,0"),
    ("B1 negative",                   f"Pminus1=1,2,{E},-1,-5,0"),
    ("B1 nan",                        f"Pminus1=1,2,{E},-1,nan,0"),
    ("B1 under the floor WITH a stage 2", f"Pminus1=1,2,{E},-1,104,2000"),
    ("B2 nan",                        f"Pminus1=1,2,{E},-1,200,nan"),
    ("B2 negative",                   f"Pminus1=1,2,{E},-1,200,-5"),
    ("B2 over the machine-word cap",  f"Pminus1=1,2,{E},-1,200,1e18"),
    ("B2 just over the cap",          f"Pminus1=1,2,{E},-1,200,5000000000"),
    ("tests_saved above 100",         f"Pfactor=1,2,{E},-1,70,101"),
    ("tests_saved negative",          f"Pfactor=1,2,{E},-1,70,-1"),
    ("tests_saved nan",               f"Pfactor=1,2,{E},-1,70,nan"),
    ("tests_saved missing",           f"Pfactor=1,2,{E},-1,70"),
    ("how_far_factored 0",            f"Pfactor=1,2,{E},-1,0,2"),
    ("how_far_factored above 127",    f"Pfactor=1,2,{E},-1,128,2"),
    ("how_far_factored nan",          f"Pfactor=1,2,{E},-1,nan,2"),
    ("how_far_factored negative",     f"Pfactor=1,2,{E},-1,-5,2"),
    ("how_far_factored 127.9",        f"Pfactor=1,2,{E},-1,127.9,2"),
    ("how_far_factored 1000",         f"Pfactor=1,2,{E},-1,1000,2"),
    ("how_far_factored inf",          f"Pfactor=1,2,{E},-1,inf,2"),
    ("tests_saved 100.1",             f"Pfactor=1,2,{E},-1,70,100.1"),
    ("tests_saved inf",               f"Pfactor=1,2,{E},-1,70,inf"),
    ("B2 one over the cap",           f"Pminus1=1,2,{E},-1,200,4000000001"),
    ("B1 inf",                        f"Pminus1=1,2,{E},-1,inf,2000"),
    ("B2 inf",                        f"Pminus1=1,2,{E},-1,200,inf"),
    ("exponent 4294967293 not prime", f"Pminus1=1,2,4294967293,-1,200,2000"),
    ("known factor not of this M_p",  f'Pminus1=1,2,{E},-1,200,2000,70,0,"12345"'),
    ("known factor 0",                f'Pminus1=1,2,{E},-1,200,2000,70,0,"0"'),
    ("known factor 1",                f'Pminus1=1,2,{E},-1,200,2000,70,0,"1"'),
    ("known factor not decimal",      f'Pminus1=1,2,{E},-1,200,2000,70,0,"abc"'),
    ("known factors with an empty entry", f'Pminus1=1,2,{E},-1,200,2000,70,0,"{F1},,{F2}"'),
    ("unterminated known-factors quote",  f'Pminus1=1,2,{E},-1,200,2000,70,0,"{F1}'),
    ("empty field",                   f"Pminus1=1,2,{E},-1,,2000"),
    ("spaces around fields",          f"Pminus1= 1 , 2 , {E} , -1 , 200 , 2000"),
    ("a field in quotes",             f'Pminus1=1,2,{E},-1,"200",2000'),
    ("extra field after known factors", f'Pminus1=1,2,{E},-1,200,2000,70,0,"{F1}",9'),
    ("too few fields entirely",       f"Pminus1=1,2,{E}"),
    ("unsupported keyword",           f"PRP=1,2,{E},-1"),
    ("AID of 33 hex chars",           f"Pminus1={AID}A,1,2,{E},-1,200,2000"),
    ("AID with a non-hex character",  f"Pminus1={AID[:-1]}Z,1,2,{E},-1,200,2000"),
]

# Valid lines that are too expensive to RUN -- huge bounds, the largest allowed
# exponent -- but whose parsing is worth asserting. Fed through --bounds, which
# parses and chooses without doing any of the work.
PARSE_OK = [
    ("B2 exactly at the cap",         f"Pminus1=1,2,{E},-1,200,4000000000"),
    ("B1 just under B2 at the cap",   f"Pminus1=1,2,{E},-1,3999999999,4000000000"),
    ("largest allowed exponent",      "Pminus1=1,2,4294967291,-1,200,2000"),
    ("smallest allowed exponent",     "Pminus1=1,2,786613,-1,200,2000"),
    ("AID in lower-case hex",         f"Pminus1={AID.lower()},1,2,{E},-1,200,2000"),
    ("how_far_factored 127 with real bounds", f"Pfactor=1,2,{E},-1,127,2"),
    ("tests_saved 100 with real bounds",      f"Pfactor=1,2,{E},-1,70,100"),
    ("B1 = 2 with no stage 2",        f"Pminus1=1,2,{E},-1,2,0"),
]

def main():
    scenario = sys.argv[1] if len(sys.argv) > 1 else "accept"
    if scenario == "accept":
        build_accept()
    elif scenario == "resume":
        build_resume()
    else:
        raise SystemExit(f"unknown scenario {scenario!r}: use accept or resume")
    with open("worktodo.txt", "w", newline="") as f:
        f.write(f"# {scenario} scenario, generated by tools/genmatrix.py\r\n")
        for c in CASES:
            f.write(c["line"] + "\r\n")
    with open("expected.json", "w") as f:
        json.dump(dict(scenario=scenario, cases=CASES, reject=REJECT,
                   parse_ok=PARSE_OK), f, indent=1)
    print(f"{len(CASES)} lines written for the {scenario} scenario")

if __name__ == "__main__":
    main()
