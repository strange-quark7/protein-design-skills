import os, sys, json, subprocess, tempfile

# Standard SELEX library layout: 20-nt forward primer + N-core + 20-nt reverse-primer site.
# Swap these for your own primer set.
BASES = "ACGT"
DEFAULT_FWD = "AGCAGCACAGAGGTCAGATG"
DEFAULT_REV = "CCTATGCGTGCTACCGTGAA"

# Standalone folding/characterization script. Run in a subprocess so its worker
# function lives at module scope (kernel-defined functions cannot be pickled for
# multiprocessing.Pool). Reads cores.json -> writes folded.json.
FOLD_SCRIPT = """import json, re, sys
from multiprocessing import Pool
from seqfold import dg
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

d = json.load(open(sys.argv[1]))
fwd, rev, cores, temp = d['fwd'], d['rev'], d['cores'], d['temp']
full = [fwd + c + rev for c in cores]
gq = re.compile(r'G{3,}\\w{1,7}G{3,}\\w{1,7}G{3,}\\w{1,7}G{3,}')

def analyze(args):
    i, s, core = args
    try: g = round(float(dg(s, temp=temp)), 2)
    except Exception: g = float('nan')
    try: tmv = round(float(mt.Tm_NN(Seq(s), nn_table=mt.DNA_NN4)), 1)
    except Exception: tmv = float('nan')
    gc = round(100*(s.count('G')+s.count('C'))/len(s), 1)
    gcc = round(100*(core.count('G')+core.count('C'))/len(core), 1)
    return i, g, tmv, gc, gcc, (1 if gq.search(core) else 0)

if __name__ == '__main__':
    args = list(zip(range(len(full)), full, cores))
    with Pool(d['workers']) as p:
        res = p.map(analyze, args, chunksize=64)
    res.sort()
    json.dump([{'idx':i,'dG_kcal_mol':g,'Tm_C':t,'GC_full_pct':gc,
                'GC_core_pct':gcc,'Gquad_motif':q} for i,g,t,gc,gcc,q in res],
              open(sys.argv[2], 'w'))
"""


def generate_candidate_library(n, core_len=40, fwd=None, rev=None,
                               gc_range=(0.40, 0.60), max_homopolymer=4,
                               max_base_frac=0.50, require_all_bases=True,
                               seed=2024):
    """Generate `n` unique, filtered random ssDNA candidates (fast; no folding).

    Returns a DataFrame with columns: core, full_sequence, length.
    Filters each randomized core on GC range, max homopolymer run, max single-base
    fraction, and (optionally) presence of all four bases; enforces full uniqueness.
    """
    import random, re, pandas as pd
    if fwd is None: fwd = DEFAULT_FWD
    if rev is None: rev = DEFAULT_REV
    run = max_homopolymer + 1
    homo = re.compile("(A{%d,}|C{%d,}|G{%d,}|T{%d,})" % (run, run, run, run))
    lo, hi = gc_range
    def valid(s):
        gc = (s.count("G") + s.count("C")) / len(s)
        if not (lo <= gc <= hi): return False
        if homo.search(s): return False
        for b in BASES:
            if s.count(b) / len(s) > max_base_frac: return False
            if require_all_bases and s.count(b) == 0: return False
        return True
    rng = random.Random(seed)
    cores, seen = [], set()
    while len(cores) < n:
        c = "".join(rng.choice(BASES) for _ in range(core_len))
        if c in seen: continue
        if valid(c):
            seen.add(c); cores.append(c)
    df = pd.DataFrame({"core": cores})
    df["full_sequence"] = fwd + df["core"] + rev
    df["length"] = df["full_sequence"].str.len()
    return df


def characterize_library(df, fwd=None, rev=None, core_col="core", temp=37.0,
                         workers=None, add_score=True):
    """Fold + characterize each candidate (needs seqfold + biopython).

    Adds dG_kcal_mol, Tm_C, GC_full_pct, GC_core_pct, Gquad_motif and (if
    add_score) a structure-informed `score`, then returns df sorted by score
    with rank + aptamer_id. Folding runs in parallel via a subprocess script.
    """
    import pandas as pd
    if fwd is None: fwd = DEFAULT_FWD
    if rev is None: rev = DEFAULT_REV
    if workers is None: workers = min(24, (os.cpu_count() or 4))
    cores = df[core_col].tolist()
    tmpd = tempfile.mkdtemp(prefix="seqgen_")
    cin, cout = os.path.join(tmpd, "cores.json"), os.path.join(tmpd, "folded.json")
    spath = os.path.join(tmpd, "fold_lib.py")
    json.dump({"fwd": fwd, "rev": rev, "cores": cores, "temp": temp,
               "workers": workers}, open(cin, "w"))
    open(spath, "w").write(FOLD_SCRIPT)
    r = subprocess.run([sys.executable, spath, cin, cout],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("folding failed: " + r.stderr[-2000:])
    fold = pd.DataFrame(json.load(open(cout))).sort_values("idx").drop(columns="idx")
    out = pd.concat([df.reset_index(drop=True), fold.reset_index(drop=True)], axis=1)
    if add_score:
        out = rank_library(out)
    return out


def rank_library(df, dg_col="dG_kcal_mol", gc_core_col="GC_core_pct",
                 gc_target=50.0, gc_tol=10.0, gquad_col="Gquad_motif",
                 id_prefix="Apt"):
    """Score + sort by a structure-informed heuristic.

    score = 0.70*fold_stability + 0.25*GC-closeness-to-target + 0.05 (+0.10 if Gquad).
    Prioritization only — not a measured binding affinity.
    """
    import numpy as np
    d = df.copy()
    g = d[dg_col].astype(float)
    dg_score = (g.max() - g) / (g.max() - g.min()) if g.max() != g.min() else 0.0
    gc_close = (1 - (d[gc_core_col].astype(float) - gc_target).abs() / gc_tol).clip(0, 1)
    bonus = d[gquad_col].astype(float) * 0.10 if gquad_col in d else 0.0
    d["score"] = (0.70 * dg_score + 0.25 * gc_close + 0.05).round(4) + bonus
    d = d.sort_values("score", ascending=False).reset_index(drop=True)
    d.insert(0, "rank", np.arange(1, len(d) + 1))
    d.insert(1, "aptamer_id", ["%s_%05d" % (id_prefix, i) for i in range(1, len(d) + 1)])
    return d
