"""Helpers for designing and validating pH-selective protein binders.

Thermodynamics, replicate statistics, interface descriptors and developability flags.
Structural functions need Biopython and NumPy.
"""

import math

RT_KCAL_310K = 0.6163313981374999  # R*T at 310.15 K, kcal/mol
LIVE_LINKAGE_MIN = 1.05
INERT_LINKAGE_LO = 0.95
INERT_LINKAGE_HI = 1.05


def ph_window_ceiling(ph_lo=6.5, ph_hi=7.4, max_sites=4):
    """Maximum achievable affinity window for 1..max_sites ideal titratable sites.

    One ideal site contributes at most 10**(ph_hi - ph_lo): fully protonated when
    bound, fully deprotonated when free. Returns a list of dicts.
    """
    per = 10.0 ** (ph_hi - ph_lo)
    out = []
    for n in range(1, max_sites + 1):
        w = per ** n
        out.append({"n_ideal_sites": n,
                    "max_window": round(w, 1),
                    "ddG_kcal": round(RT_KCAL_310K * math.log(w), 2)})
    return out


def sites_needed(target_window=100.0, ph_lo=6.5, ph_hi=7.4):
    """How many IDEAL titratable sites a target window implies.

    Fractional by design: 2.22 means three real histidines, since no real site
    reaches the ideal limit.
    """
    per = 10.0 ** (ph_hi - ph_lo)
    return math.log(target_window) / math.log(per)


def site_linkage(pka_free, pka_bound, ph_lo=6.5, ph_hi=7.4):
    """Ka(ph_lo)/Ka(ph_hi) contributed by ONE ionizable site.

    >1 favours binding at the acidic pH; ~1 is thermodynamically inert;
    <1 anti-selects (favours the alkaline pH).
    """
    def f(pk, ph):
        return 1.0 + 10.0 ** (pk - ph)
    return (f(pka_bound, ph_lo) / f(pka_free, ph_lo)) / \
           (f(pka_bound, ph_hi) / f(pka_free, ph_hi))


def classify_site(linkage):
    """Label one site's linkage: 'live', 'inert' or 'inverted'."""
    if linkage > LIVE_LINKAGE_MIN:
        return "live"
    if linkage < INERT_LINKAGE_LO:
        return "inverted"
    return "inert"


def combined_window(sites, ph_lo=6.5, ph_hi=7.4):
    """Overall window from a list of (pka_free, pka_bound) pairs for ONE complex.

    Only valid within a single measured structure: per-site linkages do NOT
    transfer between sequence variants, so never combine across designs.
    """
    lks = [site_linkage(f, b, ph_lo, ph_hi) for f, b in sites]
    total = 1.0
    for v in lks:
        total *= v
    return {"window": total,
            "ddG_kcal": RT_KCAL_310K * math.log(total) if total > 0 else float("nan"),
            "per_site_linkage": lks,
            "per_site_class": [classify_site(v) for v in lks],
            "n_live": sum(1 for v in lks if v > LIVE_LINKAGE_MIN),
            "n_inert": sum(1 for v in lks if INERT_LINKAGE_LO < v < INERT_LINKAGE_HI),
            "n_inverted": sum(1 for v in lks if v < INERT_LINKAGE_LO),
            "direction": "acid-preferring" if total > 1 else "inverted"}


def replicate_verdict(windows, require_all_above=1.0):
    """Summarise one design's window across independent predicted poses.

    Pass the window measured on EVERY pose the predictor emitted. Reproducibility,
    not the maximum, is what distinguishes signal from pose noise.
    """
    ws = [float(w) for w in windows if w is not None and w > 0]
    if not ws:
        return {"n": 0, "reproducible": False}
    logs = [math.log(w) for w in ws]
    gm = math.exp(sum(logs) / len(logs))
    mean = sum(logs) / len(logs)
    if len(logs) > 1:
        var = sum((x - mean) ** 2 for x in logs) / (len(logs) - 1)
        sd = math.sqrt(var)
    else:
        sd = 0.0
    return {"n": len(ws),
            "window_geomean": gm,
            "window_min": min(ws),
            "window_max": max(ws),
            "fold_spread": max(ws) / min(ws),
            "log_sd": sd,
            "reproducible": all(w > require_all_above for w in ws),
            "ddG_kcal": RT_KCAL_310K * math.log(gm)}


def noise_floor(free_pkas, bound_pkas, ph_lo=6.5, ph_hi=7.4, lo_pct=10, hi_pct=90):
    """Noise band implied by an INVARIANT control residue across many complexes.

    Pass the predicted pKa of a titratable residue that is chemically identical in
    every complex (a target-side histidine is ideal). A tight free-state spread with
    a wide bound-state spread means the pKa predictor is self-consistent and the
    variation is pose-driven. If your best design's window sits inside the returned
    band, your ranking is substantially noise.
    """
    def pct(xs, p):
        s = sorted(xs)
        if not s:
            return float("nan")
        k = (len(s) - 1) * (p / 100.0)
        lo_i = int(math.floor(k))
        hi_i = int(math.ceil(k))
        if lo_i == hi_i:
            return s[lo_i]
        return s[lo_i] + (s[hi_i] - s[lo_i]) * (k - lo_i)

    def sd(xs):
        if len(xs) < 2:
            return 0.0
        m = sum(xs) / len(xs)
        return math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))

    fmed = pct(free_pkas, 50)
    blo, bhi = pct(bound_pkas, lo_pct), pct(bound_pkas, hi_pct)
    w_lo = site_linkage(fmed, blo, ph_lo, ph_hi)
    w_hi = site_linkage(fmed, bhi, ph_lo, ph_hi)
    return {"n": len(bound_pkas),
            "free_pka_sd": sd(list(free_pkas)),
            "bound_pka_sd": sd(list(bound_pkas)),
            "free_pka_range": (min(free_pkas), max(free_pkas)) if free_pkas else None,
            "bound_pka_range": (min(bound_pkas), max(bound_pkas)) if bound_pkas else None,
            "implied_window_band": (min(w_lo, w_hi), max(w_lo, w_hi)),
            "band_width_fold": max(w_lo, w_hi) / min(w_lo, w_hi)}


def developability_flags(seq, parent_seq=None, max_met=6, avoid_positions=None):
    """Sequence-level liability flags. Apply as a HARD gate ahead of the window."""
    import re
    n_cys = seq.count("C")
    flags = {"length": len(seq),
             "n_Met": seq.count("M"),
             "n_Cys": n_cys,
             "unpaired_Cys": n_cys % 2,
             "non_canonical_Cys": max(0, n_cys - 2),
             "n_glyc_sequons": len(re.findall(r"N[^P][ST]", seq)),
             "max_hydrophobic_run": max((len(m) for m in re.findall(r"[MIVLFWY]+", seq)),
                                        default=0),
             "n_DG_isomerisation": len(re.findall(r"DG", seq)),
             "n_DP_cleavage": len(re.findall(r"DP", seq))}
    if parent_seq is not None:
        flags["n_mut_vs_parent"] = sum(1 for a, b in zip(seq, parent_seq) if a != b)
        flags["excess_Met_vs_parent"] = seq.count("M") - parent_seq.count("M")
    if avoid_positions:
        flags["hits_avoid_positions"] = [p for p in avoid_positions
                                         if 1 <= p <= len(seq)]
    flags["clean"] = (flags["unpaired_Cys"] == 0
                      and flags["non_canonical_Cys"] == 0
                      and flags["n_glyc_sequons"] == 0
                      and flags["n_Met"] <= max_met)
    return flags


def interface_descriptors(pdb_path, target_chain="A", binder_chain="B",
                          contact_cut=4.5, hbond_cut=3.5, saltbridge_cut=4.0):
    """Classical interface descriptors from a complex structure.

    Compare against the SAME quantities computed on a natural binder for the same
    target. Natural antibody-antigen interfaces bury roughly 600-900 A^2 per side at
    an apolar fraction near 0.45; designed interfaces drifting to 0.65-0.75 carry
    aggregation and off-target risk.
    """
    import copy
    import numpy as np
    from Bio.PDB import PDBParser, ShrakeRupley
    from Bio.PDB.Polypeptide import is_aa

    apolar = set("ALA VAL LEU ILE MET PHE TRP PRO GLY CYS TYR".split())
    pos_atoms = {"LYS": {"NZ"}, "ARG": {"NH1", "NH2", "NE"}, "HIS": {"ND1", "NE2"}}
    neg_atoms = {"ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}}

    model = PDBParser(QUIET=True).get_structure("x", pdb_path)[0]
    sr = ShrakeRupley()

    cx = copy.deepcopy(model)
    sr.compute(cx, level="R")
    t_only = copy.deepcopy(model)
    t_only.detach_child(binder_chain)
    sr.compute(t_only, level="R")
    b_only = copy.deepcopy(model)
    b_only.detach_child(target_chain)
    sr.compute(b_only, level="R")

    def total(struct, ch):
        return sum(r.sasa for r in struct[ch] if is_aa(r, standard=True))

    bsa = (total(t_only, target_chain) + total(b_only, binder_chain)
           - total(cx, target_chain) - total(cx, binder_chain))

    tres = [r for r in model[target_chain] if is_aa(r, standard=True)]
    bres = [r for r in model[binder_chain] if is_aa(r, standard=True)]
    batoms = [(r, a) for r in bres for a in r]

    n_contact = 0
    n_hb = 0
    n_sb = 0
    iface_t = set()
    iface_b = set()
    for rt in tres:
        for at in rt:
            for rb, ab in batoms:
                d = float(np.linalg.norm(at.coord - ab.coord))
                if d > contact_cut:
                    continue
                n_contact += 1
                iface_t.add(rt.id[1])
                iface_b.add(rb.id[1])
                if d <= hbond_cut and at.element in ("N", "O") and ab.element in ("N", "O"):
                    n_hb += 1
                if d <= saltbridge_cut:
                    tn = rt.get_resname()
                    bn = rb.get_resname()
                    if (at.get_id() in pos_atoms.get(tn, ()) and
                            ab.get_id() in neg_atoms.get(bn, ())):
                        n_sb += 1
                    elif (at.get_id() in neg_atoms.get(tn, ()) and
                            ab.get_id() in pos_atoms.get(bn, ())):
                        n_sb += 1

    free_b = {r.id[1]: r.sasa for r in b_only[binder_chain] if is_aa(r, standard=True)}
    bound_b = {r.id[1]: r.sasa for r in cx[binder_chain] if is_aa(r, standard=True)}
    buried_ap = 0.0
    buried_tot = 0.0
    for r in bres:
        d = free_b.get(r.id[1], 0.0) - bound_b.get(r.id[1], 0.0)
        if d > 1.0:
            buried_tot += d
            if r.get_resname() in apolar:
                buried_ap += d

    return {"BSA_total": round(bsa, 1),
            "BSA_per_side": round(bsa / 2.0, 1),
            "n_contacts": n_contact,
            "n_interface_residues_target": len(iface_t),
            "n_interface_residues_binder": len(iface_b),
            "n_hbonds": n_hb,
            "n_saltbridges": n_sb,
            "apolar_fraction_buried": round(buried_ap / buried_tot, 3) if buried_tot else None,
            "interface_residues_target": sorted(iface_t)}


def find_carboxylate_clusters(pdb_path, chain="A", max_span=20.0, min_sasa=20.0,
                              cluster_size=3, exclude_paired=True):
    """Find clusters of solvent-exposed, unpaired carboxylates on a target surface.

    These are the only viable partners for protonated switch histidines, so the
    clusters this returns bound how many live sites the target can host. A cluster
    being present does not mean your binder can reach it - check CDR distances too.

    Residue numbers are reported in the COORDINATE FILE's own numbering, which for a
    truncated construct is usually not database numbering. Map to your reference
    numbering before comparing with an epitope definition (see SKILL.md section 8).
    """
    import copy
    import itertools
    import numpy as np
    from Bio.PDB import PDBParser, ShrakeRupley
    from Bio.PDB.Polypeptide import is_aa

    model = PDBParser(QUIET=True).get_structure("x", pdb_path)[0]
    only = copy.deepcopy(model)
    for ch in list(only):
        if ch.id != chain:
            only.detach_child(ch.id)
    sr = ShrakeRupley()
    sr.compute(only, level="R")
    sasa = {r.id[1]: r.sasa for r in only[chain] if is_aa(r, standard=True)}

    acid_atoms = ("OD1", "OD2", "OE1", "OE2")
    base_atoms = ("NZ", "NH1", "NH2", "NE")
    residues = [r for r in model[chain] if is_aa(r, standard=True)]
    basics = [(r, [a for a in r if a.get_id() in base_atoms])
              for r in residues if r.get_resname() in ("LYS", "ARG")]

    cands = []
    for r in residues:
        if r.get_resname() not in ("ASP", "GLU"):
            continue
        ats = [a for a in r if a.get_id() in acid_atoms]
        if not ats:
            continue
        if sasa.get(r.id[1], 0.0) < min_sasa:
            continue
        paired = False
        if exclude_paired:
            for rb, bats in basics:
                if rb.id[1] == r.id[1]:
                    continue
                for a in ats:
                    for b in bats:
                        if float(np.linalg.norm(a.coord - b.coord)) < 4.0:
                            paired = True
                            break
                    if paired:
                        break
                if paired:
                    break
        if paired:
            continue
        cands.append({"resid": r.id[1], "resname": r.get_resname(),
                      "sasa": round(sasa.get(r.id[1], 0.0), 1),
                      "centroid": np.mean([a.coord for a in ats], axis=0)})

    out = []
    for combo in itertools.combinations(range(len(cands)), cluster_size):
        span = 0.0
        for i, j in itertools.combinations(combo, 2):
            d = float(np.linalg.norm(cands[i]["centroid"] - cands[j]["centroid"]))
            span = max(span, d)
        if span <= max_span:
            out.append({"cluster": [f"{cands[k]['resname']}{cands[k]['resid']}"
                                    for k in combo],
                        "resids": [cands[k]["resid"] for k in combo],
                        "max_span_A": round(span, 1),
                        "min_sasa": min(cands[k]["sasa"] for k in combo),
                        "total_sasa": round(sum(cands[k]["sasa"] for k in combo), 1)})
    out.sort(key=lambda d: d["max_span_A"])
    return {"n_candidate_carboxylates": len(cands),
            "candidates": [{k: v for k, v in c.items() if k != "centroid"} for c in cands],
            "n_clusters": len(out), "clusters": out}
