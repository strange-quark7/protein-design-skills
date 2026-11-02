"""Full-atom interface measurements for targeted binder campaigns.

Every function here expects a FOLDED, full-atom complex. Measuring contacts,
burial or reachability on a generated (poly-glycine, backbone-only) structure
is the error this skill exists to prevent -- see SKILL.md.
"""

import os
import glob

BACKBONE_ATOMS = ("N", "CA", "C", "O")
REACH_WINDOW = (3.5, 7.5)


def load_complex(pdb_path, model_index=0):
    """Parse a PDB/mmCIF complex and return the requested model."""
    from Bio.PDB import PDBParser, MMCIFParser
    parser = MMCIFParser(QUIET=True) if pdb_path.lower().endswith(
        (".cif", ".mmcif")) else PDBParser(QUIET=True)
    return parser.get_structure("x", pdb_path)[model_index]


def residue_map(model, chain_id):
    """{resid: Residue} for standard residues of one chain."""
    return {r.id[1]: r for r in model[chain_id] if r.id[0] == " "}


def virtual_cb(residue):
    """CB position from N/CA/C (works for Gly). None if backbone incomplete."""
    import numpy as np
    try:
        n, ca, c = (np.array(residue[a].coord) for a in ("N", "CA", "C"))
    except KeyError:
        return None
    b, cc = ca - n, c - ca
    a = np.cross(b, cc)
    return -0.58273431 * a + 0.56802827 * b - 0.54067466 * cc + ca


def binder_interface(pdb_path, target_chain="A", binder_chain="B", cutoff=5.0):
    """Per-target-residue minimum heavy-atom distance to the binder.

    Returns {'distances': {resid: float}, 'contacts': [resid, ...],
             'n_contacts': int}.
    """
    import numpy as np
    model = load_complex(pdb_path)
    tgt = residue_map(model, target_chain)
    binder = np.array([a.coord for r in model[binder_chain]
                       if r.id[0] == " " for a in r])
    out = {}
    for resid, res in tgt.items():
        at = np.array([a.coord for a in res])
        out[resid] = float(np.linalg.norm(
            at[:, None, :] - binder[None, :, :], axis=-1).min())
    hits = sorted(q for q, d in out.items() if d < cutoff)
    return {"distances": out, "contacts": hits, "n_contacts": len(hits)}


def epitope_fidelity(pdb_path, epitope, hotspots=None, target_chain="A",
                     binder_chain="B", cutoff=5.0):
    """Did the prediction bind the epitope it was designed against?

    High interface confidence does not mean the intended epitope -- a design
    can score well with its partner far off target. Check `off_epitope`: if
    those residues are numerically adjacent to the epitope the footprint merely
    grew; scattered values mean relocation.
    """
    if hotspots is None:
        hotspots = []
    iface = binder_interface(pdb_path, target_chain, binder_chain, cutoff)
    contacts = set(iface["contacts"])
    epi = set(epitope)
    on = contacts & epi
    engaged = [q for q in hotspots if iface["distances"].get(q, 1e9) < cutoff]
    return {
        "n_contact": len(contacts),
        "on_epitope": sorted(on),
        "on_epitope_frac": round(len(on) / max(len(contacts), 1), 3),
        "off_epitope": sorted(contacts - epi),
        "hotspots_engaged": engaged,
        "n_hotspots_engaged": len(engaged),
        "hotspot_distances": {q: round(iface["distances"].get(q, float("nan")), 2)
                              for q in hotspots},
    }


def interface_conservation(pdb_path, conserved_resids, target_chain="A",
                           binder_chain="B", cutoff=5.0):
    """Conserved fraction of the interface the fold actually made.

    Measure cross-species reactivity here, not on a generated backbone -- a
    backbone-level 'avoids residue X' filter was 0% predictive in one campaign.
    Compare against a known cross-reactive binder folded the same way.
    """
    iface = binder_interface(pdb_path, target_chain, binder_chain, cutoff)
    contacts = iface["contacts"]
    cons = set(conserved_resids)
    kept = [q for q in contacts if q in cons]
    return {
        "n_contact": len(contacts),
        "n_conserved": len(kept),
        "conserved_frac": round(len(kept) / max(len(contacts), 1), 3),
        "variable_contacts": [q for q in contacts if q not in cons],
    }


def contact_occupancy(pose_paths, target_chain="A", target_resid=None,
                      binder_chain="B", target_atoms=None, binder_resname=None,
                      binder_atoms=None, threshold=3.5):
    """How often a specific contact is present across resampled poses.

    The gate that separates a converged feature from a lucky pose. A single
    prediction's distance does not predict occupancy (rho = +0.21 in one
    campaign). Report occupancy plus spread, never one distance.

    `target_atoms` / `binder_atoms` restrict to named atoms, e.g.
    ("OD1", "OD2") for an Asp carboxylate and ("ND1", "NE2") for a His ring.
    """
    import numpy as np
    dists = []
    for path in pose_paths:
        model = load_complex(path)
        tres = residue_map(model, target_chain).get(target_resid)
        if tres is None:
            continue
        tat = np.array([a.coord for a in tres
                        if target_atoms is None or a.get_name() in target_atoms])
        best = float("inf")
        for res in model[binder_chain]:
            if res.id[0] != " ":
                continue
            if binder_resname is not None and res.get_resname() != binder_resname:
                continue
            bat = [a.coord for a in res
                   if binder_atoms is None or a.get_name() in binder_atoms]
            if not bat:
                continue
            d = float(np.linalg.norm(
                np.array(bat)[:, None, :] - tat[None, :, :], axis=-1).min())
            best = min(best, d)
        if best < float("inf"):
            dists.append(round(best, 2))
    if not dists:
        return {"n_poses": 0, "occupancy": None, "distances": []}
    arr = np.array(dists)
    held = arr[arr <= threshold]
    return {
        "n_poses": len(dists),
        "n_present": int((arr <= threshold).sum()),
        "occupancy": round(float((arr <= threshold).mean()), 3),
        "median": round(float(np.median(arr)), 2),
        "spread_when_present": round(float(held.max() - held.min()), 2) if len(held) else None,
        "d_min": round(float(arr.min()), 2),
        "d_max": round(float(arr.max()), 2),
        "distances": dists,
    }


def reachable_positions(pdb_path, target_chain="A", target_resid=None,
                        binder_chain="B", target_atoms=None, window=None):
    """Binder positions whose CB could reach a required target contact.

    A side chain extends ~4 A past CB, so positions with CB in the 3.5-7.5 A
    window can form a hydrogen bond or salt bridge to the target atom if the
    right residue is placed there. Use this to decide whether a required
    contact needs a per-residue bias (cheap) or new backbones (expensive).
    """
    import numpy as np
    if window is None:
        window = REACH_WINDOW
    lo, hi = window
    model = load_complex(pdb_path)
    tres = residue_map(model, target_chain).get(target_resid)
    if tres is None:
        return []
    tat = np.array([a.coord for a in tres
                    if target_atoms is None or a.get_name() in target_atoms])
    out = []
    for res in model[binder_chain]:
        if res.id[0] != " ":
            continue
        cb = res["CB"].coord if "CB" in res else virtual_cb(res)
        if cb is None:
            continue
        d = float(np.linalg.norm(np.array(cb)[None, :] - tat, axis=-1).min())
        if lo <= d <= hi:
            out.append({"resid": res.id[1], "resname": res.get_resname(),
                        "cb_distance": round(d, 2)})
    return sorted(out, key=lambda r: r["cb_distance"])


def composition_report(sequences, reference=None):
    """Residue fractions versus a reference set, to catch composition collapse.

    Default sequence-design settings on idealised backbones return poly-alanine
    (38% Ala, 58% hydrophobic in one campaign) -- invisible by eye, obvious
    against real expressible proteins of the same class. Pass `reference` as a
    list of such sequences.
    """
    from collections import Counter

    def profile(seqs):
        joined = "".join(seqs)
        c = Counter(joined)
        n = max(len(joined), 1)
        return {
            "n_seqs": len(seqs),
            "mean_length": round(sum(len(s) for s in seqs) / max(len(seqs), 1), 1),
            "Ala": round(c["A"] / n, 3),
            "Gly": round(c["G"] / n, 3),
            "Cys": round(c["C"] / n, 3),
            "hydrophobic": round(sum(c[a] for a in "AVLIMFWY") / n, 3),
            "charged": round(sum(c[a] for a in "DEKR") / n, 3),
            "Glu": round(c["E"] / n, 3),
            "net_charge": c["K"] + c["R"] - c["D"] - c["E"],
            "distinct": len(c),
        }

    out = {"designs": profile(list(sequences))}
    if reference:
        out["reference"] = profile(list(reference))
        out["flags"] = [k for k in ("Ala", "hydrophobic", "charged")
                        if out["designs"][k] > 1.5 * max(out["reference"][k], 1e-6)]
    return out


def count_predictions(root, pattern=None):
    """Count written predictions -- never trust a co-folder's exit code.

    Boltz exits 0 on input errors while writing nothing, so an exit-code loop
    reports success for an empty run. Implausibly fast runs are the other tell.
    """
    if pattern is None:
        pattern = "confidence_*_model_0.json"
    return len(glob.glob(os.path.join(root, "**", pattern), recursive=True))
