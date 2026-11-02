"""Geometry helpers for concave epitopes: frame, rim, cavity profile, insertion."""
import numpy as np

PROBE_WATER = 1.4
VDW_HEAVY = 1.7
GRID_STEP = 0.5


def heavy_atoms(model, chains, resid_range=None):
    out = []
    for ch in chains:
        if ch not in model:
            continue
        for r in model[ch]:
            if r.id[0] != ' ':
                continue
            if resid_range is not None and r.id[1] not in resid_range:
                continue
            for a in r:
                if a.element != 'H':
                    out.append((ch, r.id[1], a.coord))
    return out


def groove_frame(model, floor_resids, wall_a_resids, wall_b_resids,
                 ligand_chain=None, wall_chain='A', floor_chain='A'):
    """Return (origin, u, v, n). All residue ranges are target-specific."""
    if ligand_chain is not None and ligand_chain in model:
        pts = np.array([r['CA'].coord for r in model[ligand_chain]
                        if r.id[0] == ' ' and 'CA' in r])
    else:
        pts = np.array([a[2] for a in heavy_atoms(model, [wall_chain])])
    c0 = pts.mean(0)
    u = np.linalg.svd(pts - c0)[2][0]
    u = u / np.linalg.norm(u)
    idx = {r.id[1] for r in model[wall_chain]}
    a1 = np.array([model[wall_chain][i]['CA'].coord for i in wall_a_resids if i in idx])
    a2 = np.array([model[wall_chain][i]['CA'].coord for i in wall_b_resids if i in idx])
    v = a2.mean(0) - a1.mean(0)
    v = v - (v @ u) * u
    v = v / np.linalg.norm(v)
    n = np.cross(u, v)
    n = n / np.linalg.norm(n)
    fid = {r.id[1] for r in model[floor_chain]}
    fl = np.array([model[floor_chain][i]['CA'].coord for i in floor_resids if i in fid])
    if (fl.mean(0) - c0) @ n > 0:
        n = -n
    return c0, u, v, n


def rim_height(model, frame, wall_a_resids, wall_b_resids, wall_chain='A', u_half=12.0):
    """Height of the LOWER wall crest within +-u_half along the long axis."""
    c0, u, v, n = frame
    at = heavy_atoms(model, [wall_chain])
    X = np.array([a[2] for a in at])
    h = (X - c0) @ n
    near = np.abs((X - c0) @ u) < u_half
    m1 = np.array([a[1] in wall_a_resids for a in at])
    m2 = np.array([a[1] in wall_b_resids for a in at])
    return float(min(h[near & m1].max(), h[near & m2].max()))


def cavity_profile(model, frame, rim, target_chains, v_lo, v_hi,
                   u_half=12.0, depths=None, u_bin=3.0):
    """Open volume (A^3) below the rim, by depth and by position along u."""
    if depths is None:
        depths = [0.0, 2.0, 4.0, 6.0, 8.0]
    from scipy.spatial import cKDTree
    c0, u, v, n = frame
    X = np.array([a[2] for a in heavy_atoms(model, target_chains)])
    tree = cKDTree(X)
    gu = np.arange(-u_half, u_half + 1e-9, GRID_STEP)
    gv = np.arange(v_lo, v_hi + 1e-9, GRID_STEP)
    gh = np.arange(-6.0, rim + 1e-9, GRID_STEP)
    G = np.stack(np.meshgrid(gu, gv, gh, indexing='ij'), -1).reshape(-1, 3)
    pts = c0 + G[:, 0:1] * u + G[:, 1:2] * v + G[:, 2:3] * n
    free = tree.query(pts)[0] > (VDW_HEAVY + PROBE_WATER)
    vox = GRID_STEP ** 3
    rows = []
    edges = np.arange(-u_half, u_half + 1e-9, u_bin)
    for d in depths:
        sel = free & (G[:, 2] < rim - d)
        row = {'depth_below_rim_A': float(d), 'open_volume_A3': round(float(sel.sum() * vox), 1)}
        for i in range(len(edges) - 1):
            b = sel & (G[:, 0] >= edges[i]) & (G[:, 0] < edges[i + 1])
            row['u_%d_%d' % (edges[i], edges[i + 1])] = round(float(b.sum() * vox), 1)
        rows.append(row)
    return rows


def insertion(model, binder_chains, frame, rim, v_lo, v_hi, u_half=12.0):
    """Depth stats for binder_chains. Compare against a reference binder."""
    c0, u, v, n = frame
    at = heavy_atoms(model, list(binder_chains))
    if not at:
        return {'max_depth': float('nan'), 'n_below_rim': 0, 'u_of_deepest': float('nan')}
    X = np.array([a[2] for a in at])
    h = (X - c0) @ n
    uu = (X - c0) @ u
    vv = (X - c0) @ v
    inside = (vv > v_lo) & (vv < v_hi) & (np.abs(uu) < u_half)
    if not inside.any():
        return {'max_depth': float('nan'), 'n_below_rim': 0, 'u_of_deepest': float('nan')}
    j = np.where(inside)[0][h[inside].argmin()]
    return {'max_depth': round(float(rim - h[inside].min()), 2),
            'n_below_rim': int(((h < rim) & inside).sum()),
            'u_of_deepest': round(float(uu[j]), 1)}


def pocket_lining(model, frame, rim, target_chains, v_lo, v_hi,
                  min_depth=4.0, u_min=None, u_half=12.0, cutoff=6.0):
    """Residues within `cutoff` of the deep sub-pocket, nearest first."""
    from scipy.spatial import cKDTree
    c0, u, v, n = frame
    at = heavy_atoms(model, target_chains)
    X = np.array([a[2] for a in at])
    tree = cKDTree(X)
    gu = np.arange(-u_half if u_min is None else u_min, u_half + 1e-9, GRID_STEP)
    gv = np.arange(v_lo, v_hi + 1e-9, GRID_STEP)
    gh = np.arange(-6.0, rim - min_depth + 1e-9, GRID_STEP)
    G = np.stack(np.meshgrid(gu, gv, gh, indexing='ij'), -1).reshape(-1, 3)
    pts = c0 + G[:, 0:1] * u + G[:, 1:2] * v + G[:, 2:3] * n
    deep = pts[tree.query(pts)[0] > (VDW_HEAVY + PROBE_WATER)]
    if not len(deep):
        return []
    dt = cKDTree(deep)
    best = {}
    for ch, ri, co in at:
        d = float(dt.query(co)[0])
        if d < cutoff:
            k = (ch, ri)
            best[k] = min(best.get(k, 99.0), round(d, 1))
    return sorted(({'chain': k[0], 'resid': k[1], 'closest_A': d}
                   for k, d in best.items()), key=lambda r: r['closest_A'])
