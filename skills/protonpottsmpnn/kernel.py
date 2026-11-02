"""Thermodynamic helpers for pH-selective binder design (ProtonPottsMPNN)."""

RT_KCAL_310K = 0.6163313981374999  # R*T at 310.15 K, kcal/mol


def ph_window_ceiling(ph_lo=6.5, ph_hi=7.4, max_sites=4):
    """Max achievable affinity window Ka(ph_lo)/Ka(ph_hi) for n ideal titratable sites.

    One ideal site (fully protonated when bound, fully deprotonated when free)
    gives 10**(ph_hi - ph_lo). Returns a list of dicts, one row per site count.
    """
    import math
    per_site = 10.0 ** (ph_hi - ph_lo)
    rows = []
    for n in range(1, int(max_sites) + 1):
        w = per_site ** n
        rows.append({"n_ideal_sites": n,
                     "max_window": round(w, 1),
                     "ddG_kcal": round(RT_KCAL_310K * math.log(w), 2)})
    return rows


def sites_needed(window, ph_lo=6.5, ph_hi=7.4):
    """How many ideal titratable sites a target affinity window implies."""
    import math
    return math.log(window) / math.log(10.0 ** (ph_hi - ph_lo))


def site_linkage(pka_free, pka_bound, ph_lo=6.5, ph_hi=7.4):
    """Linkage factor for ONE site: Ka(ph_lo)/Ka(ph_hi) from its free/bound pKa.

    >1 acid-preferring (what a tumour-pH switch needs), ~1 inert, <1 inverted.
    """
    def frac(pka, ph):
        return (1.0 + 10.0 ** (pka - ph))
    return (frac(pka_bound, ph_lo) / frac(pka_free, ph_lo)) / \
           (frac(pka_bound, ph_hi) / frac(pka_free, ph_hi))


def combined_window(sites, ph_lo=6.5, ph_hi=7.4):
    """Total window from several independent sites.

    `sites` is an iterable of (pka_free, pka_bound) pairs. Returns a dict with the
    product window, ddG, per-site linkages, and which sites were effectively inert.
    """
    import math
    lks = [site_linkage(f, b, ph_lo, ph_hi) for f, b in sites]
    total = 1.0
    for v in lks:
        total *= v
    return {"window": total,
            "ddG_kcal": RT_KCAL_310K * math.log(total) if total > 0 else float("nan"),
            "per_site_linkage": lks,
            "n_inert_sites": sum(1 for v in lks if 0.95 < v < 1.05),
            "direction": "acid-preferring" if total > 1 else "inverted"}
