---
name: protonpottsmpnn
description: "Design pH-dependent (acid- or base-selective) protein binders with ProtonPottsMPNN, a PottsMPNN variant with an explicit protonation alphabet (HIS-P/HIS-S, ASP-P/ASP-D, GLU-P/GLU-D) that optimises the selectivity gap E_P - E_D directly. Reach for this skill to engineer a pH switch into a binder, nanobody, or interface - e.g. a binder that engages at tumour pH 6.5 but not at blood pH 7.4 - to place protonated histidine/carboxylate centres at an interface, or to trace a stability-versus-pH-selectivity Pareto front. Covers install on a wheel-only or Windows sandbox, the HBPLUS dependency, the three placement routes and how each of them fails, enumeration settings that decide the wall clock, and the validation gates that a raw design score does not cover."
---

# ProtonPottsMPNN — pH-selective binder design

ProtonPottsMPNN (Jacobsen et al. 2026, github.com/christian-creator/ProtonPottsMPNN)
is a Potts/MPNN inverse-folding model whose alphabet carries **explicit protonation
states**: `HIS-P`/`HIS-S`, `ASP-P`/`ASP-D`, `GLU-P`/`GLU-D`. Because protonated and
deprotonated forms are separate tokens, it can score a sequence under each and
optimise the **selectivity gap**

    E_selective = E_P - E_D

i.e. it designs *for pH-dependent binding* rather than biasing histidine content and
hoping. `combined_lambda` interpolates: `0.0` = pure Potts stability, `1.0` = pure
selectivity.

> **Methodology and findings:** the validation protocol, the nine failure modes and the
> quantitative evidence behind them are tool-agnostic and live in the
> `ph-switch-binder-validation` skill. Load that one for *how to design and validate a pH
> switch*; this one covers *how to install and drive this particular model*. Sections
> 1-9 below summarise the overlapping rules so this skill stands alone.

## When to use it

Use it to put a pH switch into an interface — acid-selective tumour-targeting
binders, pH-triggered release, endosomal-escape switches. For plain
(pH-independent) inverse folding use **ProteinMPNN** instead; this model is slower
and its objective is the pH gap, not sequence recovery.

## Inverting the polarity: designing for binding at the HIGHER pH

The default objective maximises `E_P - E_D`, which favours the **protonated**
histidine and therefore designs **acid-preferring** binders. A brief that asks for
binding at pH 7.4 and release at 6.0 (release-on-acidification, extracellular
retention, high-pH-selective capture) needs the opposite, and running the model
as-shipped optimises in precisely the wrong direction.

**The objective inverts without touching the engine** — swap which protonation
state is the centre and which is its dependent. The same maximised gap then
favours the neutral form:

| | `center_protonation_types` | `dep_map` | favours |
|---|---|---|---|
| standard | `["HIS-P", ...]` | `{"HIS-P": ["HIS-S"]}` | protonated → acid-preferring |
| **inverted** | `["HIS-S", ...]` | `{"HIS-S": ["HIS-P"]}` | neutral → **binds at high pH** |

**The v6 checkpoint's neutral-His token is `HIS-S`.** The `DEFAULT_DEP_MAP` in the
source says `{"HIS-P": ["HID", "HIE"]}`, which is stale — HID/HIE exist only in
v3/v4 and v6 raises `KeyError: ... absent from the model vocabulary` on both the
centre and the dep_map. The full v6 alphabet is `HIS`, `HIS-A`, `HIS-P`, `HIS-S`
plus `ASP-{A,D,P}` / `GLU-{A,D,P}`.

Verified contrast on one backbone, same engine, sign flipped — the inverted run
produces exactly the chemistry the mechanism predicts:

| near the placed centres | standard | inverted |
|---|---|---|
| acidic (D/E) | 13 per 24 residues | **0** |
| basic (K/R/H) | 3 | **8** |
| histidines placed | 4.5 | 8.5 |

The inverted objective surrounds the switch histidines with Arg/His clusters and
**no** carboxylates — positive charge that destabilises the protonated imidazolium
and drives bound-state pKa down, which is what makes binding favoured at the
higher pH. That is the one lever that attacks **free-state** pKa; see §1a on why
that matters more than adding sites.

Everything else in this skill still applies unchanged — above all §3: the model's
`selective_energy` does not predict the validated window in *either* polarity, so
these are candidate sequences only. Fold them and measure the linkage.

### 1a. Which pKa to push — the window saturates in one and not the other

Computing `1/L` across the (pka_free, pka_bound) plane shows the window is far
more sensitive to free-state pKa than to bound-state pKa. Over 6.0 → 7.4:

| pKa free | pKa bound | window per site | two sites |
|---|---|---|---|
| 6.5 | 4.5 | 3.6x | 13x |
| 6.8 | 4.5 | 5.7x | 32x |
| 6.8 | **3.0** | 5.8x | 34x |
| **7.2** | 4.5 | **10.0x** | **101x** |

Driving bound pKa from 4.5 to 3.0 buys almost nothing (5.67 → 5.84). Raising free
pKa from 6.8 to 7.2 nearly doubles the per-site window. **Burial alone caps a site
near 5.7x**, because burial acts on the bound state. If the arithmetic says you
need more than burial can give, the answer is a local environment that raises the
*free* pKa — not a third site.

## Read this before you trust a design

Hard-won on a real campaign — these are the failure modes, in the order they bite.

### 1. The physics sets a floor on how many sites you need

Over a pH interval `pH_hi - pH_lo`, **one ideal titratable site** can change affinity
by at most `10^(pH_hi - pH_lo)`. For 6.5 -> 7.4 that is `10^0.9` = **7.9x**
(1.28 kcal/mol). So:

| ideal sites | max window | meets "no detectable binding" (>=100x) |
|---|---|---|
| 1 | 7.9x | no |
| 2 | 63x | no |
| 3 | 501x | yes |

**A >=100x window needs >=2.2 ideal sites — practically three histidines, each near
its full theoretical shift.** Use `center_types=["HIS-P"]*3` (an *exact multiset*, see
below). Asking for one His plus two carboxylates typically yields one live site: at
a typical interface Asp/Glu sit near pKa 3.4-4.9, ~2 units below a 6.5-7.4 window, so
they are fully ionised at both ends and contribute a linkage factor of 1.00 — inert.

Call `ph_window_ceiling()` to get this table for your own pH pair.

### 2. Placement is the hardest part to get right — and on some systems it cannot be

With the default `placement_region=["all"]` the model will happily put protonated
centres in the scaffold framework, where they cannot gate binding. Pass
`placement_region=["interface"]` (CA-CA within 6 A of the target). In one campaign
this moved 2 of 3 centres out of the framework and into CDR3.

**But check that the interface mask is actually big enough first.** On a later campaign
all three placement routes failed on the same system:

| setting | outcome |
|---|---|
| `placement_region=["interface"]` | **zero designs returned**, no error. The CA-CA <= 6 A mask held only **1-5 binder positions** on these poses - fewer than the 2-3 centres requested, so no valid plan exists and the engine silently yields nothing. |
| `placement_region=["interface","surface"]` | designs produced, but centres landed **18.8-24.7 A** from the nearest carboxylate on both a VHH and an scFv scaffold - zero partnered sites |
| `explicit_centers=[{"res_id": r, "protonation_type": "HIS-P"}, ...]` | **hangs** - no output in 20 min on a complex whose normal placement path completes a solve in ~15 min |

Measure the mask before you trust a run: count binder residues whose CA lies within 6 A
of any target CA. A useful calibration point - a *native* antibody-antigen complex had
only **3** positions by that criterion, so a small mask is normal, not a sign of a bad
pose. If the mask is smaller than your centre count, `["interface"]` cannot work.

The widened region reproduces failure mode 5 exactly: the selectivity-gap objective does
not require an electrostatic partner, so the model places centres where nothing can
gate. **Always measure each returned centre's distance to the nearest conserved,
unpaired, exposed carboxylate and discard designs with no partnered site** - that filter
is cheap and it rejected 12 of 12 designs in that campaign.

**The fallback when placement cannot be steered:** mutate the partnered positions to
histidine directly and validate with the same fold-plus-pKa gate. It is a cruder design
than the model's, but it tests the physics without depending on the placement
machinery - and in that campaign it gave the decisive result the tool could not.

### 3. `selective_energy` does NOT predict the validated window

Measured on one campaign: the design with the best `selective_energy` (-15.1) had the
**worst** independently computed window (0.58x, i.e. inverted), while a middling
-12.8 gave the best. Treat `selective_energy` as a search objective, never as a
selection criterion. Every candidate needs an independent gate:

1. **Fold** the complex (Boltz-2, AF2-multimer, or RF3) — reject low ipTM.
2. **Epitope check** — reject anything that relocated off the intended face.
3. **pKa / linkage** (PROPKA on free vs bound states) — the actual window.

### 4. Gate on epitope BEFORE interpreting a window

The single most important ordering rule. A design that drifts off the intended
epitope can show an apparently favourable window that is just a property of a
different surface. In one campaign, predicted window was **anti-correlated** with
on-target epitope contact (r = -0.84, n = 7): every complex that genuinely engaged the
target face was pH-inverted, and the only two with window > 1 had **0%** on-target
contact. Both were false positives. Discard designs below ~50% on-target contact
before reading their pH numbers.

### 5. Check which centres actually ended up live

Compute per-site linkage from free/bound pKa (`site_linkage()`); a site with linkage
~1.00 contributed nothing. Buried histidines lose pKa to desolvation and invert the
switch (one case: pKa 5.92 free -> 0.28 bound). Favour **interface-rim** histidines
that retain solvation when bound. Useful empirical correlations from one interface
(n = 12 sites): dpKa vs retained solvation r = +0.71; dpKa vs burial r = -0.74.
Also avoid anchoring a switch on a *target* histidine without checking it — a target
His whose own pKa drops on binding anti-selects for acid.

### 6. Per-site linkage does NOT transfer between sequence variants

Do not measure per-site linkage, keep the live sites, mutate the inverted ones away,
and predict the new window as the product of the kept linkages. This was tested
directly and failed by ~40-fold: predicted 61.9x, measured 0.30-1.48x across six
repaired variants. Every retained site lost its linkage (one fell 4.88 -> 0.61). Even
a histidine on the *target* protein, never mutated, dropped from 3.91 to 0.69 after two
mutations on the binder 15+ A away. The window is an emergent property of the entire
bound interface. Fold and re-measure every variant; never arithmetic-combine linkages
across variants, and never present such a product as a predicted window.

### 7. Establish your noise floor before trusting any ranking

Pick a titratable residue that is identical in every complex (a target-side His is
ideal) and inspect the spread of its predicted pKa across the whole set. Across 46
complexes in one campaign its *free*-state pKa was tight (SD 0.05 - the pKa predictor
is self-consistent) while its *bound*-state pKa spanned 2.58-10.25 (SD 1.47), driven
entirely by predicted-pose variation. Propagating that spread yields an apparent
single-site window of 0.69-2.63x from noise alone. If the best design sits inside or
barely above that band, the ranking is substantially noise and further sampling
optimises against the noise rather than toward a real switch. Stop and either change
the physics engine (constant-pH MD) or go to experiment.

### 8. Count the sites the geometry can host, before designing

A >=100x window needs >=3 live histidines, hence >=3 carboxylate partners within one
paratope span. Search the target surface for that cluster first: mouse-conserved (if
cross-reactivity matters), not already salt-bridged, solvent-exposed, mutually within
~20 A. In one case the best cluster spanned 18.9 A but two of its three members sat
18.7-18.8 A from the nearest CDR position - beyond side-chain reach - so that pose
could host exactly one live site and no sequence design could fix it. That is a
re-docking problem, not an inverse-folding one, and this tool cannot solve it.

**Do not read that as "fix the geometry and the switch works".** A later cycle fixed it
and the switch still failed. Pocket-constrained co-folding plus a format change (122 aa
VHH -> 247 aa scFv, because a single-domain paratope cannot straddle an 18.3 A cluster)
brought all three carboxylates within reach; histidines placed beside them then gave
**0 of 21 live sites, 18 of 21 inverted**, with the 3-His variant more inverted than its
unmutated parent. Every engineered histidine lost pKa on binding, because burial
desolvates the imidazole faster than the neighbouring carboxylate stabilises it.

So the geometry count rules a site **out**, never in. The cheap test that should come
first: put one histidine at the single partner you can already reach and measure the
*sign* of its dpKa on binding. If it inverts, reaching more partners will not help.

### 9. Score EVERY pose, never just the first

Structure predictors emit several diffusion samples per job; scoring only `model_0`
throws away the error bar you already paid for. Re-scoring 49 designs across 3 poses
each (147 complexes) gave a median within-design window spread of **3.3x** - as large
as the entire signal range in that campaign. Only 3 of 49 designs were acid-preferring
in ALL three poses, and the window's *direction* flipped between the single-pose value
and the 3-pose geometric mean in **11 of 49** designs (8 false positives, 3 false
negatives). One design reported at 6.56x from a single pose gave 0.14 / 4.54 / 6.31
across its three - a 45.6x spread with one pose strongly inverted.

Single-pose estimates are not biased (median log-ratio +0.004) - they scatter ~2x in
either direction. So this is not a correction you can apply after the fact; you need
the replicates. Report the geometric mean and the full range, require
`all poses > 1` as a reproducibility gate before any candidate is called selective, and
prefer a design with a tight spread over one with a higher single-pose maximum.

## Install

The design checkpoint ships in the repo (~21 MB,
`checkpoints/potts_v6_afdb_edge_his0.3_acid0.06/epoch-0125.ckpt`) — no weight
download needed. RF3 folding weights do **not** ship; fold with whatever you have.

`install.sh` assumes Linux + `uv`. On a wheel-only or Windows sandbox, build an
equivalent env (Python 3.12) and put the packaged source trees on `PYTHONPATH`
instead of an editable install:

```
PYTHONPATH = <repo>/foundry/models/mpnn/src : <repo>/foundry/src : <your shims>
```

Dependencies: `atomworks rootutils environs jaxtyping beartype loralib einops einx
opt_einsum dm-tree zstandard toolz threadpoolctl propka rich lightning ipdb flaml
xgboost lightgbm`.

Version pins that matter (wrong ones fail at checkpoint-unpickle or model-load):

- `omegaconf==2.3.0` + `antlr4-python3-runtime==4.9.3` + `hydra-core==1.3.2`, all
  from **pip**. omegaconf is needed to unpickle the checkpoint even though the
  source never imports it. Conda-forge builds of omegaconf/hydra mismatched here
  (`ImportError: cannot import name 'ConfigIndexError'`, then an 0.x hydra).
- `scikit-learn==1.8.0` — matches the bundled FLAML labeller pickles.
- `assertpy` is **sdist-only**; a wheel-only pip cannot install it. It is used in a
  handful of foundry-core files and none on the design path, so a small shim module
  is enough (`assert_that` returning a fluent object whose checks raise
  `AssertionError`).

The design engine runs on **CPU** — set `CUDA_VISIBLE_DEVICES="-1"` and keep the GPU
for folding. Use the literal `-1`, not an empty string: on Windows an empty value is
not honoured, `torch.cuda.is_available()` stays True while `device_count()` is 0, and
a checkpoint load then fails with *"Attempting to deserialize object on CUDA device 0
but torch.cuda.device_count() is 0"*.

### HBPLUS

The v6 *edge* checkpoint computes H-bond edge features by shelling out to **HBPLUS**,
a Unix C binary; there is no fallback path. Point `HBPLUS_PATH` at the executable and
keep `standard.data` and `vdw.radii` beside it.

The EBI thornton-srv `hbplus.tar.gz` URL now returns an HTML landing page rather than
the source. The page redirects to `github.com/RomanLas/HBPLUS`, which carries both
the real source tarball and **`winhbplus.zip` — a prebuilt Windows `hbplus.exe`**
(v3.2). That binary runs natively, which is the way out when no compiler is
available.

Two assembly details that each cost a debug cycle:

- **`winhbplus.zip` extracts one level deeper than expected** — the binary lands at
  `<dest>/hbplus/hbplus.exe`, so a `HBPLUS_PATH` pointed at `<dest>/hbplus.exe`
  fails with `FileNotFoundError: [WinError 2]` from inside `subprocess.run`, which
  reads as a code bug rather than a path bug. Locate it with a recursive search.
- **`standard.data` and `vdw.radii` are NOT in the Windows zip** — extract them from
  the `hbplus.tar.gz` in the same repo and place them beside the executable. The
  caller runs HBPLUS with `cwd` set to a temp directory, so also set **`HBDIR`** to
  the executable's directory or it will not find them.

### Sandboxes that block multiprocessing

The FLAML protonation labeller calls `joblib.Parallel`. Where the sandbox denies pipe
creation (`PermissionError: [WinError 5]` from `CreateNamedPipe`), wrap every engine
call:

```python
import joblib
with joblib.parallel_config(backend="sequential"):
    engine = PottsMPNNPHEngine(...)
    ...
```

## Usage

```python
import os, joblib
from dataclasses import replace
from biotite.structure.io.pdb import PDBFile
from mpnn.inference_engines.potts_mpnn_ph import PottsMPNNPHEngine, PHDesignCriteria

with joblib.parallel_config(backend="sequential"):
    engine = PottsMPNNPHEngine(checkpoint_path=CKPT, extended_vocab="v6",
                               out_directory=None, write_fasta=False,
                               write_structures=False)
    aa = PDBFile.read("complex.pdb").get_structure(model=1)   # target + binder
    crit = PHDesignCriteria(
        method="block_descent", backend="potts",
        center_types=["HIS-P", "HIS-P", "HIS-P"],   # EXACT multiset -> 3 His centres
        center_count=3,
        placement_region=["interface"],             # centres must gate binding
        placement_by="scan_potts",
        dep_map={"HIS-P": ["HIS-S"], "ASP-P": ["ASP-D"], "GLU-P": ["GLU-D"]},
        forbidden_tokens=["HIS-A", "ASP-A", "GLU-A", "UNK"],
        seed_source="native", combined_lambda=0.85,
        temperature=0.05, block_size=3, neighbour_k=16, max_mutations=24,
    )
    # sweep lambda to trace the stability / selectivity trade-off
    crits = [replace(crit, combined_lambda=l) for l in (0.55, 0.7, 0.85, 1.0)]
    designs = engine.run_ph_redesign(atom_array=aa, binder_chain="B",
                                     criteria_list=crits, seed=0, n_jobs=1)

for d in designs:
    print(d.combined_lambda, d.final_potts_energy, d.selective_energy,
          d.center_res_ids, d.center_protonation_types, d.canonical_sequence)
```

### Criteria fields worth knowing

- `center_types` — **exact protonation-type multiset** to place (e.g.
  `["ASP-P","ASP-P","HIS-P"]`); placement still chooses *where*. Produces one plan.
  This is the field that controls how many live sites you get.
- `center_count` — number of centres when using the pool/combo logic.
- `explicit_centers` — `[{res_id, protonation_type}, ...]` to pin positions by hand;
  its length overrides `center_count`.
- `placement_region` — subset of `{interface, core, surface, all}`.
- `placement_by` — `scan_potts` (rank by Potts dE_p - dE_d), `scan_mpnn`, `random`.
- `combined_lambda` — 0 = stability, 1 = selectivity.
- `neighbour_k`, `max_mutations` — size of the redesigned neighbourhood / mutation budget.
- `seed_source="native"` starts from the input sequence.

### Enumeration cost - the setting that decides your wall clock

The defaults `n_plan_samples=64`, `max_plans_per_seed=8` run up to **8 full solves per
criterion**, so an 8-criteria sweep is 64 solves. Nothing is printed between them, which
makes a slow run indistinguishable from a hang. Set `max_plans_per_seed=1` (and
`n_plan_samples=16`) for an interpretable sweep, and run one criterion per
`run_ph_redesign` call so you get a progress line per design.

Cost is **highly machine-dependent**. The original note of ~75 s per design for a
~320-residue complex did not reproduce on a 24-core Windows CPU, where a single solve
took **~15 min** (904 s for a 327-residue complex, 1245 s for a 452-residue one). Budget
from a measured single solve before launching a sweep, not from this number. A 16-solve
sweep is a several-hour job on that hardware, and folding no longer dominates the wall
clock - design does.

### Running it on Windows

Four host-level traps, each of which cost a full debug cycle:

- **Sanitise `PATH` before launching the env's interpreter.** Invoking one conda env's
  `python.exe` while another env is on `PATH` fails with a delay-load DLL error
  (`0xC06D007F`). Set `PATH` to the target env's own directories plus `C:\Windows\system32`.
- **A bare `python` may not be your env.** On one host it resolved to an unrelated
  system Python 2.7 install. Always call the env's `python.exe` by absolute path.
- **PowerShell pipelines buffer child output until the process exits**, so progress
  polling shows nothing and a working run looks hung. Redirect to a file
  (`> log.txt 2>&1`) and poll the file. Do not kill a quiet run on the assumption it has
  hung - check CPU occupancy first.
- **Detached children are reaped.** `Start-Process ... -NoNewWindow` children are killed
  when the launching shell exits; keep the process in the foreground of a backgrounded
  cell instead.

The co-folding side has its own: pass `--num_workers 0` (the sandbox denies
`CreateNamedPipe`, which the dataloader's worker pool needs) and `--no_kernels` (the
fast triangle-attention kernels ship Linux-only binaries). Keep run and output names
short - the Windows 260-character path limit is reachable with nested
`out_<job>/boltz_results_<job>/predictions/<job>/<job>_model_N.pdb`.

## Helpers

`kernel.py` provides, in the python kernel:

- `ph_window_ceiling(pH_lo, pH_hi, max_sites)` — the sites-vs-window physics table.
- `site_linkage(pka_free, pka_bound, pH_lo, pH_hi)` — per-site linkage factor
  (`>1` acid-preferring, `~1` inert, `<1` inverted).
- `combined_window(sites)` — product of per-site linkages, plus ddG.
- `sites_needed(window, pH_lo, pH_hi)` — how many ideal sites a target window implies.
