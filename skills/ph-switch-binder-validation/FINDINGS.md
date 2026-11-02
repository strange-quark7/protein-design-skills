# Campaign findings — the evidence behind this skill

Every quantitative claim in `SKILL.md` traces to this document. The source was a
four-cycle campaign to design a pH-selective binder against human EGFR domain III,
targeting binding at pH 6.5 with no detectable binding at pH 7.4 (≥100-fold window),
plus mouse cross-reactivity.

**Outcome: the specification was not met.** Reporting that plainly is the point — the
transferable value is in *why*, and in the measurements that established it.

## Scale of the campaign

| | |
|---|---|
| design cycles | 4 |
| designs built and folded | 56 |
| predicted structures (3 diffusion poses each) | 168 |
| complexes fully scored (epitope + conservation + per-site pKa) | 147 |
| titratable sites measured | >500 |
| GPU time (single 8 GB laptop GPU) | ≈11 h |
| designs meeting the ≥100× specification | **0** |
| designs reproducibly acid-preferring across all poses | **3 / 49** |
| …and engaging the intended epitope | **1 / 49** |
| …and also developability-clean | **0 / 49** |

## Finding 1 — The physics ceiling

Over pH 6.5 → 7.4 (0.9 units), one ideal titratable site contributes at most
10^0.9 = **7.9×**. Therefore 1 site → 7.9×, 2 → 63×, 3 → 501×. A 100× window needs
**2.22 ideal sites**, i.e. three real histidines near their theoretical limit.

Best achieved in the campaign: **1.25 ideal-site equivalents**.

## Finding 2 — Epitope drift manufactures false positives

Cycle 1 (manual histidine–carboxylate design, n = 7): predicted window correlated
**r = −0.84** with on-target epitope contact. Every complex genuinely engaging the
intended face was pH-inverted; the only two above 1.0× had **0 %** on-target contact.
Both were false positives.

After restricting placement to interface positions and pinning composition explicitly
(cycle 2, n = 40) the correlation became **+0.17**, and 8 designs were both
epitope-valid and acid-preferring. The confound is fixable, but only if you measure it.

## Finding 3 — Two intuitive mechanisms, two distinct failures

Cycle 1 tested (A) binder histidines paired to conserved target carboxylates, and
(B) binder carboxylates paired to the target's own His433. Windows: 0.304×, 1.180×,
0.538× — all inverted or off-epitope.

- **Mechanism B was anti-selective by construction.** The target histidine *loses*
  ≈2 pKa units on binding, biasing it toward neutral in the complex, so a switch
  anchored to it selects for alkaline pH.
- **Mechanism A sites went inert.** Engineered carboxylates settled into contexts with
  pKa far below the operating window. Two engineered positions had ΔSASA = 0.0 — they
  never contacted the interface at all.

Across all measured sites, ΔpKa on binding correlated **+0.71 with retained solvation**
and **−0.74 with burial**: this interface stabilises anions and destabilises protonated
histidine.

## Finding 4 — Placement scans choose thermodynamically poor positions

The protonation-aware design model's placement scan selected only four binder positions
across all 40 designs — 101, 102, 106, 108 — with live-site rates of **13–20 %** and a
*median* linkage below 1.0 at every one.

Measured on the lead complex:

| position | distance to nearest carboxylate | SASA when bound | live rate |
|---|---|---|---|
| 101 | 15.9 Å | 57.5 Å² | 13 % |
| 102 | 11.1 Å | 111.2 Å² | 18 % |
| 106 | **4.6 Å (D379)** | 30.4 Å² | 15 % |
| 108 | 8.8 Å | **0.0 Å² (buried)** | 20 % |

Three of four had no electrostatic partner. The model's selectivity-gap objective does
not require one.

## Finding 5 — More pinned sites does not help

Compositions swept at 2His / 3His / 4His / 2His+Asp, five trade-off settings, two seeds:

| composition | max window | median window | mean live sites | mean inverted sites |
|---|---|---|---|---|
| 2His | **13.32×** | 0.88× | 0.8 | 2.5 |
| 3His | 6.56× | 0.69× | 0.9 | 2.6 |
| 4His | 6.47× | 0.67× | 0.8 | 2.3 |
| 2His+Asp | 1.10× | 0.48× | 0.5 | 2.5 |

Correlation between requested centre count and log window: **+0.02** — no relationship.
Maximum live sites in any design: **3**. Designs with **zero** live sites: **18 / 40**.
Interface carboxylates stayed inert throughout, confirming Finding 3.

## Finding 6 — Per-site linkage does not transfer between variants

The most expensive lesson. Top designs' *inverted* sites cost a mean factor of
**2.36×**; for the best design, live sites alone multiplied to **61.9×**. The obvious
intervention — keep the live histidines, neutralise the inverted ones — was built as six
surgical repairs plus three model-driven variants.

| family | predicted | measured |
|---|---|---|
| keep His102 + His108 | 61.9× | 0.45×, 0.71×, 1.48× |
| keep His101 | 11.6× | 0.30×, 0.49×, 1.12× |
| model-driven, pinned live sites | — | 0.43×, 0.58×, 0.58× |

Every retained site lost its linkage (His102: 4.88 → 0.99 / 2.04 / 0.61; His101:
2.07 → ≈0.65). Decisively, **target His433 — on the target protein, never mutated —
fell from 3.91× to 0.69–0.97×** after two binder mutations >15 Å away.

The window is an emergent property of the whole bound interface. The
product-of-kept-linkages arithmetic is invalid.

## Finding 7 — The noise floor, measured

Target His433 is chemically identical in all 46 complexes — an internal precision
control.

| quantity | range | SD |
|---|---|---|
| free-state pKa (binder absent) | 6.16 – 6.38 | **0.05** |
| bound-state pKa | 2.58 – 10.25 | **1.47** |

The free state is reproducible — the pKa predictor is self-consistent given the same
structure. The bound state varies **29× more**, driven entirely by predicted-pose
variation. Propagating the 10–90 % bound range gives an apparent single-site window of
**0.69 – 2.63×** from method noise alone.

Against that band: best window anywhere **13.32×** (≈5× above the top), best
developability-clean candidate **6.56×** (≈2.5×), and all nine cycle-3 designs
(0.30–1.48×) **inside** it.

## Finding 8 — Single-pose scoring flips verdicts

All analysis initially used only the first of three diffusion poses per design.
Re-scoring all three (49 designs, 147 complexes) cost 15 min of CPU and no new GPU time:

- median within-design window spread: **3.3×**;
- designs acid-preferring in **all** three poses: **3 / 49**;
- window **direction** flips between single-pose and three-pose mean in **11 / 49**
  designs (8 false positives, 3 false negatives);
- single-pose values are **not biased** (median log-ratio +0.004) but scatter ≈**2×**
  either way.

This forced a retraction: a design reported at 6.56× and recommended for synthesis gave
**0.14× / 4.54× / 6.31×** across its poses — a **45.6×** spread with one pose strongly
inverted, geometric mean 1.58×.

## Finding 9 — The design objective does not predict the validated window

The model's internal `selective_energy` ranked designs in essentially the wrong order:
its best-scoring design (−15.1) had the **worst** measured window (0.58×), while a
mid-scoring one (−12.8) gave the best. Use it to search, never to select.

## Finding 10 — Designed interfaces drift hydrophobic

Interface descriptors against a folded natural antibody control for the same target:

| | natural control | candidate A | candidate B |
|---|---|---|---|
| buried surface area per side | 792 Å² | 797 Å² | 758 Å² |
| interface contacts (≤4.5 Å) | 249 | 202 | 177 |
| interface hydrogen bonds | 15 | 9 | 12 |
| interface salt bridges | 8 | 4 | 4 |
| apolar fraction of buried area | **0.449** | **0.645** | **0.755** |

Buried area is normal, but both designs make fewer polar contacts and bury a markedly
more hydrophobic interface — a recognised non-specific-binding and aggregation
signature. Absolute values look acceptable; the comparison is what reveals the problem.

## Finding 11 — Developability must gate ahead of the objective

The highest-window design of the entire campaign (13.32×) carried an **unpaired
cysteine** and **14 methionines** against one in the parent scaffold. Ranking on the
design objective alone would have promoted an unmanufacturable molecule.

This creates a trap: repairing a liability changes the sequence, which per Finding 6
invalidates the measured window. Repair, then re-fold and re-gate from scratch.

## Finding 12 — Geometry can cap live sites below the specification

Searching the target surface for orthologue-conserved, unpaired, solvent-exposed
carboxylate triads within one paratope span returned a clear best candidate:
**D379 + D388 + E391** — all mouse-identical, span 18.9 Å, combined SASA 154 Å², patch
centroid only 12.6 Å from the existing epitope.

But in the existing binding mode, **D388 and E391 sit 18.7 Å and 18.8 Å from the nearest
CDR position** — beyond side-chain reach — while seven CDR positions lie within 12 Å of
D379. **That pose supports exactly one live site**, capping the achievable window at
7.9× regardless of sequence. Reaching the triad is a re-docking problem; an
inverse-folding model cannot reposition a binder.

## Finding 13 — Check numbering and construct truncation

**Numbering.** The reference crystal structure used mature-protein numbering, offset
**24 residues** from database numbering. Read directly, only **50 of 609** residues
matched; at the +24 offset, **607 of 609** matched. The first pass of an epitope
occlusion analysis compared entirely wrong residues and had to be redone.

**Truncation.** Folding used an isolated 205-residue domain while the affinity objective
specified the full 621-residue ectodomain. Comparing epitope SASA between the isolated
domain and the intact structure showed only **3.2 %** occluded (1568 → 1518 Å²), with the
switch-relevant residues all exposed (D379 35.3 Å², D347 136.3 Å², H433 145.1 Å²). The
fragment was a fair proxy — but that was a *result*, not an assumption.

## Scope limits of this evidence

Honesty about what the campaign did *not* establish:

- **No classical protein–protein docking.** Poses came from co-folding; no rigid-body
  sampling or docking scoring function. No docking software was available.
- **No binding affinity was ever computed.** Fold-confidence metrics (ipTM, pLDDT) are
  not affinity. The structure predictor's affinity module was small-molecule only
  (≤128 atoms) and could not score a 122-residue binder. The affinity objective was
  therefore never scored.
- **No molecular dynamics.** All pKa values are single-snapshot calculations on static
  predicted structures. Constant-pH MD — the method that would test whether the large
  bound-state shifts are real — was not run.
- **No experimental data.** Every number here is a prediction.
- **Not assessed:** immunogenicity, expression yield, thermal stability, aggregation
  beyond the sequence and hydrophobicity proxies, and off-target binding against
  homologous receptors.

One consequence is worth stating: because the noise floor (0.69–2.63×) sits within a
small factor of the best reproducible signal (2–5×), **the per-design rankings in this
campaign are substantially noise-driven**. The generalisable findings above rest on
aggregate statistics and controlled comparisons, which are far more robust than any
individual design's score — but no single predicted window from this work should be
treated as a measurement.

---

# Cycle 5 — removing the geometry limitation

Cycles 1–4 concluded that the binding mode was the binding constraint: the paratope
reached one of the three carboxylates a ≥100× window needs. Cycle 5 removed that
limitation and re-tested. **The switch still failed**, which is what makes it the most
informative cycle of the campaign.

## Finding 10 — A single-domain paratope cannot span a three-carboxylate cluster

The target cluster (three conserved, unpaired, solvent-exposed carboxylates; SASA
45.0 / 84.6 / 59.1 Å²) is essentially co-facial — depth spread along the surface normal
only **2.5 Å** — but spans **18.3 Å** between the outer two carboxylate groups.

Pocket-constrained co-folding (`pocket` constraint, `force: true`, 6 Å) was used to
place a binder over the cluster. Distances are each carboxylate's Cα to the nearest
binder Cα; ≤11 Å is side-chain reach:

| binding mode | format | d1 | d2 | d3 | within reach | ipTM |
|---|---|---|---|---|---|---|
| original mode | VHH 122 aa | 6.5 | 16.3 | 19.5 | **1** | 0.84 |
| constrained | VHH 122 aa | 13.4 | 10.7 | 9.4 | **2** | 0.73 |
| constrained | scFv 247 aa | 6.1 | 5.4 | 9.4 | **3** | 0.28–0.55 |

Constrained onto the cluster, the VHH **slides across** rather than straddling it —
trading the first carboxylate for the other two. Two sites, never three. Reaching all
three required changing format to a 247-residue scFv with a six-CDR paratope.

**Transferable:** when a specification needs *n* partnered sites, check the paratope
span against the cluster span before choosing a format. Roughly, a single-domain VHH
handles a cluster up to ~12–13 Å; beyond that an scFv or Fab is required. This is a
format decision, and it is cheaper to make before designing than after.

## Finding 11 — Geometry is necessary but not sufficient (the key result)

With all three partners in reach, histidines were placed directly at the binder
positions adjacent to each carboxylate and measured (folded unconstrained, three poses
each, PROPKA free vs bound).

| variant | switch sites | window (geo. mean) | range | all poses > 1? |
|---|---|---|---|---|
| scFv-WT (unmutated control) | — | 0.799 | 0.688–0.870 | no |
| scFv-2His | 2 | 0.944 | 0.629–1.859 | no |
| VHH-2His | 2 | 0.636 | 0.632–0.641 | no |
| scFv-3His | 3 | **0.586** | 0.341–1.256 | no |

**Of 21 switch-site measurements: 0 live, 3 inert, 18 inverted.** Adding histidines made
the binder *more* alkaline-preferring, not less — the three-histidine variant is below
its own unmutated control.

The mechanism is in the pKa shifts; every engineered histidine lost pKa on binding:

| variant · position | free pKa | bound pKa | linkage |
|---|---|---|---|
| scFv-3His · 103 | 6.30 | 5.31 | 0.61–0.91 |
| scFv-3His · 233 | 6.71 | 4.16 | 0.43–0.48 |
| scFv-3His · 234 | 3.99 | −0.36 | 1.00 (inert) |
| scFv-2His · 103 | 6.25 | 5.61 | 0.71–0.92 |
| VHH-2His · 44 | 6.23 | 5.07 | 0.71–0.73 |
| VHH-2His · 45 | 5.69 | 4.00 | 0.88–0.89 |

Burial desolvates the imidazole faster than the adjacent carboxylate stabilises its
protonated form, so the bound state favours neutral histidine and therefore the higher
pH. Position 234 shows the degenerate case: free pKa already 3.99, far below the window,
driven to −0.36 on binding, contributing nothing.

## Finding 12 — This result is not noise

Noise floor measured the same way as Finding 7, on the unmutated target histidine across
12 complexes:

| | range | SD |
|---|---|---|
| free-state pKa | 6.22 – 6.26 | **0.01** |
| bound-state pKa | 3.84 – 6.23 | 0.87 |

Implied single-site band from noise alone: **0.69–0.98×, only 1.43-fold wide** — much
tighter than the 0.69–2.63× of cycles 1–4, because a locally-built MSA replaced the
public MSA server and gave more consistent poses.

Against that band, per-site linkages repeat to two decimal places across independent
poses: 0.71 / 0.73 / 0.72 at one position, 0.89 / 0.88 / 0.89 at another. **A tight
spread is usually reassuring; here it confirms a real and unwanted result.**

## Finding 13 — Affinity was adequate throughout; selectivity was always the limit

A contact-based affinity predictor was run on all 33 complexes, calibrated against a
natural antibody for the same target (−11.38 kcal/mol, Kd 5.4 nM). Every candidate fell
in the control's range (−10.4 to −12.3 kcal/mol, 0.5–20 nM), while un-redesigned
scaffolds forced onto the new site scored weakest (184–542 nM) — the expected ordering,
and a check on the method.

**Transferable:** score affinity early even if the specification is about selectivity. A
pH window is a *ratio* and says nothing about absolute binding; confirming affinity is
adequate tells you which objective is actually limiting, and it is cheap.

## Finding 14 — Direct orthologue co-folding separates candidates that sequence identity does not

Cross-reactivity had been scored as epitope sequence identity. Co-folding both final
candidates against the mouse orthologue (87.8 % identical construct), three poses each:

| | human ipTM | mouse ipTM | epitope overlap (Jaccard) |
|---|---|---|---|
| candidate A | 0.82 | 0.69 | **0.80** |
| candidate B | 0.49 | 0.60 | 0.41 |

Candidate B engages the orthologue but at a **shifted epitope** — its binding is a
different interaction, not the designed one. Sequence identity of the predicted epitope
(80.7 % vs 77.7 %) had ranked these two as near-equivalent.

**Transferable:** if orthologue cross-reactivity is an objective, fold against the
orthologue. Epitope-identity is a screen, not a measurement, and it cannot see a binder
that relocates.

## Cycle 5 scope limits

- The scFv arm never completed a full inverse-folding design cycle; its switch variants
  were built by direct histidine placement (see the `protonpottsmpnn` skill, placement
  section, for why).
- Mouse cross-reactivity was measured for the two incumbents only, not the new variants.
- All structures aglycosyl; all complexes used an isolated domain, not the full
  ectodomain; no molecular dynamics; no experimental data.
