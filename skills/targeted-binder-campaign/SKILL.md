---
name: targeted-binder-campaign
description: "Run a de novo protein binder campaign against a defined epitope under constraints - pH-gating, cross-species reactivity, a required contact, developability, or any multi-objective brief. Covers backbone generation with hotspot conditioning, sequence design settings that avoid the composition collapse, and the co-fold gate stack (interface confidence against a control, epitope fidelity, pose resampling) that decides what is real. Use this skill whenever the user wants to design a binder, nanobody, minibinder or scFv against a named target or epitope, wants to add a constraint to an existing binder (pH switch, species cross-reactivity, a specific contact), asks which of several designs to advance, or asks why a design campaign is not converging - even if they do not name a specific tool. The central lesson it encodes: almost every wrong conclusion in such a campaign comes from measuring an interface property on a structure that has no side chains."
---

# Targeted binder campaign

A constrained binder campaign is three generative stages feeding a filter, and
nearly all of its failures are measurement failures rather than design failures.
This skill is about where the measurements go wrong and which gate catches each
one.

The pipeline: **backbone generation** (RFdiffusion or similar, conditioned on
target + hotspots) -> **sequence design** (ProteinMPNN family) -> **co-folding**
(Boltz-2, Chai-1, OpenFold3) -> **selection**. Sibling skills cover the tools:
`proteinmpnn`, `ligandmpnn`, `boltz`, `chai1`. For pH-switch physics and its
nine failure modes, load `ph-switch-binder-validation`; for a protonation-aware
design model, `protonpottsmpnn`.

## The error that will cost you the campaign

**Generated backbones have no side chains.** RFdiffusion emits N, CA, C, O only
— 4 atoms per residue — and designed positions are poly-glycine until sequence
design. A side chain reaches ~4 A past CB, so any contact, burial, or
reachability measured on a backbone is wrong by roughly that much on *both*
sides of an interface.

This single mistake produced three confident, wrong conclusions in one campaign:

| measured on backbones | conclusion drawn | truth on folded structures |
|---|---|---|
| paratope burial by a designed lid | "0/50 designs mask the target" | burial numbers were mostly missing side chains |
| contact with a non-conserved hotspot | "5 designs avoid it, use them for cross-reactivity" | **all 18** engaged it at 1.7-3.5 A |
| distance to the nearest target carboxylate | "a His-carboxylate switch is geometrically impossible here" | 15/18 designs already within 4.0 A |

That third one nearly killed the project's best mechanism.

**The rule:** any quantity that depends on side-chain atoms — contacts, burial,
SASA, salt bridges, pocket composition, hydrogen bonds — is measured on a
**folded, full-atom** structure, never on a generated backbone. Use backbones
only for what they can support: chain length, fold topology, CA-level geometry,
and centroid distances.

If you must compare a backbone to a full-atom reference, strip the reference to
backbone atoms first so both sides match. A proxy computed from CB positions
alone does not rescue this: tested against measured pKa it gave rho = -0.09
(p = 0.51) where the full-atom measurement gave rho = -0.40 (p = 0.003).

## Stage 1 - backbones

Use binder/PPI mode with explicit hotspots. Without hotspot conditioning the
generator packs the binder against whatever surface is convenient: one campaign
got 32% of contacts on-epitope unconditioned versus 84% with hotspots.

**Hotspots must fit inside one paratope.** A 110-130 aa binder presents roughly
a 20-25 A face. Check the pairwise geometry of your chosen hotspots *before*
committing to a long run — a set spread across a 32 A epitope cannot all be
engaged, and the generator will silently drop the outlier. In one case a hotspot
25-31 A from the other three was missed by every design; the remaining three sat
within 11.7 A and were engaged reliably. Replacing it with a residue 8.8 A from
the cluster centroid fixed the run.

**Run a 2-5 design canary first.** It costs minutes and catches contig errors,
wrong chain layout, the wrong checkpoint being selected, and unreachable
hotspots. Confirm from the canary: chain count and lengths, which checkpoint
loaded, and on-epitope contact fraction against the random expectation
(epitope size / target size).

**Generation converges — screen before scaling.** Doubling a batch from 42 to 84
backbones produced nothing better on any axis; the all-hotspots-engaged rate
actually fell from 40% to 26%, and the best new backbone only tied one already
in hand. Score a tranche and compare it to the earlier tranche before spending
GPU-hours on the rest. Stop when a new tranche stops producing leaders.

## Stage 2 - sequence design

Two settings decide whether the sequences are usable at all.

**Composition collapse.** On idealised backbones (noise scale 0, very regular
geometry), default ProteinMPNN weights at low temperature return poly-alanine:
one run gave **38% Ala and 58% hydrophobic** against 5% and 36% in a natural
nanobody reference. Those sequences would not express. The fix is soluble
weights, temperature ~0.25, Cys omitted, and a global alanine bias around -1.5
— which brought the same backbones to 8.6% Ala and 34.2% hydrophobic.

Watch the side effect: the alanine penalty displaces into glutamate. Charged
fraction went to 42% against 24% in the reference. Select low-charge designs out
of the pool rather than re-running; regenerating trades one composition artefact
for another.

**Always compare composition to a reference set of real, expressible proteins of
the same class** — Ala, Gly, hydrophobic and charged fractions, plus distinct
residue count. The collapse is invisible by eye and obvious against a reference.

**Per-residue bias needs to be strong.** To install a residue at a chosen
position, a bias of 2.5 achieved only 17% uptake; **6.0 achieved 239/240**. For
a single critical position (a required contact) use ~8.0. Verify uptake in the
output sequences rather than assuming the bias worked.

**Design on the best structure you have.** Once you have folded complexes,
redesign on those rather than on the original backbones — the folded pose is the
validated one and it carries side chains.

## Stage 3 - the gate stack

Run these in order. Each catches something the others cannot.

### Gate 1: interface confidence, against a control

ipTM > 0.5 is the generic pass line, but it is far more informative to fold a
**known binder for the same target** and read designs against that. One campaign
had a wild-type control at ipTM 0.891, which turned "above 0.5" into "5 of 18
designs match or beat a real antibody."

### Gate 2: epitope fidelity on the folded complex

ipTM says two chains form a confident interface, **not that it is the designed
one**. A design once scored 0.889 with its partner chain 45 A off target.
Measure the interface the prediction actually made: on-epitope contact fraction
and hotspot engagement, in the fold's own frame.

Expect the footprint to expand beyond a narrow epitope definition. Judge whether
off-epitope contacts are *contiguous* with the epitope (flanking residues =
footprint growth, fine) or elsewhere on the target (relocation, fail).

**Validate this gate on a known-answer control BEFORE reading any design number
through it.** Fold a natural binder whose epitope is crystallographically known
and confirm your code recovers it. Running the control alongside the designs
rather than ahead of them is not the same thing: in one campaign a bad number was
reported from the design set before the control's own 2.9% score revealed the
scoring was broken, and the claim had to be retracted.

#### Homo-oligomeric targets: score every equivalent interface

If the target is a homodimer, homotrimer or higher oligomer, the epitope recurs
at each equivalent subunit interface, and a correct binder may occupy **any** copy
of it. An epitope defined on the A/B pair is satisfied just as well on B/C or C/A,
but literal chain-label matching scores those as total misses.

This is not a corner case — it inverted a campaign's headline statistic:

| | literal A/B match | symmetry-aware |
|---|---|---|
| clinical-antibody control recovery | **0.029** | **0.912** |
| designs passing gate 2 | 1/10 | 7/10 |
| ipTM vs epitope recovery | rho = -0.57 | rho = +0.10 |

The `-0.57` was reported as a real anti-correlation before the control exposed it
as an artefact. **The control is what catches this**: a clinical antibody scoring
2.9% on its own published epitope indicts the scoring code, not the prediction.

Score the epitope over every subunit permutation and keep the best-fitting
assignment, recording which interface the binder actually occupied:

```python
best = None
for first, second in [("A","B"),("B","C"),("C","A"),("B","A"),("C","B"),("A","C")]:
    mapped = remap_epitope_chains(epitope, {"A": first, "B": second})
    ef = epitope_fidelity(pdb, mapped, mapped_hotspots, ...)
    key = (ef["on_epitope_frac"], ef["hotspots_engaged"])
    if best is None or key > best[0]:
        best = (key, ef, f"{first}/{second}")
```

**Do not truncate the oligomer to save compute.** Deleting a protomer exposes
faces that are buried in the real assembly, and adjacent exposed half-sites can
fuse into a groove that does not exist. In the same campaign a two-protomer model
of a trimer created exactly such an artificial groove, and the co-folder docked a
clinical antibody into it preferentially at ipTM 0.77-0.82 with 0-20% recovery of
the true epitope. Model the full oligomer; the extra residues cost fold time, not
correctness.

### Gate 3: pose resampling — what separates real from lucky

This is the gate most campaigns skip, and it invalidated half of one campaign's
headline claims. Re-fold the designs that matter with many independent diffusion
samples (25 is enough) and ask how often the feature you care about is present.

For a required salt bridge across four designs at 25 poses each:

| design | first prediction | across 25 poses |
|---|---|---|
| A | 2.83 A | **25/25**, spread 0.28 A |
| B | 3.12 A | **25/25**, spread 0.84 A |
| C | 2.71 A | 17/25, outliers to 35.6 A |
| D | 2.84 A | **7/25**, outliers to 21.9 A |

**The single-prediction distance did not predict robustness** (rho = +0.21).
Design D looked second-best on one structure and holds the feature in barely a
quarter of poses — in the rest the whole binder adopts a different pose. Any
claim resting on one predicted structure would have been wrong.

Report occupancy ("2.78 A in 25/25 poses"), not a single distance. A converged
feature shows a tight spread; a fragile one shows bimodal distances with far
outliers.

### Gate 4: developability, as a hard gate

Apply before ranking, not as a tradeable term: no unpaired Cys (simplest to omit
Cys entirely), no N-glycosylation sequons, bounded Met, Gly below ~15%,
hydrophobic fraction 0.25-0.50, enough distinct residues, no long hydrophobic
runs.

## Never rank on geometry when you can afford to fold

Geometric proxies computed from generated backbones were **overturned three
times** in one campaign by the co-fold gate:

- the affinity ranking **inverted** — the proxy's top pick scored 0.683 (last),
  the design it demoted on a charge penalty scored 0.935 (first)
- the category predicted weakest by a composite score delivered 0.838-0.894
- the cross-reactivity criterion was invalid outright (the poly-glycine error)

So: fold **more designs than you intend to ship**, then select on measured
values. Folding 18 to pick 7 is cheap insurance against shipping the wrong one.
Use geometry only to decide *what to fold*, never what to advance.

## Constrained and multi-objective briefs

**Verify the constraint axis actually varies before stratifying on it.** One
campaign built affinity-optimised and cross-reactivity-optimised slots around a
tension that did not exist in folded structures — the epitope was 81% conserved
between species, so every binder to that face was cross-reactive (0.76-0.90 of
contacts conserved, all above a 0.74 wild-type control). The honest presentation
was a panel-wide property with per-design numbers, not a contested category.

When objectives genuinely conflict, prefer **purpose-optimised slots** over one
composite ranking — a design optimised for everything is best at nothing. Give
each slot only the terms it is judged on, and let slots differ deliberately
(a pure-affinity design should carry no switch residue, because an interfacial
titratable residue costs binding energy and that slot has no pH requirement).

**Diversify backbones under uncertainty.** If the designs are not yet validated,
cap designs per backbone (one is cleanest) so a single bad scaffold cannot take
several slots. Report pairwise sequence identity to show the panel is distinct
scaffolds rather than one fold resampled.

## Required-contact constraints (including pH switches)

A common brief is "the binder must contact residue X" — a catalytic residue, an
allosteric site, or a carboxylate that will anchor a titratable switch.

Measure reachability on **folded** complexes. In one campaign the nearest target
carboxylate appeared 14.4 A away on backbones (declared impossible) but was
within 4.0 A in 15 of 18 folded structures.

If designs reach the position but the desired residue is not there, the cause is
usually that sequence design chose a better partner: Arg and Lys outcompete His
for an aspartate, so an unbiased design fills that position with them. The fix is
a per-residue bias at that one position (~8.0), not new backbones. Before
spending GPU on it, check that a position exists whose **CB sits 3.5-7.5 A** from
the target atom — the window where a side chain can reach bonding distance.

For a pH switch specifically, pairing a histidine to a target carboxylate *in the
bound state* raises its pKa on binding, which makes binding acid-preferring. This
is the opposite of burying a titratable residue at an interface, which desolvates
the protonated form, drops its pKa, and runs the switch backwards. Prefer a
**conserved** target carboxylate so the switch costs nothing in cross-species
reactivity. Load `ph-switch-binder-validation` before committing to a pH design.

## What a prediction cannot tell you

Say this plainly in any write-up:

- **ipTM is interface confidence, not affinity.** 0.94 means the model is sure
  these chains form this interface. It is not a Kd.
- **A predicted geometry is not a measured energy.** A salt bridge present in
  25/25 poses establishes that a mechanism is structurally in place; it does not
  give the size of the effect.
- **Contact-based affinity predictors mostly measure interface size.** PRODIGY's
  dG across 34 complexes was 60% explained by raw contact count alone
  (rho = -0.656, R2 = 0.60), it is calibrated on natural complexes rather than de
  novo designs, and it barely tracked ipTM (rho = -0.254, p = 0.15). Report such
  values as secondary with the caveat; do not rank on them.
- **Cross-model consensus is a separate question from pose robustness.**
  Resampling tests one model's convergence. A second, independent co-folder
  (`chai1`) tests whether the finding is model-specific. Say which you did.

## Operational traps

| Symptom | Cause / fix |
|---|---|
| Co-folder exits 0 but no prediction written | Boltz exits 0 on input errors. **Count output files, never trust the exit code.** Implausibly fast runs are the tell — a 300-residue complex cannot fold in 25 s. |
| `Error: '\x00'` from a co-folder | ColabFold null-terminates its a3m. One trailing byte in a valid 900 KB alignment kills every run. Strip nulls after fetching. |
| MSA server fails after ~15 requests | Rate limiting. Fetch the target alignment **once**, cache it, and point every run at the file. A de novo binder has no homologs — single-sequence is correct for that chain, and it also removes the network dependency. |
| Co-folder refuses mixed MSA sources | Some require all chains custom or all auto. Pair a cached target alignment with an explicit empty/single-sequence setting for the binder. |
| Interrupted generation run | Check for an output file lacking its companion metadata file — that design was killed mid-write and must be deleted before resuming, or the skip-existing logic will treat it as complete. Stop runs with a plain signal, never SIGKILL. |
| Design environments fighting over one GPU | Sequence design is CPU-bound; co-folding is GPU-bound. Run them concurrently rather than serially — but **force the design step onto CPU explicitly**. ProteinMPNN grabs CUDA whenever it is visible: left to itself it shares the card with the generator, slowing diffusion 1.3-2.4x and eventually killing it with `CUDA error: unspecified launch failure`. Set `CUDA_VISIBLE_DEVICES=-1` (the value `-1`, not an empty string, which is not honoured on Windows — `torch.cuda.is_available()` stays True while `device_count()` is 0, so the checkpoint then loads onto a device that is not there). |
| Generation crashed partway through a tranche | Make the runner resumable and make it **exit non-zero**. Count completed designs, delete any structure file lacking its companion metadata file (killed mid-write — otherwise the resume count treats a truncated file as finished), and pass the remainder plus a start offset. A runner that swallows the generator's return code and exits 0 will let a chained `cmd1 && cmd2` march on with the first tranche incomplete. |
| Long campaign, weights vanish between stages | Idle sweeps reclaim large files. Multi-GB model checkpoints living in a session workspace are not durable across a long campaign; keep them outside it, or expect to re-fetch. Re-fetching assumes the network still resolves — which is not guaranteed mid-campaign. |

## Reporting

Lead with the gate table: how many designs passed each gate, against which
control. Give per-design measured values, state which quantities are predicted
versus measured, and list what was not validated. When a proxy was overturned by
a gate, say so — it tells the reader which numbers to trust, and a campaign where
no proxy was ever overturned usually means the gates were not run.

## Helpers

`kernel.py` ships the measurements this skill keeps insisting on, all full-atom
and folded-structure based:

- `binder_interface(pdb, target_chain, binder_chain, cutoff)` — per-target-residue
  minimum distance to the binder, plus the contact set
- `epitope_fidelity(pdb, epitope, hotspots, ...)` — on-epitope fraction, hotspots
  engaged, off-epitope residues
- `interface_conservation(pdb, conserved_resids, ...)` — conserved fraction of the
  realised interface
- `contact_occupancy(pose_paths, target_chain, target_resid, binder_chain, ...)` —
  occupancy of a specific contact across resampled poses, with spread
- `reachable_positions(pdb, target_chain, target_resid, binder_chain)` — binder
  positions whose CB lies in the 3.5-7.5 A window for a required contact
- `composition_report(sequences)` — residue fractions against a reference, to catch
  the collapse
- `count_predictions(root, pattern)` — guard against the exit-0-no-output trap
