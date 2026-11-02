---
name: ph-switch-binder-validation
description: "Design and validate pH-selective (acid- or base-switching) protein binders - a binder that engages at tumour pH 6.5 but not blood pH 7.4, or any histidine/carboxylate-gated interface. Reach for this skill to set a realistic selectivity target from first principles, to decide how many titratable sites a specification actually requires, to place switch residues where they can work, and to run the validation gates (epitope fidelity, replicate poses, per-site pKa linkage, interface quality, developability) that separate a real pH switch from prediction noise. Covers the nine failure modes that invalidate naive designs, each with the measurement that catches it, including the one that survives every other check: a geometrically perfect switch site whose interface still destabilises protonated histidine. Tool-agnostic: applies to any inverse-folding or sequence-design model plus any structure predictor. For the specific install and API of a protonation-aware design model, see the protonpottsmpnn skill."
---

# pH-switch binder design and validation

Engineering a binder that works at one pH and not another sounds like a composition
problem - add histidines, they protonate at low pH, done. It is not. This skill is the
distilled result of a five-cycle campaign (68 designs, 201 predicted structures, 159
fully scored complexes) in which **every** intuitive approach failed, each for a
different and measurable reason. `FINDINGS.md` holds the campaign numbers behind every
claim here.

Read Section 1 before designing anything. It will often tell you the specification is
physically out of reach, which is worth knowing before spending GPU-days.

## 0. Which polarity does the brief actually ask for?

Decide this before anything else, because it determines whether the physics is
working for you or against you.

- **Acid-preferring** (bind at low pH, release at high) — the usual tumour-pH or
  endosomal-capture brief. This is what §§1-9 below are about, and it is **hard**:
  the dominant failure mode is that interface burial desolvates the protonated
  histidine, drops its pKa on binding, and runs the switch backwards.
- **Inverted / base-preferring** (bind at the HIGHER pH, release on
  acidification) — release-on-acidification, extracellular retention. For this
  brief **the dominant failure mode above is the design objective**, and the
  catalogue below reads as a recipe rather than a warning list.

A campaign targeting the inverted direction confirmed this directly: every
engineered histidine lost pKa on binding (mean dpKa -1.35), **12 of 15 sites
inverted, 1 inert, 2 acid-preferring**, and **6 of 6 designs held the inverted
direction across all 5 resampled poses** — against a measured noise band of
0.85-1.00x. For contrast, the acid-preferring campaign in §5 obtained 0 of 21 live
sites with 18 of 21 inverted, after a re-dock and a format change specifically
intended to fix it.

So when the brief is inverted, invert the design rules too:

| | acid-preferring brief | inverted brief |
|---|---|---|
| switch-site search | cluster of free carboxylates | conserved, **carboxylate-free**, near a conserved Lys/Arg |
| burial | keep the site solvated at the rim (§3) | **moderate interface burial is the mechanism** |
| His-carboxylate pair in the bound state | the goal | **avoid** — raises bound pKa, flips polarity |
| the physics | fighting you | working for you |

Everything downstream — the validation protocol (§4), the noise floor (§5.4),
pose resampling (§5.2), non-transferable linkage (§5.3) — applies **unchanged** to
both polarities. Only the sign of the selection changes: select `L > 1.05` for
acid-preferring, `L < 0.95` for inverted, and report `1/L` as the achieved window
in the inverted case.

## 1. Set the target from first principles, before designing

A single ionisable site can only do so much. Over a pH interval from `lo` to `hi`, the
best possible affinity ratio from one perfectly-behaved site is `10^(hi - lo)`:

| pH interval | max window per ideal site | sites needed for 100x |
|---|---|---|
| 6.5 -> 7.4 (tumour vs blood) | 7.9x | 2.22 |
| 6.0 -> 7.4 (endosomal vs blood) | 25x | 1.43 |
| 5.5 -> 7.4 (late endosome) | 79x | 1.05 |

"No detectable binding" conventionally means at least 100-fold. So a tumour-pH switch
over 0.9 pH units needs **three histidines each achieving close to its full theoretical
shift** - and "close to full" is the hard part (Section 5). Over the wider endosomal
interval, one or two good sites suffice.

Call `ph_window_ceiling(lo, hi)` and `sites_needed(target, lo, hi)` before committing.
If the arithmetic says four sites, say so to whoever set the specification rather than
designing toward it.

**Consider negotiating the specification.** A graded 5-10x window is a different
engineering problem from 100x, and for many diagnostic or imaging readouts a graded
window is sufficient. That is a question about the end application, and it is far
cheaper to settle it up front.

## 2. Check the target geometry can host the sites

Each switch histidine needs an electrostatic partner - a carboxylate it can salt-bridge
to when protonated. So *n* live sites requires *n* carboxylates within reach of one
paratope. Before designing, search the target surface for that cluster:

- mouse/orthologue-conserved, if cross-reactivity is an objective;
- **not already salt-bridged** to a target lysine or arginine (a paired carboxylate is
  unavailable);
- solvent-exposed (SASA above roughly 20 A^2);
- mutually within one paratope span, about 20 A.

Then check your binder can actually reach them in its existing pose. In the campaign the
best cluster spanned 18.9 A and looked ideal, but two of its three members sat 18.7 and
18.8 A from the nearest CDR position - beyond side-chain reach. That pose could host
exactly **one** live site, and no sequence design could change it. That is a re-docking
problem, not an inverse-folding one.

Measuring this costs minutes and can save a campaign.

**But geometry is necessary, not sufficient — and this is the trap.** A later cycle
removed exactly this limitation and the switch still failed. Pocket-constrained
co-folding put a binder over the full three-carboxylate cluster, and because a
single-domain VHH paratope cannot straddle an 18.3 A span, the format was changed to a
247-residue scFv whose six-CDR surface reaches all three (6.1 / 5.4 / 9.4 A). With every
partner in reach, histidines placed beside them gave **0 of 21 live sites and 18 of 21
inverted**, and the three-histidine variant (window 0.586) was *more* inverted than its
unmutated parent (0.799).

So treat the geometry check as a screen that can only rule a site **out**. Passing it
tells you the site is worth testing, not that a switch will work there. The decisive
property is whether the interface stabilises or destabilises the protonated form once
bound (Section 3, last rule, and failure mode 9).

## 3. Place switch residues where they can work

Three placement rules, each learned from a failure:

**Restrict to the interface.** Design models with a free placement scan will put
protonated centres in the scaffold framework, where they change folding energetics and
cannot gate binding. Restrict placement to positions contacting the target.

**Require an electrostatic partner per centre.** A selectivity-gap objective does not
require one, and in the campaign a model repeatedly chose four positions of which three
had none: two were 11-16 A from the nearest carboxylate, one was fully buried. For every
proposed centre, measure the distance to the nearest conserved unpaired exposed
carboxylate, and reject centres without one.

**Keep the site solvated when bound.** Across all measured sites, the pKa change on
binding correlated **+0.71 with retained solvation** and **-0.74 with burial**. A
histidine buried at the interface core is pushed toward its neutral form and
anti-selects for acid. Aim for interface rim positions that retain roughly 10-25 A^2 when
bound.

**Do not anchor the switch on a target histidine that loses pKa on binding.** Check it
first; in the campaign the obvious target histidine lost about two pKa units on binding,
which selects for alkaline rather than acidic pH - the opposite of the goal.

## 4. The validation protocol

Run these gates **in this order**. Order matters: an early gate failing makes later
numbers meaningless.

1. **Fold confidence.** Predict the complex. Reject obviously failed interfaces.
   Keep every pose the predictor emits - you need them in step 3.
2. **Epitope fidelity.** Does the binder contact the intended face? Require a majority
   of predicted epitope residues to be on-target. **A pH window from a binder that has
   drifted elsewhere is not your designed mechanism** (Section 5).
3. **Replicate agreement.** Score every pose, not just the first. Require all poses to
   agree in direction before calling a design selective (Section 5).
4. **Per-site pKa linkage.** For each titratable site, compute free-state and bound-state
   pKa, and from those its linkage. Identify which sites are actually live.
5. **Orthologue conservation**, if cross-reactivity is an objective: fraction of the
   predicted epitope identical in the orthologue.
6. **Interface quality** against a natural binder control (Section 6).
7. **Developability** as a hard gate (Section 7).

The per-site linkage for one ionisable site, free pKa `pf` and bound pKa `pb`:

    L = [(1 + 10^(pb-lo)) / (1 + 10^(pf-lo))] / [(1 + 10^(pb-hi)) / (1 + 10^(pf-hi))]

`L > 1.05` is live and acid-preferring, `L` near 1.0 is thermodynamically inert, and
`L < 0.95` actively anti-selects. Call `site_linkage()` and `combined_window()`.

**Validate the pipeline on a known complex first.** Fold a natural binder with a
published epitope and confirm you recover it. This calibrates the scale before you
interpret any design number, and it is how you find out your epitope-scoring code has an
off-by-one.

## 5. The nine failure modes

Each of these produced a wrong conclusion in the campaign. Each has a cheap measurement
that catches it.

### 1. Apparent selectivity from epitope drift
Early on, predicted window was **anti**-correlated with on-target epitope contact
(r = -0.84): every complex that genuinely engaged the intended face was pH-inverted, and
the only ones showing a positive window had relocated elsewhere. Every positive was a
false positive.
**Catch it:** gate on epitope fidelity before reading any pH number.

### 2. Single-pose estimates are unreliable
Structure predictors emit several poses; scoring only the first discards the error bar
you already paid for. Median within-design window spread across three poses of the
*same sequence* was **3.3x** - as large as the entire signal range. The window's
*direction* flipped between the single-pose value and the three-pose mean in **11 of 49**
designs. One design reported at 6.56x from one pose gave 0.14 / 4.54 / 6.31 across its
three. Single-pose values are not biased (median log-ratio +0.004) - they scatter about
2x either way, so this is not correctable after the fact.
**Catch it:** score every pose; report geometric mean and range; require all poses above
1.0; prefer a tight spread over a high single-pose maximum. Call `replicate_verdict()`.

### 3. Per-site linkage does not transfer between variants
It is tempting to keep the sites that came out live, mutate the inverted ones away, and
predict the new window as the product of the kept linkages. Tested directly, this failed
by about **40-fold**: predicted 61.9x, measured 0.30-1.48x across six repaired variants.
Every retained site lost its linkage. Most strikingly, a histidine on the **target**
protein - never mutated - fell from 3.91x to 0.69x after two mutations on the binder more
than 15 A away.
**Catch it:** treat the window as an emergent property of the whole bound interface.
Never arithmetic-combine linkages across variants; fold and re-measure each one.

### 4. The noise floor can swamp the signal
Pick a titratable residue identical in every complex - a target-side histidine is ideal -
and measure the spread of its predicted pKa across your whole set. In the campaign, across
46 complexes, its **free**-state pKa was tight (SD 0.05: the pKa predictor is
self-consistent) while its **bound**-state pKa spanned 2.58-10.25 (SD 1.47), driven
entirely by pose variation. That spread alone implies an apparent single-site window of
**0.69-2.63x**. The best signal obtained anywhere was 13.32x, and the best reproducible
signal 2-5x.
**Catch it:** measure this before trusting any ranking. If your best design sits inside
or barely above the band, the ranking is substantially noise and more sampling optimises
against it. Call `noise_floor()`.

### 5. Pinning more sites does not help monotonically
Compositions were swept from two to four histidines. The two-histidine composition
produced both the **best** window and a median of 0.88x; adding a carboxylate capped the
family at 1.10x. Across all designs, correlation between requested centre count and log
window was **+0.02** - no relationship. At most three sites ever came alive, and 18 of 40
designs had **zero**.
**Catch it:** select on *measured* live-site count, not requested composition.

### 6. A design model's internal objective is a search heuristic, not a selection criterion
The protonation-aware model's own selectivity score did not predict the validated window:
its best-scoring design had the **worst** measured window.
**Catch it:** never rank final candidates by the design objective. Independent
validation is mandatory.

### 7. Fold confidence is not affinity
ipTM, pLDDT and similar are *model-confidence* metrics. A confidently predicted complex
is not a tightly bound one. A pH *window* is a ratio of association constants and says
nothing about absolute affinity - a binder with a perfect window and negligible affinity
scores identically.
**Catch it:** state plainly that affinity is unscored unless you have actually computed
it. Note that some structure predictors' affinity modules are small-molecule only and
cannot score a protein-protein interface.

### 8. Designed interfaces drift hydrophobic
Both final candidates buried a far more hydrophobic interface than a natural antibody
control (apolar fraction of buried area 0.65 and 0.76 versus 0.45) and made fewer polar
contacts (9 and 12 hydrogen bonds versus 15; 4 and 4 salt bridges versus 8), at
comparable buried area. Excess interface hydrophobicity is a recognised non-specific
binding and aggregation signature.
**Catch it:** compute interface descriptors against a natural control, not in isolation.
Call `interface_descriptors()`.

### 9. A geometrically perfect site can still be electrostatically hostile
This is the one that survives every other check. After the geometry limitation of
failure mode 8's campaign was removed — all three conserved, free, exposed carboxylates
brought within side-chain reach by re-docking plus a format change — histidines placed
beside them still did not switch: **0 of 21 sites live, 3 inert, 18 inverted**, and the
three-histidine variant was more inverted than its unmutated parent.

The cause is visible in the pKa values: *every* engineered histidine **lost** pKa on
binding (free 5.7-6.7 -> bound 4.0-5.6). Burial desolvates the imidazole faster than the
neighbouring carboxylate stabilises its protonated form, so the bound state is biased
toward neutral histidine and binding is favoured at the *higher* pH. One site started at
free pKa 3.99 and reached -0.36 when bound, contributing exactly nothing.

Critically, this was **reproducible, not noise** — one site read 0.71 / 0.73 / 0.72
across three independent poses, another 0.89 / 0.88 / 0.89, against a measured noise
band only 1.43-fold wide. A tight spread is usually reassuring; here it confirmed a real
and unwanted result.
**Catch it:** before investing in re-docking or a format change to reach more partners,
test the cheap version first — put one histidine at the single partner you *can* already
reach and measure its linkage. If it inverts, more partners will not rescue it, because
the problem is the local dielectric environment rather than the count of sites. The
diagnostic is the sign of dpKa on binding, not the number of carboxylates.

## 6. Interface quality against a natural control

Absolute interface numbers are hard to judge; ratios against a real binder are not. Fold
a natural antibody or native partner for the same target and compare:

- **buried surface area per side** - roughly 600-900 A^2 is normal for a VHH-antigen
  interface; much less suggests a weak interface, much more suggests an implausible one;
- **polar contact count** - hydrogen bonds and salt bridges across the interface;
  substantially fewer than the natural control is a specificity concern;
- **apolar fraction of buried area** - natural interfaces sit near 0.45; designed ones
  drifting to 0.65-0.75 carry aggregation and off-target risk.

## 7. Developability as a hard gate, ahead of the window

Rank on developability **first**, then window. In the campaign the highest-window design
of the entire effort carried an **unpaired cysteine** and **14 methionines** against one
in the parent scaffold - a disulfide-scrambling and oxidation liability that would have
compromised expression and shelf life. Ranking on the objective alone promotes
unmanufacturable molecules.

Minimum flags: unpaired or non-canonical cysteines, N-glycosylation sequons (`N[^P][ST]`),
methionine load versus the parent, long hydrophobic runs, and known target-specific sites
to avoid. Call `developability_flags()`.

Note the trap this creates: repairing a liability changes the sequence, and per failure
mode 3 that **invalidates the measured window** rather than preserving it. Repair, then
re-fold and re-gate from scratch.

## 8. Check the target construct and its numbering

Two cheap checks that are easy to skip and expensive to get wrong.

**Numbering convention.** Deposited structures often use mature-protein numbering while
databases use precursor numbering. In the campaign the reference structure was offset by
**24 residues**; read directly, only 50 of 609 residues matched, and an entire epitope
analysis compared the wrong positions. **Always confirm numbering by sequence identity
before mapping an epitope onto a structure**, testing plausible offsets and taking the
one that matches.

**Construct truncation.** If affinity is specified against a full ectodomain but you are
folding an isolated domain, check whether the epitope is occluded in the intact protein.
Compare epitope SASA in the isolated domain against the intact structure. In the campaign
only 3.2% was occluded, so the fragment was a fair proxy - but that was a result, not an
assumption, and the reverse outcome would have invalidated the whole target choice. Note
which conformation your reference structure represents; a receptor in an extended
ligand-bound state may bury quite different surface from its autoinhibited form.

## 9. Knowing when to stop

A campaign should have a stopping rule, and "we have not hit the target yet" is not one.
Stop computational iteration when:

- your best **reproducible** signal sits within a small factor of the measured noise
  floor (failure mode 4), and
- a targeted intervention derived from your own data has already failed to improve things
  (in the campaign, cycle 3 regressed), and
- **either** the geometry caps live sites below what the specification needs
  (Section 2), **or** the geometry is adequate and the sites still invert reproducibly
  (failure mode 9).

Those last two are different diagnoses with different next steps, and distinguishing
them is worth one cycle. If geometry is the cap, a new binding mode can help — but
verify it by *measuring* a site, not by counting partners. If sites invert with partners
in reach, the route is closed: no amount of re-docking, sequence design or sampling
fixes an interface that destabilises protonated histidine.

At that point further sweeps optimise against noise. The informative next steps are
**experimental measurement** or a **better physics engine** — constant-pH MD gives
pose-averaged protonation instead of single-snapshot predictions, and is the one
calculation that could overturn an inversion verdict, since empirical pKa predictors are
least reliable at exactly the charged, partly-buried interfaces where switches live.

When you hand over candidates, prefer a pair that fails on **complementary** criteria
over a single best-scoring one. Two candidates - one on-epitope but chemically liable,
one clean but marginally off-epitope - discriminate "the predictions are directionally
real" from "the ranking is noise" in a way neither does alone.

## Helpers

`kernel.py` auto-loads with this skill and needs Biopython and NumPy for the structural
functions:

- `ph_window_ceiling(lo, hi)` - max achievable window for 1-4 ideal sites
- `sites_needed(target_window, lo, hi)` - ideal sites a specification implies
- `site_linkage(pka_free, pka_bound, lo, hi)` - one site's contribution
- `classify_site(linkage)` - live / inert / inverted
- `combined_window(sites, lo, hi)` - product over sites, with an inert-site count
- `replicate_verdict(windows)` - geometric mean, range, spread, reproducibility flag
- `noise_floor(free_pkas, bound_pkas, lo, hi)` - implied window band from an invariant
  control residue
- `interface_descriptors(pdb, target_chain, binder_chain)` - BSA, contacts, hydrogen
  bonds, salt bridges, apolar fraction
- `developability_flags(seq, parent_seq)` - liability flags
- `find_carboxylate_clusters(pdb, chain, max_span, min_sasa)` - candidate switch-partner
  clusters on a target surface
