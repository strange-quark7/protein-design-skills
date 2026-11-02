---
name: concave-epitope-geometry
description: "Map the accessible volume of a concave epitope - an MHC peptide groove, enzyme active-site cleft, receptor canyon or antibody-inaccessible pocket - and measure how deeply a designed binder inserts into it, calibrated against a natural binder of the same target. Reach for this skill when designing against a target whose epitope is a cleft rather than a convex patch, when deciding whether a protruding loop or finger could reach a buried residue, when a reviewer asks which surface a design engages and how deeply, or when choosing hotspot residues for a pocket. It answers three questions the usual binder metrics cannot: is there room in there at all, where along the cleft is the room, and is this design actually using it."
---

# Concave epitope geometry

Standard binder metrics - buried area, contact counts, ipTM - treat an
interface as a surface. For a cleft they hide the question that matters: **is
there space inside, where along the cleft is it, and does the design reach
it?** This skill supplies a reproducible coordinate frame, a volume-vs-depth
profile, and an insertion-depth metric calibrated on a reference binder.

Load `kernel.py` helpers by loading this skill; they work on any Bio.PDB model.

---

## 1. Build the frame from the structure, never hard-code it

A cleft has a long axis, a transverse axis and a normal. Derive all three from
the structure so the numbers survive a different model, allele or crystal:

* **long axis `u`** - first principal component of the ligand/substrate Calpha
  (or, with no ligand, of the cleft-lining residues' centroids)
* **transverse `v`** - vector between the two flanking-wall centroids,
  orthogonalised against `u`
* **normal `n`** - `u x v`, signed to point AWAY from the floor (check the
  floor centroid's projection and flip if needed)

`groove_frame(model, floor_resids, wall_a_resids, wall_b_resids, ligand_chain)`
returns `(origin, u, v, n)`. Heights are measured along `n` from the origin.

## 2. The rim is the LOWER wall

Access is limited by the lower of the two walls, not the mean and not the
higher one. Using the higher wall flatters every design. Take the maximum
height over wall atoms within a window along `u` (the walls rise away from the
cleft, so an unwindowed max picks up a distant shoulder).

## 3. Profile volume against depth, binned along the long axis

The single most useful output is **open volume as a function of depth below
the rim, binned along `u`**. Grid the region at 0.5 A with a 1.4 A water
probe; a point is open if it clears every heavy atom by vdW + probe.

**Do not report maximum depth per bin.** One narrow crevice dominates it and
invents a pocket that nothing can enter. Report the *volume* beyond a depth
threshold - that is what a side chain or loop needs.

On an HLA-A\*11:01 peptide groove this profile showed ~620 A3 open below the
rim overall, but only ~84 A3 more than 4 A down, and that residual volume sat
over peptide positions 8-10 (64 A3) rather than over the position-6 mutation
the campaign was targeting (2 A3). The design question changed completely.

Report the *localisation*, not just the totals. Absolute volumes move by about
10% depending on whether the lateral window's boundary planes are counted;
the C-terminal-to-central ratio does not, and that is what a design decision
should rest on.

## 4. Calibrate on a natural binder before judging a design

Insertion depth means nothing in the abstract. Measure it for a known binder
of the same target - a TCR, a natural substrate, a published antibody - and
quote every design against that. Depth alone is not evidence of a good
interface; a reference value makes it interpretable.

`insertion(model, binder_chains, ...)` returns `max_depth` (rim height minus
the binder's lowest atom, counting only atoms laterally between the walls),
`n_below_rim`, and `u_of_deepest` - the position along the cleft where the
binder reaches deepest. That last field is often the informative one: it tells
you *which end* of the cleft the design is exploiting.

## 5. Check whether your designs already insert

Before building a hotspot set to force insertion, measure the designs you
already have. In the campaign that produced this skill, 14 of 20 binders
generated with ordinary surface-directed hotspots already inserted deeper than
the natural receptor, and their deepest atoms sat in the same sub-pocket the
receptor uses. Hotspot conditioning on an epitope that *contains* a cleft
tends to find the cleft on its own.

## 6. Measure it; do not select on it

Insertion depth is a descriptor, not a validated objective. Measured across
two design rounds it predicted:

* **nothing** about target-state dependence (Spearman rho -0.08, p = 0.74)
* **nothing** about pose reproducibility across resampling (rho -0.07, p = 0.60)
* **nothing** about tag accessibility at the opposite terminus (rho +0.004)

A hotspot set aimed explicitly at the deep sub-pocket produced no significant
advantage over a length-matched control on any measured quantity, including
the burial it was designed to increase (p = 0.89). Report depth as a per-design
column and a tie-break. Giving it weight in a ranking is not supported.

## 7. Pitfalls

* **Hard-coded helix ranges.** The wall residue ranges are target-specific
  arguments, not constants. Pass them explicitly; a different allele or a
  different fold needs different numbers.
* **Backbone-only models.** Generated backbones carry no side chains, so depth
  measured on them is understated by several angstroms. Compare like with
  like - backbone-only against a backbone-only reference, full-atom against
  full-atom.
* **Lateral window.** Count only binder atoms between the two walls. Without
  it, a binder lying alongside the cleft scores as if it were inside.
* **Centroid is not contact.** A binder much larger than the cleft has its
  centroid well outside it - the natural TCR in the calibrating campaign sits
  at u = +14.8, past the ligand's C-terminus. Judge engagement by hotspot
  contacts, not by where the centre of mass lands.

## 8. Reporting

For each design: `max_depth`, `n_below_rim`, `u_of_deepest`, each alongside the
reference binder's value. For the target once: the volume-depth profile table
and the residues lining any deep sub-pocket (`pocket_lining`) - the latter
often names the mechanistically interesting residues, including ones that are
invisible from the outer surface.
