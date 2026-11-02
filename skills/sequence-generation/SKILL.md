---
name: sequence-generation
description: Generate large candidate libraries of ssDNA aptamer sequences against a protein target, then characterize and rank them. Use when the user wants to generate N candidate aptamer/oligo sequences (e.g. "generate 10000 aptamer candidates for ferritin"), build a de novo SELEX starting library with a randomized core flanked by primer regions, or fold/score/rank existing candidate sequences by predicted structure stability (dG), Tm, and GC content. Ships kernel.py helpers: generate_candidate_library, characterize_library, rank_library.
---

# Sequence generation — candidate aptamer libraries

Generate, characterize, and rank de novo ssDNA aptamer candidate libraries for a
protein target. This is a computational **starting library** for an in-silico
screen / SELEX round — scores are structure-informed prioritizations, **not**
measured binding affinities. Candidates still require experimental validation.

## When to use

- "Generate N candidate aptamer sequences for <target>" (ferritin, a venom toxin, etc.)
- Build a SELEX starting library: randomized core flanked by fixed primer regions
- Fold / score / rank existing candidate sequences by predicted stability, Tm, GC

## Library architecture (standard SELEX layout)

Each candidate is `5'-[20-nt forward primer]-[N-core randomized]-[20-nt reverse-primer site]-3'`.
Default core length 40 → 80-nt oligo. Primer regions are **placeholders** — swap
in the user's own primer set via `fwd=`/`rev=`. The randomized core is the
target-binding region and carries the library diversity.

There is rarely a validated public aptamer that binds a given target directly, so
this generates a principled random library rather than seeding from a known binder.
If the user *has* a known binder or motif, mention that seeding/doping the core
around it is an alternative not covered by the default generator.

## Environment

Needs `numpy`, `pandas`, `biopython`, and `seqfold` (pip). matplotlib for the
summary figure. Create/pick an env with these before running.

## Workflow

The `kernel.py` sidecar auto-loads three helpers when this skill is loaded.

1. **Generate** the filtered, unique library (fast, no folding):

   ```python
   df = generate_candidate_library(10000, core_len=40)   # fwd=..., rev=... to swap primers
   ```
   Filters each core: GC 40-60%, no homopolymer run >4, no single base >50%,
   all four bases present, full-library uniqueness. Tune via `gc_range`,
   `max_homopolymer`, `max_base_frac`, `require_all_bases`, `seed`.

2. **Characterize + rank** (parallel folding in a subprocess):

   ```python
   ranked = characterize_library(df)   # adds dG_kcal_mol, Tm_C, GC_*, Gquad_motif, score, rank, aptamer_id
   ranked.to_csv("candidates.csv", index=False)
   ```
   Folding is ~0.13 s/seq; parallelized across cores it is ~1-2 min for 10k.
   `characterize_library` calls `rank_library` internally. To re-rank with
   different weights/target, call `rank_library(df, ...)` directly.

3. **Summarize** with a publication-grade figure (load `figure-style` first):
   ΔG distribution + GC-core-vs-ΔG scatter colored by score.

4. **Save** the CSV, figure, and a short README as artifacts.

## Ranking score

`score = 0.70*(normalized fold stability, more-negative dG better) +
0.25*(closeness of core GC to 50%) + 0.05 (+0.10 if a G-quadruplex motif is present)`.
Sorted descending; rank 1 = most promising fold. Note: the composition filters
suppress G-rich runs, so random N40 cores rarely carry a G-quadruplex motif —
that bonus mainly matters for seeded/relaxed libraries.

## Notes & caveats

- Always state to the user that scores are computational prioritization, not
  measured affinity, and that experimental SELEX/binding validation is required.
- For a diverse experimental panel, cluster top candidates by sequence/structure
  before picking — a top-N slice by score is often near-identical folds.
- seqfold gives ssDNA MFE folding; for pseudoknot-aware or RNA folding, a
  different tool is needed.
