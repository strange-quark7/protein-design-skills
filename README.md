# Protein-design skills

Personal [Claude Science](https://claude.ai) skills for de novo binder design and
structural biology, authored while running design campaigns on a local
Windows + RTX 4060 workstation.

| skill | what it is for |
|---|---|
| `targeted-binder-campaign` | end-to-end de novo binder campaign against a defined epitope: hotspot-conditioned backbone generation, sequence-design settings that avoid composition collapse, and the co-fold gate stack (interface confidence vs a control, epitope fidelity, pose resampling) |
| `concave-epitope-geometry` | map the accessible volume of a cleft epitope (MHC groove, enzyme active site, receptor canyon) and measure how deeply a binder inserts, calibrated against a natural binder |
| `ph-switch-binder-validation` | design and validate pH-selective binders; the failure modes that invalidate naive designs and the measurement that catches each |
| `protonpottsmpnn` | protonation-aware sequence design with an explicit protonation alphabet |
| `sequence-generation` | sequence generation utilities |
| `windows-local-bio-stack` | operating conventions for a local Windows GPU workstation: where weights and repos live, per-tool launch flags, and how to tell a sandbox network block from a host DNS failure |

Each skill is a directory containing `SKILL.md` (the instructions) and, where
useful, a `kernel.py` sidecar of helper functions.

## Provenance

`concave-epitope-geometry` was derived from a KRAS G12D / HLA-A\*11:01 binder
campaign. Its central results are negative and are stated as such in the skill:
insertion depth predicts neither target-state dependence nor pose
reproducibility, and a hotspot set aimed explicitly at the deep sub-pocket gave
no significant advantage over a length-matched control. It is a measurement
worth making, not an objective worth optimising.
