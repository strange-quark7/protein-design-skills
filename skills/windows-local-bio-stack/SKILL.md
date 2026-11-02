---
name: windows-local-bio-stack
description: "Conventions for the local Windows + RTX 4060 protein-design workstation: where persistent model weights, tool repos and caches live so they survive workspace sweeps; how to tell a sandbox network block from a genuine host DNS failure; and the per-tool launch flags (Boltz-2, RFdiffusion, ProteinMPNN, ProtonPottsMPNN, OpenMM) that this machine needs. Load at the START of any structural-biology or binder-design session on this machine, before installing, downloading weights, or concluding that something is unavailable."
---

# Windows local bio stack

Machine-specific operating conventions. Load this **first** in any
structural-biology session on this workstation — most of it is knowledge that
costs an hour to rediscover and thirty seconds to read.

## 0. The two rules that waste the most time

**1. The session workspace is swept. Nothing you download belongs in it.**
Conda environments live outside the workspace and persist. Model weights, git
clones and caches written to the workspace **do not** — an idle sweep reclaimed
4.3 GB of Boltz-2 checkpoints in the middle of a campaign, and they could not be
re-fetched because DNS was down at the time. Persistent assets live under
`C:\Users\drani\cs-assets` (see §1). Point tools at that path *before*
downloading, never after.

**2. A proxy 403 does NOT mean the user is offline.** See §2. Diagnose before
reporting; the usual cause is on the host, not the sandbox, and the user's
browser will look perfectly fine while every system-resolver lookup fails.

## 1. Persistent asset layout

```
C:\Users\drani\cs-assets\
  README.md  REGISTRY.json  SETUP-LOG.md   <- entry point, inventory, history
  env\cs_assets.py         the one import that wires everything up
  env\activate.ps1         same, for a PowerShell session
  bin\doctor.ps1           what is working? assets, envs, GPU, egress
  bin\fetch-weights.ps1    idempotent + resumable weight downloader
  weights\boltz\          BOLTZ_CACHE target: boltz2_conf.ckpt (2.3 GB),
                            boltz2_aff.ckpt, mols\ (CCD, ~45k entries)
  weights\rfdiffusion\     Complex_base_ckpt.pt and friends
  repos\RFdiffusion\       entry point is scripts\run_inference.py
  repos\ProteinMPNN\       clone of dauparas/ProteinMPNN (+ soluble weights)
  repos\ProtonPottsMPNN\   clone of christian-creator/ProtonPottsMPNN
                            (v6 checkpoint ships in-repo, ~21 MB)
  repos\hbplus\            hbplus.exe + standard.data + vdw.radii
  data\targets\  data\scratch\  projects\<name>\  logs\
```

**Do this first, in any project:**

```python
import sys; sys.path.insert(0, r"C:\Users\drani\cs-assets\env")
import cs_assets; cs_assets.activate()
```

It sets `BOLTZ_CACHE`/`HBPLUS_PATH`/`HBDIR` and gives you `python_for(env)`,
`project_dir(name)`, `require("boltz2_conf")` (fail fast instead of 20 min in),
`status()` and `LAUNCH_NOTES`. Stdlib-only, so it imports in every env.
When anything looks broken, run `bin\doctor.ps1` BEFORE theorising — it
separates "missing", "needs a launch flag" and "blocked by NAT64", which all
look identical from a traceback.

Note `ccd.pkl` is NOT needed: Boltz-2 reads the CCD from `mols\`, and
`ccd.pkl` belongs to the Boltz-1 code path only.

Request access once per project with
`request_host_access("~/cs-assets", mode="rw")`; the grant persists. The
directory must already exist — a grant on a non-existent path is refused with
"Directory could not be accessed", which reads like a permissions problem and
is not one.

Point tools at it rather than copying into the workspace:

```python
os.environ["BOLTZ_CACHE"] = r"C:\Users\drani\cs-assets\weights\boltz"
```

## 2. Network: is it the sandbox, or the host?

Egress leaves through a host-side filtering relay; in-sandbox DNS always fails
by design, so `socket.getaddrinfo` proving nothing is expected. Read the error:

| symptom | meaning | action |
|---|---|---|
| `403 ... host X **not on the allowlist**` | genuinely off-allowlist | `request_network_access(domain=...)` — usually granted |
| `403 ... host X **in a private/reserved range**` for a host you know is public | **DNS64/NAT64** — see below. Not a block, and not grantable | host-side fix |
| `401`, or a plain `403` with no sandbox/policy wording | you reached the server; it refused you | ordinary auth error |

Read the exact wording: "not on the allowlist" and "in a private/reserved
range" are *different checks* with different fixes, and conflating them wastes
the most time of anything in this file.

**The private/reserved case is DNS64, not poisoning** (diagnosed 2026-10-11;
the earlier "a DNS filter is returning private IPs" reading was WRONG). The
host's resolver performs **DNS64 synthesis**: alongside correct public A
records it returns AAAA records in `64:ff9b::/96`, the RFC 6052 well-known
NAT64 prefix. The egress proxy prefers the AAAA, sees a special-purpose prefix
and rejects it as reserved — although the embedded IPv4 is a perfectly good
public address. Nothing is filtered; the DNS is *correct*.

Expect it to be **partial, and to flip mid-session**: hosts without synthesized
AAAAs (e.g. `conda.anaconda.org`) keep working while `huggingface.co`,
`pypi.org` and `files.ipd.uw.edu` fail, and a cached good entry can expire into
a failing one. "It worked ten minutes ago" does not refute this diagnosis.

The user's browser keeps working regardless, because Chrome and Firefox default
to DNS-over-HTTPS with their own resolver — so "my internet is fine" and "the
agent cannot reach anything" are both true at once. One command settles it, run
by the user in a host terminal:

```powershell
Resolve-DnsName huggingface.co | Select-Object Name, IPAddress
```

Decode what comes back rather than pattern-matching it:

* `64:ff9b::de1:564` -> low 32 bits `0d.e1.05.64` = `13.225.5.100`, which will
  match one of the A records. That is **DNS64**.
* a `10.x` / `172.16-31.x` / `192.168.x` / `0.0.0.0` answer with **no** matching
  public A record would be genuine filtering (Pi-hole / AdGuard / NextDNS).

Fix for DNS64, any one of: disconnect the VPN (**Cloudflare WARP uses DNS64 by
default** and is the usual culprit), pin the adapter's DNS to `1.1.1.1`, or
disable IPv6 on the adapter. Then restart the backend and run
`bin\fetch-weights.ps1`. **Never** work around it with a mirror, cache or
alternate endpoint.

## 3. Interpreter and path traps

- Bare `python` in a shell resolves to a **Python 2.7 from an MGLTools install**.
  Always call an explicit env interpreter:
  `C:\Users\drani\.claude-science\conda\envs\<env>\python.exe`.
- `$HOME` in the powershell tool expands to a sandbox-local `.home`, not
  `C:\Users\drani`. Write user paths out literally.
- `wsl.exe` cannot be invoked from the sandbox at all (`Wsl/E_ACCESSDENIED`).
  To use a WSL-only tool, write a self-diagnosing bash script into the
  workspace and have the user run it in their own WSL terminal.
- Managed conda env site-packages are mounted **read-only** in kernels: a
  `pip install` into one appears to succeed and writes nothing. Use
  `manage_packages`.

## 4. Per-tool launch requirements

**Boltz-2** — needs a `fairscale` shim on `PYTHONPATH` (the real package cannot
be installed here; boltz imports exactly one symbol, `checkpoint_wrapper`). Run
with `--no_kernels` (triangle kernels are Linux-only), `--num_workers 0`, and
**`--preprocessing-threads 1`** — the sandbox denies `CreateNamedPipe`, so
`multiprocessing.Pool` dies with `WinError 5`. That only triggers with more than
one input, so a batch exits 0 with **no output** in ~27 s and looks like a
silent success. Budget ~5 min per design for a ~540-residue complex.

**RFdiffusion** — runs natively on Windows (it does *not* need WSL); the
`rfdiff` env already has torch 2.3, DGL, e3nn, hydra. Only the weights are
missing. ~8.5 min/design on a 536-residue trimer.

**ProteinMPNN** — must be pinned to CPU with `CUDA_VISIBLE_DEVICES=-1`. The
literal `-1`, not an empty string: an empty value is not honoured on Windows and
the process will contend with whatever is on the GPU, producing
`CUDA error: unspecified launch failure` in the *other* job. Three Windows bugs:
weights are located with `rfind("/")` (returns -1 on backslash paths) so pass
`--path_to_model_weights`; `parse_multiple_chains.py` writes full paths into
the record `name`, breaking bias lookup and output filenames — rewrite names to
bare stems; and the same path bug hits the **single `--pdb_path` route**, where
an absolute path builds `out//seqs/C:\...\x.fa` and dies with `OSError:
[Errno 22] Invalid argument` — `cd` to the PDB's directory and pass a **bare
filename**. Also: `%TEMP%` is unusable as a working dir (the sandbox denies
`Push-Location`/`Remove-Item` there) — use the workspace or `data\scratch\`.

**ProtonPottsMPNN** — needs an `assertpy` shim (sdist-only), HBPLUS at the
**nested** `hbplus/hbplus/hbplus.exe`, both `HBPLUS_PATH` *and* `HBDIR` set, and
`PATH` sanitised to just the env plus `System32` (otherwise a delay-load DLL
error `0xC06D007F`). See the `protonpottsmpnn` skill for the polarity-inversion
API.

**OpenMM** — the **CUDA platform is unusable**: the driver exposes CUDA 13.1
while the env ships NVRTC 13.4, so OpenMM emits PTX the driver rejects
(`CUDA_ERROR_UNSUPPORTED_PTX_VERSION`). Use **OpenCL**, which runs on the same
RTX 4060 — and pass no `DeviceIndex` unless you also pass
`OpenCLPlatformIndex`. Tune it: implicit-solvent GBn2 on a 8292-atom system goes
37 ns/day naive -> 55 with a 1.8 nm `CutoffNonPeriodic` -> **178 ns/day with
hydrogen-mass repartitioning (`hydrogenMass=4*amu`, `constraints=HBonds`) at a
4 fs timestep**. A 1.2 nm cutoff is *slower* than 1.8 nm. Boltz-predicted PDBs
lack C-terminal `OXT`, so run PDBFixer `findMissingAtoms`/`addMissingAtoms`
first (set `missingResidues={}` so it does not build loops) or OpenMM fails with
"No template found for residue N (LEU) ... missing 1 C atom".

## 5. Environments that already exist

`binder-design` (py3.11) is the analysis workhorse: biopython, freesasa, propka,
**prodigy_prot** (protein-protein affinity, has a `prodigy` CLI), openmm 8.6.1,
pdbfixer, scipy, matplotlib — and it runs ProteinMPNN. `folding` has boltz
2.2.1. `rfdiff` has torch 2.3 + DGL. `protonpotts` is py3.12.

Check with `manage_environments(mode="list", dependencies=[...])` before
creating anything new.
