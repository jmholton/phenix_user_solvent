# phenix_user_solvent

Adds support for a user-supplied bulk solvent electron density map to
`phenix.refine`, equivalent to the `FPART`/`PHIPART` feature in Refmac.

## Motivation

`phenix.refine` models the total structure factor as (Afonine et al., 2013):

```
F_model = k_total * (F_calc + k_mask * F_mask)
```

where `F_mask` is computed from a flat (0/1) solvent mask derived from the
atomic model, and `k_mask` is a resolution-dependent scale factor determined
analytically per resolution bin ([Afonine et al., 2013](https://doi.org/10.1107/S0907444913000462)).

This patch allows the user to supply their own bulk solvent structure
factors (amplitude + phase columns from an MTZ file) in place of `F_mask`.
`k_mask` continues to be determined analytically against the user-supplied
map. The mask is not recomputed from the atomic model during refinement.

A key use case is supplying a bulk solvent model computed from an explicit-solvent
MD simulation (e.g. AMBER) or other sources, where the solvent
structure may be better defined than the flat-mask approximation.

## Files changed

| File | Repo | Lines | Change |
|---|---|---|---|
| `phenix/refinement/__init__.params` | phenix | +23 | New `refinement.input.bulk_solvent_map` PHIL scope |
| `phenix/programs/phenix_refine.py` | phenix | +46 | `_inject_user_bulk_solvent()` helper; called in `validate()` |
| `mmtbx/f_model/f_model.py` | cctbx_project | +25 −12 | `set_user_f_masks()` method; guard in `update_xray_structure()`; propagation through `select()` |

## Usage

```
phenix.refine model.pdb data.mtz \
  refinement.input.bulk_solvent_map.file_name=solvent.mtz \
  refinement.input.bulk_solvent_map.amplitudes_label=Fpart \
  refinement.input.bulk_solvent_map.phases_label=PHIpart
```

Default column labels are `FPART` and `PHIPART`. The solvent MTZ can be
the same file as the data MTZ if the columns are present there.

## Example

The solvent model in `example/solvent_Fpart.mtz` was computed from an AMBER
MD trajectory of the same crystal, followed by density editing in which strong
difference features from Fo-Fc maps (computed without free-flagged reflections)
were added back. `example/starthere.pdb` is the starting multi-conformer model;
`example/refme.mtz` contains the experimental data with the original R-free flags.

To reproduce the user-supplied solvent result (log: `example/usersolvent_001.log`):

```
phenix.refine \
  example/starthere.pdb \
  example/refme.mtz \
  refinement.input.bulk_solvent_map.file_name=example/solvent_Fpart.mtz \
  refinement.input.bulk_solvent_map.amplitudes_label=Fpart \
  refinement.input.bulk_solvent_map.phases_label=PHIpart \
  refinement.main.number_of_macro_cycles=3 \
  refinement.refine.strategy="individual_sites+individual_adp" \
  output.prefix=usersolvent
```

To reproduce the default flat-mask result (log: `example/defaultsolvent_001.log`):

```
phenix.refine \
  example/starthere.pdb \
  example/refme.mtz \
  refinement.main.number_of_macro_cycles=3 \
  refinement.refine.strategy="individual_sites+individual_adp" \
  output.prefix=defaultsolvent
```

| Bulk solvent model | R-work | R-free |
|---|---|---|
| Default flat mask | 0.128 | 0.147 |
| Ensemble Fpart/PHIpart, cycle 1 | 0.093 | 0.109 |
| Ensemble Fpart/PHIpart + density editing, iterated | **0.084** | **0.105** |

The improvement compounds across iterations (see [Iterative workflow](#iterative-workflow) below).

## Applying the patch

First, clone this repository somewhere convenient:

```bash
git clone https://github.com/jmholton/phenix_user_solvent.git
```

Then apply the patches. The target directory depends on your phenix release:

**phenix 2.x (2.0-5761, 2.1rc2-6037, …)**

The files to patch live under `lib/python3.*/site-packages/`:

| Patch | File patched |
|---|---|
| `f_model.patch` | `lib/python3.*/site-packages/mmtbx/f_model/f_model.py` |
| `phenix_refine.patch` | `lib/python3.*/site-packages/phenix/programs/phenix_refine.py` |
| `__init__.params.patch` | `lib/python3.*/site-packages/phenix/refinement/__init__.params` |

```bash
cd /programs/phenix-2.x-NNNN/lib/python3.*/site-packages
patch -p0 < ~/phenix_user_solvent/f_model.patch
patch -p0 < ~/phenix_user_solvent/phenix_refine.patch
patch -p0 < ~/phenix_user_solvent/__init__.params.patch
```

**phenix 1.21 and earlier**

These releases use a source-tree layout; the files live under `modules/`:

| Patch | File patched |
|---|---|
| `f_model.patch` | `modules/cctbx_project/mmtbx/f_model/f_model.py` |
| `phenix_refine.patch` | `modules/phenix/phenix/programs/phenix_refine.py` |
| `__init__.params.patch` | `modules/phenix/phenix/refinement/__init__.params` |

```bash
cd /programs/phenix-1.21-NNNN/modules
patch -p0 < ~/phenix_user_solvent/f_model.patch
patch -p0 < ~/phenix_user_solvent/phenix_refine.patch
patch -p0 < ~/phenix_user_solvent/__init__.params.patch
```

> **Note:** All releases contain other files named `phenix_refine.py`
> (e.g. `phenix/command_line/phenix_refine.py`) that must **not** be patched.
> The patch targets `phenix/programs/phenix_refine.py` only.

## Iterative workflow

The per-bin `k_mask` in phenix exposes a limitation of preprocessing pipelines
(e.g. Refmac) that fit a single global B-factor to the MD solvent map. That
global B-factor blurs the map uniformly, suppressing high-resolution solvent
content while over-weighting the dominant low-frequency terms. Genuine
ordered-water peaks are suppressed, appear as Fo-Fc features, and get patched
with explicit coordinate waters — a compensatory loop that inflates model
complexity without improving the underlying solvent description.

**phenix breaks this loop.** With `k_mask` fit independently per resolution
bin, it can correct the resolution-dependent miscalibration:

| Resolution (Å) | k_mask (flat mask) | k_mask (MD Fpart) |
|---|---|---|
| 14–6.7 | 0.33 | 1.05 |
| 5.4–4.4 | 0.32 | 1.12 |
| 2.9–2.4 | 0.00 | 1.12 |
| 2.4–1.9 | 0.00 | 1.32 |
| 1.5–1.3 | 0.00 | 1.45 |

The flat mask contributes nothing above 3 Å. The MD map contributes across
the full resolution range, with k_mask > 1 at high resolution — meaning the
global B-factor had undersold the map amplitude there. Once phenix corrects
each bin independently, the explicit coordinate waters that were patching the
blurred-map deficit become less necessary, the model simplifies, and R-free
improves further.

**After phenix refinement, any remaining Fo-Fc features in the solvent are
genuine missing structure**, not blurring artifacts. These become input for the
next density-editing cycle. The recommended iteration is:

1. Refine with user solvent map (`individual_sites occupancies`, 3–10 cycles)
2. Inspect Fo-Fc map (`output.export_final_f_model=True` to get `FMASK` and `K_MASK`)
3. Density-edit: add positive solvent blobs back into the solvent MTZ
4. Repeat from step 1 until solvent-region Fo-Fc features fall below ~3σ

Convergence indicators: k_mask profile stabilizes; number of >5σ solvent peaks
decreases monotonically; R-free decreases monotonically.

## Implementation notes

The user-supplied MTZ columns are read once in `validate()`, mapped to ASU,
and aligned to `f_obs`'s exact index set using `miller.match_indices`.
Reflections present in `f_obs` but absent from the user map receive zero
contribution (in practice this is rare: the map is typically computed to
at least the resolution of the data). The aligned complex Miller array is
stored as `fmodel._user_f_masks` and used in place of the mask-manager
output for the lifetime of the run.

`_user_f_masks` is propagated through `manager.select()` so that outlier
removal does not silently reset it to `None` at the start of each scaling
cycle.

## Known limitations

- Joint X-ray/neutron refinement path in `validate()` is not wired up
- Twin refinement untested (user f_mask_twin is not set)

## Reference

Afonine, P. V., Grosse-Kunstleve, R. W., Adams, P. D. & Urzhumtsev, A. (2013).
*Acta Cryst.* D**69**, 625–634.
https://doi.org/10.1107/S0907444913000462
