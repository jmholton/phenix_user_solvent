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
analytically per resolution bin.

This patch allows the user to supply their own bulk solvent structure
factors (amplitude + phase columns from an MTZ file) in place of `F_mask`.
`k_mask` continues to be determined analytically against the user-supplied
map. The mask is not recomputed from the atomic model during refinement.

A key use case is supplying a bulk solvent model computed from an ensemble
or explicit-solvent MD simulation, where the solvent structure is better
defined than the flat-mask approximation.

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

See `example/README.txt`. The solvent model in `example/solvent_Fpart.mtz`
was computed by Refmac using a 48-copy ensemble of the same crystal.
Results after 1 macro cycle (no RSR):

| Structure | Bulk solvent | R-work | R-free |
|---|---|---|---|
| 1AHO (single chain) | Input | 0.152 | 0.146 |
| 1AHO (single chain) | Default (flat mask) | 0.151 | 0.147 |
| 1AHO (single chain) | Ensemble Fpart/PHIpart | 0.151 | 0.146 |
| Multi-conf | Default (flat mask) | 0.177 | 0.180 |
| Multi-conf | Ensemble Fpart/PHIpart | **0.108** | **0.119** |

For a simple single-chain structure the two models give essentially
identical results. The benefit is larger for the multi-conformer crystal
form, where the flat-mask approximation is less accurate.

## Applying the patch

From the `site-packages` directory of your phenix installation:

```bash
cd /path/to/phenix/lib/python3.9/site-packages
patch -p0 < f_model.patch
patch -p0 < phenix_refine.patch
patch -p0 < __init__.params.patch
```

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

## Reference

Afonine, P. V., Grosse-Kunstleve, R. W., Adams, P. D. & Urzhumtsev, A. (2013).
*Acta Cryst.* D**69**, 625–634.
https://doi.org/10.1107/S0907444913000462
