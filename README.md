# phenix_user_solvent

Adds support for a user-supplied bulk solvent electron density map to
`phenix.refine`, equivalent to the `FPART`/`PHIPART` feature in Refmac.

## Motivation

`phenix.refine` currently models bulk solvent by computing a flat (0/1)
mask from the atomic model and scaling it:

```
F_bulk = k_sol * exp(-B_sol * s^2/4) * F_mask
```

This patch allows the user to supply their own bulk solvent structure
factors (amplitude + phase columns from an MTZ file) in place of `F_mask`.
`k_sol` and `B_sol` continue to refine against the user-supplied map.
The mask is not recomputed from the atomic model during refinement.

A key use case is supplying a bulk solvent model computed from an ensemble
or explicit-solvent MD simulation, where the solvent structure is better
defined than the flat-mask approximation.

## Files changed

| File | Repo | Change |
|---|---|---|
| `phenix/refinement/__init__.params` | phenix | New `refinement.input.bulk_solvent_map` PHIL scope |
| `phenix/programs/phenix_refine.py` | phenix | `_inject_user_bulk_solvent()` helper; called in `validate()` |
| `mmtbx/f_model/f_model.py` | cctbx_project | `set_user_f_masks()` method; guard in `update_xray_structure()` |

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
Refining the single-chain 1AHO model against the ensemble solvent:

```
phenix.refine \
  example/1aho.pdb \
  example/1aho.mtz \
  refinement.input.bulk_solvent_map.file_name=example/solvent_Fpart.mtz \
  refinement.input.bulk_solvent_map.amplitudes_label=Fpart \
  refinement.input.bulk_solvent_map.phases_label=PHIpart
```

| | R-work | R-free |
|---|---|---|
| Input | 0.165 | 0.160 |
| After 1 cycle | 0.151 | 0.147 |

## Implementation notes

The user-supplied MTZ columns are read once in `validate()`, mapped to ASU,
and aligned to `f_obs`'s exact index set using `miller.match_indices`.
Reflections present in `f_obs` but absent from the user map receive zero
contribution (in practice this is rare: the map is typically computed to
at least the resolution of the data). The aligned complex Miller array is
stored as `fmodel._user_f_masks` and used in place of the mask-manager
output for the lifetime of the run.
