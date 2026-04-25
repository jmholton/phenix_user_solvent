# phenix_user_solvent

Feature branch adding user-supplied bulk solvent map support to phenix.refine.
Changes are applied to `/programs/phenix-2.1rc2-6037/` (writable).

**Important:** The default `phenix.refine` in PATH is version 2.0-5936 (unmodified).
Always use `/programs/phenix-2.1rc2-6037/bin/phenix.refine` to run the patched version.

## What was changed and why

Four changes in the phenix-2.1rc2-6037 installation:

**`lib/python3.9/site-packages/phenix/refinement/__init__.params`**
Added `refinement.input.bulk_solvent_map` scope inside existing `input { }` scope:
`file_name`, `amplitudes_label` (default `FPART`), `phases_label` (default `PHIPART`).

**`lib/python3.9/site-packages/phenix/programs/phenix_refine.py`**
Added `_inject_user_bulk_solvent(fmodel, params, log)` at module level.
Called in `Program.validate()` after `self.fmodel` is created (single
refinement path only; joint XN path not yet wired up).

The function reads the MTZ, finds the two columns by label, calls
`phase_transfer(deg=True)` to get a complex Miller array, maps to ASU,
then uses `miller.match_indices` to build a new complex array with exactly
`f_obs`'s index set and ordering. Missing reflections get zero amplitude.
Calls `fmodel.set_user_f_masks([f_mask])`.

**`lib/python3.9/site-packages/mmtbx/f_model/f_model.py`** (three changes):
- Added `self._user_f_masks = None` in `manager.__init__`
- Added `manager.set_user_f_masks(f_masks)` method just before `update_xray_structure`
- In `update_xray_structure`, when `update_f_mask=True` and `_user_f_masks` is set,
  uses the stored masks instead of calling `mask_manager.shell_f_masks()`
- In `manager.select()`, propagate `_user_f_masks` (selecting the same reflections)
  so that outlier removal (`remove_outliers` → `select(in_place=True)`) does not
  silently reset `_user_f_masks` to `None`

## Key design decisions

- The user map **replaces** the flat 0/1 mask entirely; there is no hybrid fallback
- `k_mask` continues to be determined analytically per resolution bin (Afonine et al. 2013,
  https://doi.org/10.1107/S0907444913000462) — the old exponential `k_sol*exp(-B_sol*s²/4)`
  is not what phenix actually uses
- The mask is frozen for the entire run (not recomputed as atoms move)
- Index alignment uses `match_indices` not `common_set` — the latter gives
  wrong ordering. The root cause of the initial AssertionError was an
  ASU ordering mismatch, not a resolution gap.
- `_user_f_masks` must be propagated through `manager.select()` (including the
  `in_place=True` path used by `remove_outliers`). Without this fix, outlier
  removal silently resets `_user_f_masks` to `None` and the flat mask is
  recomputed during the LBFGS coordinate minimisation, wrecking refinement.

## Bug hunt: why did final R-factors converge?

Symptom: the user mask gave much better post-BSS R-work (0.1087 vs 0.1768)
but final R-factors were nearly identical between user and default runs.

Root cause: `f_model_all_scales.run.compute()` calls `remove_outliers()` five
times before fitting k_mask. `remove_outliers()` calls
`self.select(selection, in_place=True)`. `select(in_place=True)` creates a new
`manager` via `__init__` (setting `_user_f_masks = None`), then copies ALL its
attributes back to `self` via dict iteration — wiping `_user_f_masks`. The
subsequent `update_xray_structure(update_f_mask=True)` then recomputed the
flat mask from the atomic model. k_mask was then fit to the flat mask,
not the user mask, so coordinate refinement diverged and reached the same final
R as the default run.

## Testing

**Quick test (1aho, small, ~30 s):**
```
/programs/phenix-2.1rc2-6037/bin/phenix.refine \
  example/1aho.pdb example/1aho.mtz \
  refinement.input.bulk_solvent_map.file_name=example/solvent_Fpart.mtz \
  refinement.input.bulk_solvent_map.amplitudes_label=Fpart \
  refinement.input.bulk_solvent_map.phases_label=PHIpart \
  refinement.main.number_of_macro_cycles=1 \
  "refinement.refine.strategy=individual_sites individual_adp occupancies" \
  output.prefix=test_1aho --overwrite
```
Expected: R-work ~0.151, R-free ~0.147.

**Large ensemble test (evenmoreconf.pdb, ~5 min):**
```
/programs/phenix-2.1rc2-6037/bin/phenix.refine \
  evenmoreconf.pdb refme_minRfree.mtz \
  refinement.input.bulk_solvent_map.file_name=solvent_Fpart.mtz \
  refinement.input.bulk_solvent_map.amplitudes_label=Fpart \
  refinement.input.bulk_solvent_map.phases_label=PHIpart \
  refinement.main.number_of_macro_cycles=1 \
  "refinement.refine.strategy=individual_sites individual_adp occupancies" \
  output.prefix=test_emc --overwrite
```
Expected: Final R-work ~0.108, R-free ~0.119 (vs default ~0.177).

**Very large ensemble test (48 copies, 45762 atoms, ~25 GB RAM, slow):**
```
/programs/phenix-2.1rc2-6037/bin/phenix.refine \
  refmacout_minRfree.pdb refme_minRfree.mtz \
  refinement.input.bulk_solvent_map.file_name=refme_minRfree.mtz \
  refinement.input.bulk_solvent_map.amplitudes_label=Fpart \
  refinement.input.bulk_solvent_map.phases_label=PHIpart \
  refinement.main.number_of_macro_cycles=1 \
  "refinement.refine.strategy=individual_sites individual_adp occupancies" \
  output.prefix=test_usermap --overwrite
```

## Known limitations

- Joint X-ray/neutron refinement path in `validate()` is not wired up
- Twin refinement untested (user f_mask_twin is not set)
