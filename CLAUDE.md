# phenix_user_solvent

Feature branch adding user-supplied bulk solvent map support to phenix.refine.
Changes are applied to `/programs/phenix-2.1rc2-6037/` (writable).

## What was changed and why

Three files in the phenix-2.1rc2-6037 installation were modified:

**`lib/python3.9/site-packages/phenix/refinement/__init__.params`**
Added `refinement.input.bulk_solvent_map` scope with three parameters:
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

**`lib/python3.9/site-packages/mmtbx/f_model/f_model.py`**
- Added `self._user_f_masks = None` in `manager.__init__`
- Added `manager.set_user_f_masks(f_masks)` method just before `update_xray_structure`
- In `update_xray_structure`, when `update_f_mask=True` and `_user_f_masks` is set,
  uses the stored masks instead of calling `mask_manager.shell_f_masks()`

## Key design decisions

- The user map **replaces** the flat 0/1 mask entirely; there is no hybrid fallback
- `k_sol` and `B_sol` still refine (they scale the user-supplied F_mask)
- The mask is frozen for the entire run (not recomputed as atoms move)
- Index alignment uses `match_indices` not `common_set` — the latter gives
  wrong ordering. The root cause of the initial AssertionError was an
  ASU ordering mismatch, not a resolution gap.

## Testing

Test files in cwd: `1aho.pdb`, `1aho.mtz`, `solvent_Fpart.mtz`
The solvent MTZ has columns `Fpart` (F type) and `PHIpart` (P type).

Quick test command (no RSR, 1 cycle):
```
phenix.refine 1aho.pdb 1aho.mtz \
  refinement.input.bulk_solvent_map.file_name=solvent_Fpart.mtz \
  refinement.input.bulk_solvent_map.amplitudes_label=Fpart \
  refinement.input.bulk_solvent_map.phases_label=PHIpart \
  refinement.main.number_of_macro_cycles=1 \
  "refinement.refine.strategy=individual_sites individual_adp occupancies" \
  output.prefix=test_1aho --overwrite
```

Expected: R-work ~0.151, R-free ~0.147 after 1 cycle.

The big ensemble test (48 copies, 45762 atoms, ~25 GB RAM, slow):
```
phenix.refine refmacout_minRfree.pdb refme_minRfree.mtz \
  refinement.input.bulk_solvent_map.file_name=refme_minRfree.mtz \
  refinement.input.bulk_solvent_map.amplitudes_label=Fpart \
  refinement.input.bulk_solvent_map.phases_label=PHIpart \
  refinement.main.number_of_macro_cycles=1 \
  "refinement.refine.strategy=individual_sites individual_adp occupancies" \
  output.prefix=test_usermap --overwrite
```

Expected: Final R-work ~0.146, R-free ~0.156.

## Known limitations

- Joint X-ray/neutron refinement path in `validate()` is not wired up
- Twin refinement untested (user f_mask_twin is not set)
