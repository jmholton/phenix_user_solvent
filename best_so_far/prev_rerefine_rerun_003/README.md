# best_so_far

Lowest-R-free refinement to date, with the user-supplied bulk solvent map.

| | R-work | R-free |
|---|---|---|
| **Final** | **0.0874** | **0.1065** |
| Start (cycle 0) | 0.0916 | 0.1082 |

## Provenance

This is `rerefine_003` re-run on the **patched phenix-2.2-6143** (user-solvent feature).
It reproduces the original 2.1rc2-6037 best-R-free run (Final 0.0880 / 0.1060) to within
noise; the start R is bit-identical, confirming the port is faithful on the real
9888-atom / 45-conformer ensemble.

- Input model: `rerefine_002.pdb`
- Data: `refme.mtz`
- Bulk solvent map: `solvent_Fpart.mtz` (Fpart / PHIpart) — replaces the computed flat mask
- Strategy: `individual_sites individual_sites_real_space individual_adp occupancies`
- 3 macrocycles, `n_gaussian`, `random_seed=2679941`
- Full parameters: `rerefine_003.eff` (input) / `rerefine_rerun_003.eff` (effective)
- Refined: 2026-08-17

## Files

| File | Contents |
|---|---|
| `rerefine_rerun_003.pdb` / `.cif` | refined model |
| `rerefine_rerun_003.mtz` | map coefficients |
| `rerefine_rerun_003_f_model.mtz` | f_model: FMODEL, FCALC, FMASK, K_ISOTROPIC, K_ANISOTROPIC, K_MASK (input for `wilsonify.py`) |
| `rerefine_rerun_003.log` | refinement log |
| `rerefine_rerun_003.geo` | geometry restraints |
| `rerefine_rerun_data.mtz` | processed input data |
