# best_so_far

Lowest-R-free structure to date — first squish→phenix solvent-update cycle.

| (phenix) | R-work | R-free |
|---|---|---|
| **Final** | **0.0819** | **0.1021** |
| Start (cycle 0) | 0.0872 | 0.1056 |

Improvement over the previous best (`prev_rerefine_rerun_003/`, R-free 0.1065), confirmed on
every metric — most strongly on free LLG (the sensitive one, refmac NCYC 0 scorer):

| metric | previous | this | Δ |
|---|---|---|---|
| phenix R-free | 0.1065 | **0.1021** | −0.0044 |
| refmac R-free | 0.1124 | 0.1095 | −0.0029 |
| refmac −LL free | 4853.9 | 4815.0 | **−38.9** |
| refmac Free FSC | 0.9858 | 0.9861 | +0.0003 |

This clears the ~0.110 R-free ceiling the 1280-cycle refmac longrun plateaued at.

## Provenance — the squish→phenix loop

Starting model `rerefine_002.pdb`, data `refme.mtz`, same `.eff`/seed as the previous best;
only the bulk-solvent map differs. The solvent map (`solvent_Fpart_used.mtz`) was produced by:

1. **squish** (`../squish_solvent_runme.com`, phenix-input mode) on the previous best's
   `_f_model.mtz` → solvent-region Fo−Fc increment (`refme_Fpart.mtz`). Free reflections
   excluded from the difference map (`exclude_freeR=1`); Wilsonification `loB≥0` guard active.
2. **combine** (`../combine_solvent.py`): `new_FMASK = K_MASK·FMASK + k·increment/(K_ISO·K_ANISO)`,
   electron units, divot re-scaled `k=0.83` (fit to fofc, work set only).
3. **phenix.refine** with that map (patched user-bulk-solvent feature).

Free set never touched the solvent-map shape or the divot scale — gains are honest cross-validation.

## Files

| File | Contents |
|---|---|
| `updatedsolvent_003.pdb` / `.cif` | refined model |
| `updatedsolvent_003.mtz` | map coefficients |
| `updatedsolvent_003_f_model.mtz` | f_model (input for the next squish cycle) |
| `solvent_Fpart_used.mtz` | the squish-updated bulk-solvent map used (Fpart/PHIpart) |
| `updatedsolvent_003.log`, `.geo`, `.eff` | refinement log, geometry, params |
| `prev_rerefine_rerun_003/` | the previous best (R-free 0.1065) |

Refined 2026-08-20.
