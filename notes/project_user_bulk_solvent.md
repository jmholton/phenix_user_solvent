---
name: User-supplied bulk solvent map feature
description: Feature letting user MTZ columns (FPART/PHIPART) replace phenix's computed flat bulk solvent mask; patched into phenix-2.2-6143 (working install)
type: project
originSessionId: 388d4cb9-8cb3-48f2-8718-9e468f444302
modified: 2026-08-17T21:28:57.011Z
---
Added support for a user-supplied bulk solvent electron density map (à la Refmac FPART/PHIPART) to phenix.refine.

**WHICH INSTALL (2026-08-17): patch now lives in the working `/programs/phenix-2.2-6143/` (Python 3.11).**
The older `/programs/phenix-2.1rc2-6037/` install where it was first developed is BROKEN — after a
stale-backup restore, its cctbx/stdlib portions were "not on tape" (phenix.python there fails with
"No module named 'encodings'"). 2.2-6143 was chowned to jamesh and is the live patched target.
`/programs` is a symlink to `/home/programs` (same files). Files patched under
`/programs/phenix-2.2-6143/lib/python3.11/site-packages/`; each has a `.preusersolvent` backup.

**Usage (`phenix.refine` in PATH now resolves to the patched 2.2-6143):**
```
phenix.refine model.pdb data.mtz \
  refinement.input.bulk_solvent_map.file_name=solvent.mtz \
  refinement.input.bulk_solvent_map.amplitudes_label=Fpart \
  refinement.input.bulk_solvent_map.phases_label=PHIpart
```

**Model:** the user map REPLACES the flat 0/1 F_mask entirely. k_mask is fit per-resolution-bin
analytically (Afonine 2013), NOT the old `k_sol*exp(-B_sol*s²/4)`. Mask frozen for the whole run.
User-mask signature = k_mask stays ~1.0 and RISES with resolution; flat mask collapses to 0.00 at high res.

**Files changed (4, all in `.../python3.11/site-packages/`):**
1. `phenix/refinement/__init__.params` — added `refinement.input.bulk_solvent_map` scope
2. `mmtbx/f_model/f_model.py` — `_user_f_masks` attribute, `set_user_f_masks()`, guard in `update_xray_structure()`, and propagation through `select()` (critical bug fix)
3. `phenix/programs/phenix_refine.py` — `_inject_user_bulk_solvent()` helper called in `validate()` (single-refinement path only; joint XN not wired)
4. `mmtbx/refinement/minimization.py` — per-atom U_iso shift clamp in `apply_shifts()` (ADP pathology fix for high-copy ensembles). NOTE: only the documented ADP clamp was ported to 2.2; the surviving 2.1 file had drifted to also include an undocumented xyz-shift clamp + a `PHENIX_SHIFT_CLAMP` env toggle — deliberately NOT ported.

**Critical bug fixed:** `manager.select(in_place=True)` (called by `remove_outliers` 5× per scaling cycle) was creating a new manager with `_user_f_masks = None` then copying all attributes back, silently wiping the user mask before LBFGS minimisation. Fix: propagate `_user_f_masks` through `select()`.

**Demonstrated on evenmoreconf.pdb (48-copy ensemble solvent):**
- Default flat mask: Final R-work 0.1768
- User mask after fix: Final R-work 0.1084

**Port validated on 2.2-6143 (2026-08-17, example/starthere.pdb+refme.mtz+solvent_Fpart.mtz, 3 cyc, individual_sites+individual_adp):**
- User solvent: R-work 0.104 / R-free 0.119, k_mask 0.97→1.05 (rises, user-mask signature)
- Default solvent baseline: R-work 0.129 / R-free 0.145 (≈ documented 2.1 default 0.128/0.147), k_mask 0.31→0.00 (collapses)
- User mask beats default by ΔRwork 0.025 / ΔRfree 0.026. User absolute is ~0.01 higher than 2.1's documented 0.093/0.109 (small version-specific scaling diff; default matches 2.1 to the digit).

**Best-Rfree run (rerefine_003) re-run on patched 2.2 (2026-08-17) — definitive port validation:**
Re-ran the actual lowest-Rfree job (rerefine_002.pdb + refme.mtz + example/solvent_Fpart.mtz, its
own rerefine_003.eff: individual_sites individual_sites_real_space individual_adp occupancies, 3 cyc,
random_seed 2679941) on 2.2. Result: Start R 0.0916/0.1082 BIT-IDENTICAL to 2.1; Final 0.0874/0.1065
vs 2.1's 0.0880/0.1060 (ΔRfree +0.0005, refinement-path noise). Confirms the port is faithful on the
real 9888-atom / 45-conformer ensemble, not just the toy example. Ran in scratch with symlinked inputs
+ new prefix so the original rerefine_003.* was never overwritten.

**Why:** Bulk solvent from MD/ensemble can dramatically improve refinement vs. flat 0/1 mask for complex crystal forms.

**GitHub:** https://github.com/jmholton/phenix_user_solvent
