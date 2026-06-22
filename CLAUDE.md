# phenix_user_solvent

Feature branch adding user-supplied bulk solvent map support to phenix.refine.
Changes are applied to `/programs/phenix-2.1rc2-6037/` (writable).

**Important:** The default `phenix.refine` in PATH is version 2.0-5936 (unmodified).
Source the phenix-2.1rc2-6037 environment so that `phenix.refine` resolves to the patched version:
`source /programs/phenix-2.1rc2-6037/phenix_env.sh`

## What was changed and why

Five changes in the phenix-2.1rc2-6037 installation:

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

**`lib/python3.9/site-packages/mmtbx/refinement/minimization.py`** (one change):
- Added per-atom U_iso shift clamping in `apply_shifts()`, inside the `if(self.refine_adp)` block,
  after `xray.ext.minimization_apply_shifts()` fills `scatterers_shifted` and before
  `replace_scatterers()`. For each atom with `0 < occ < 1` and `use_u_iso()`, the total
  U_iso shift over the LBFGS session is clamped to `max_delta_u * occ`, where
  `max_delta_u = adptbx.b_as_u(5.0)` (≈0.063 Å²). Sign is preserved with `math.copysign`.

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

**Canonical user vs default solvent comparison (~5 min each):**

User solvent (log: `example/usersolvent_001.log`):
```
phenix.refine \
  example/starthere.pdb example/refme.mtz \
  refinement.input.bulk_solvent_map.file_name=example/solvent_Fpart.mtz \
  refinement.input.bulk_solvent_map.amplitudes_label=Fpart \
  refinement.input.bulk_solvent_map.phases_label=PHIpart \
  refinement.main.number_of_macro_cycles=3 \
  refinement.refine.strategy="individual_sites individual_adp" \
  output.prefix=usersolvent
```
Expected: R-work=0.093, R-free=0.109.

Default solvent (log: `example/defaultsolvent_001.log`):
```
phenix.refine \
  example/starthere.pdb example/refme.mtz \
  refinement.main.number_of_macro_cycles=3 \
  refinement.refine.strategy="individual_sites individual_adp" \
  output.prefix=defaultsolvent
```
Expected: R-work=0.128, R-free=0.147.

## ADP refinement pathology for high-copy ensembles

For multi-conformer / high-copy ensemble structures refined with a user-supplied
bulk solvent map, unmodified ADP (B-factor) refinement is harmful: b_min collapses
to 0 in the first cycle. Root cause is 1/occ gradient amplification in LBFGS.

**Mechanism:** For a partially-occupied atom, X-ray gradient ∝ occ, curvature ∝ occ²,
so the Newton-Raphson step ∝ 1/occ. For a 5-copy ensemble (occ=0.2) every LBFGS
step is 5× too large, driving b_min to the floor instantly. The through-space ADP
restraints (sphere_radius=5 Å) are insufficient to prevent this.

**Fix (minimization.py):** Per-atom U_iso shift clamped to `adptbx.b_as_u(5.0) * occ`
per LBFGS session. For occ=0.2 the per-session B-shift limit is 1.0 Å², preventing
collapse from B=1.2 to B=0 in a single macro cycle. The clamp stabilizes ADP
refinement but does not make it strictly better than skipping ADP (see table below).

Systematic controls on multiconf.pdb (5-copy ensemble, Fpart bulk solvent, 3 cycles):

| Strategy | R-free end | b_min |
|---|---|---|
| occ only | **0.1166** | 1.2 (stable) |
| xyz only = xyz + occ | 0.1177 | 1.2 (stable) |
| BSS only (zerocyc) | 0.1181 | 1.2 (stable) |
| xyz + ADP + occ, clamp | 0.1187 | 1.3–1.4 (stable) |
| B only (any wxu_scale, no clamp) | 0.1193–0.1237 | 0.0–3.5 |
| xyz + B (any wxu_scale, no clamp) | 0.1198–0.1222 | 0.0–4.1 |
| full defaults + RSR | 0.1200 | 0.0 |

The clamped xyz+ADP+occ run (log: `example/multiconf_usersolvent_xyzoB_clamp_001.log`)
stays stable over 3 cycles but remains ~0.001 worse than xyz+occ. The signal that
would drive B-factors toward the truth is apparently too weak relative to the
ill-conditioned gradient for partially-occupied conformers to show a benefit within
3 cycles at this clamp strength (`max_delta_u = adptbx.b_as_u(5.0)`).

**Recommended strategy for high-copy ensemble refinement:**
```
"refinement.refine.strategy=individual_sites occupancies"
```
Omit `individual_adp`. Use 3–10 macro cycles; xyz+occ converges stably.
The ADP clamp in minimization.py keeps the option open if tighter tuning or
more cycles eventually shows improvement.

This is also the correct strategy for generating converged local-minimum
structures as CNN training data (conformer-swap perturbation → refine → CNN
learns to detect stuck assignments from the map).

## Iterative workflow: phenix refinement + density editing

The per-bin `k_mask` in phenix exposes a limitation of the refmac-style preprocessing
pipeline that produces the user solvent map. Refmac applies a single global B-factor
to the MD solvent map to fit it to the data. This blurs the map uniformly, which:

1. Suppresses high-resolution solvent content
2. Over-weights the already-dominant low-frequency terms
3. Makes genuine ordered-solvent peaks disappear into the diffuse background

Those suppressed peaks become the largest Fo-Fc features, so the refmac pipeline
adds them back as explicit partial-occupancy coordinate waters. The model grows
more complex to compensate for the blurring artifact.

**phenix breaks this compensatory loop.** With per-bin `k_mask` free to scale each
resolution shell independently, it reveals that the refmac-preprocessed map was
globally underscaled once the low-frequency bias is released. All k_mask values
exceed 1.0 in a well-converged run:

| Resolution (Å) | k_mask (default flat) | k_mask (user, cycle 1) | k_mask (user, later) |
|---|---|---|---|
| 14–6.7 | 0.33 | 0.95 | 1.05 |
| 5.4–4.4 | 0.32 | 0.92 | 1.12 |
| 2.9–2.4 | 0.00 | 0.98 | 1.12 |
| 2.4–1.9 | 0.00 | 1.15 | 1.32 |
| 1.5–1.3 | 0.00 | 1.30 | 1.45 |

k_mask > 1.0 everywhere means the map amplitude was universally undersold by the
refmac B-factor fitting. The coordinate waters that were patching the blurred-map
deficit become less necessary, the model simplifies, and R-free improves.

**After phenix refinement with the user solvent, Fo-Fc features in the solvent
region are genuine missing structure**, not blurring artifacts. These positive-only
features (no significant negatives, because the protein model is correct) represent
solvent ordered beyond what the MD trajectory captured. They are the input for
the next density-editing cycle.

**Convergence criterion:** iterate until solvent-region Fo-Fc features fall below
~3σ. Track across cycles:
- k_mask profile: should stabilize (values stop rising)
- Number of >5σ solvent peaks: should decrease monotonically
- R-free: should decrease monotonically

**Recommended per-cycle strategy for multi-conformer ensembles:**
```
refinement.refine.strategy="individual_sites occupancies"
```
with main-chain and side-chain occupancies in separate sum-to-unity groups.
Omit `individual_adp` (see ADP pathology section above).

**Inspecting the bulk solvent contribution across cycles:**
Add `output.export_final_f_model=True` to write a `_f_model.mtz` containing
`FMASK`, `PHIFMASK`, and `K_MASK` columns. Use `wilson_b_fmask.py` to track
the Wilson B of `k_mask * F_mask` across iterations; rising k_mask and a
stabilizing Wilson B indicate convergence.

## Known limitations

- Joint X-ray/neutron refinement path in `validate()` is not wired up
- Twin refinement untested (user f_mask_twin is not set)
- ADP refinement harmful for high-copy ensembles (see above)
