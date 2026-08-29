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

**Strategy syntax — use `+`, not a quoted space.** Multiple refinement strategies MUST be
joined with `+` (`strategy=individual_sites+individual_adp`). The quoted-space form
`strategy="individual_sites individual_adp"` becomes a single argv token that phil parses
as one (invalid) choice → *empty* selection → **0 atoms refined** (only bulk-solvent
scaling runs, R stays flat at the input value). It fails silently: no error, and the log
shows `individual_sites = False (0 atoms)`. Sanity-check any run with
`grep "individual_sites .*=" *.log` — it must read `True (N atoms)`.

User solvent (log: `example/usersolvent_001.log`):
```
phenix.refine \
  example/starthere.pdb example/refme.mtz \
  refinement.input.bulk_solvent_map.file_name=example/solvent_Fpart.mtz \
  refinement.input.bulk_solvent_map.amplitudes_label=Fpart \
  refinement.input.bulk_solvent_map.phases_label=PHIpart \
  refinement.main.number_of_macro_cycles=3 \
  refinement.refine.strategy=individual_sites+individual_adp \
  output.prefix=usersolvent
```
Expected: R-work≈0.095, R-free≈0.111 (current patch, with the ADP occ-clamp active).
The reference `example/usersolvent_001.log` value of R-work=0.0928, R-free=0.1092 predates
the clamp: the clamp caps per-session ADP excursions (final b_max ≈39 vs ≈50 unclamped),
which costs ~0.002 R here. Removing the clamp reproduces 0.0928/0.1092 exactly.

Default solvent (log: `example/defaultsolvent_001.log`):
```
phenix.refine \
  example/starthere.pdb example/refme.mtz \
  refinement.main.number_of_macro_cycles=3 \
  refinement.refine.strategy=individual_sites+individual_adp \
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
refinement.refine.strategy=individual_sites+occupancies
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
refinement.refine.strategy=individual_sites+occupancies
```
with main-chain and side-chain occupancies in separate sum-to-unity groups.
Omit `individual_adp` (see ADP pathology section above).

**Inspecting the bulk solvent contribution across cycles:**
Add `output.export_final_f_model=True` to write a `_f_model.mtz` containing
`FMASK`, `PHIFMASK`, and `K_MASK` columns. Use `wilson_b_fmask.py` to track
the Wilson B of `k_mask * F_mask` across iterations; rising k_mask and a
stabilizing Wilson B indicate convergence.

## phenix vs refmac R-factors with an identical external solvent map

Comparison study: fixed model (`rerefine_003.pdb`) and the SAME k_mask-scaled solvent
handed to both programs (to refmac as an additive partial structure, `FPART1=Fpart
PHIP1=PHIFpart`, `SOLVENT NO`). The Fpart column equals phenix's `K_MASK * FMASK`
exactly (ratio 1.0000), so the per-bin solvent is already baked into refmac's input.

**Program-self-reported R-factors (not recomputed):**

| Run | R-work | R-free | source |
|---|---|---|---|
| phenix (per-bin scaling) | 0.088 | 0.106 | `phenix_usersolvent/rerefine_003.log` ("Final R-work = 0.0880, R-free = 0.1060") |
| refmac NCYC 0, Fobs as-is | 0.099 | 0.117 | refmac stdout |
| refmac NCYC 0, rescaled Fobs — scored vs *rescaled* target | 0.098 | 0.111 | refmac stdout |
| refmac full refine (many cycles) | 0.086 | 0.117 | `refmac_phenixsolvent/converge_ramp_003.log:122152` |

Refmac R-work can match phenix after full coordinate refinement (0.086), but its
R-free never reaches 0.106 — it stays 0.111–0.117.

**The rescaled-Fobs result does NOT actually help — its lower R is a scoring artifact.**
`refme_rescaled.mtz` divides Fobs by phenix's k_iso·k_aniso; refmac then reports
R=0.098/0.111 *against that rescaled target*, which is a different (pre-softened)
denominator and is NOT comparable to phenix's 0.106. Scored fairly — back-scale
refmac's rescaled-fit model to the original scale with the exact per-reflection factor
`FP_orig/FP_rescaled` (which round-trips the data exactly) and compare to the real
Fobs — it gives **0.112/0.120**, i.e. *worse* than refmac's direct fit (0.098/0.114 on
the common set). So at fixed geometry the rescaling is a detour, not an improvement.
(An earlier note claiming rescaling "closes ~half the gap" was scored against the
rescaled target and is wrong.) Whether rescaling helps during actual coordinate
refinement (NCYC>0) is a separate, untested question.
UPDATE: this pessimism was about that EARLY attempt (extrapolated scales, no SHANNON, and
scored/back-scaled awkwardly). Done PROPERLY -- exact per-reflection K_ISO*K_ANISO from the
f_model + physical Fpart + SHANNON 6 -- "Wilsonification" reproduces phenix R-free exactly
(0.106); see the Wilsonification section below.

**Reflection-set control:** the phenix↔refmac gap is not an artifact of different
reflection sets. On a common set (phenix's 30985), R-free is phenix 0.106 vs refmac
0.114; phenix without its 14-outlier rejection is 0.089/0.108, refmac on phenix's set
is 0.098/0.114 — the ~0.008–0.010 R-free gap is stable whichever common set is used.

**The refmac log's "Overall B = 5.05" is a Wilson-scaling DIAGNOSTIC, not a
transformation applied to the model.** (An earlier analysis wrongly claimed this B
damps the model and undoes the high-resolution `k_mask`; that is corrected here.)
Evidence:
- refmac bare Fcalc (`MODE SFCALC`) matches gemmi Fcalc to relative B = −0.04, flat
  per shell — a flat overall scale offset (~5%) that cancels in R and is NOT form factors
  (Cromer-Mann vs n_gaussian/it92 differ only 0.08%; later shown FCALC_R = FCALC_P under
  SHANNON). No 5 Å² lives in the structure factors.
- refmac's *effective* applied envelope (output `|FC_ALL_LS|` / bare `|model+Fpart|`)
  is FLAT at ~1.00 ±2% across all resolution, not `exp(−5·s²/4)` (which would fall to
  0.27 at 0.96 Å).
- a literal B=5 on `(model+Fpart)` gives R=0.33 (single-B optimum is B≈0); refmac
  reports 0.099 — so the 5 is not applied.
- the reported B is unchanged by `make hydr no` (→5.18) or by dropping `SIGFP`
  (→identical 5.0459); it is a property of the data's Wilson plot, which is why it is
  "locked" and does not move when the solvent map is perturbed.

**Over-sharpening the solvent does not help.** Because the applied envelope is flat and
non-adaptive, pre-sharpening Fpart by `exp(+B·s²/4)` is not re-damped — it adds
high-resolution solvent power the data lacks, so R worsens (over-sharpen +5 → 0.133;
plain → 0.099; a mild blur −5 → 0.093, still short of phenix and just hand-fitting a
second global Gaussian).

**Real source of the residual ~0.01 R-free gap:** not an overall B. It is the finer
scaling/target differences — refmac's ML/σA (Rice) target and per-bin D scaling
(`k_mask` clamped ≤1.5, `b_ls_part` clamped ≥0, blur-only) vs phenix's unweighted
per-bin R-minimization. (Form factors are NOT a factor -- 0.08%.) No refmac keyword closes it:
`SCAL MLSC FIXB B 0.0` is dead code (`MLBFIX`/`BML_DEFINED` are parsed at
`rcard_tor1.f:5226-5227` but never read by any scaling routine, and
`B_LS_OVER_REFINE_FLAG` is never set false). Matching phenix would require patching
refmac to accept per-bin scale arrays — the refmac analog of the phenix
`set_user_f_masks` patch. Note this is distinct from the solvent-map *preprocessing*
pipeline discussed above, where refmac's partial-structure B (`b_ls_part`, genuinely
applied and clamped ≥0) does blur the MD map.

### Reconstructing the logged R-factors from output MTZ columns (validation recipe)

Both programs' logged R-factors are reproducible from their output MTZ columns to the
last digit. The traps are the free-flag convention and which Fc column refmac reports.

**Phenix** (`_f_model.mtz`): R = Σ|FOBS − |FMODEL|| / Σ|FOBS|, **k=1** (FMODEL is the
final scaled model, already carrying k_iso·k_aniso·k_mask). Free set = `R_FREE_FLAGS
== 1`, work = `== 0`. Reproduces 0.0880 / 0.1059 vs logged 0.0880 / 0.1060.

**Refmac** (output `.mtz`): R = Σ|FP − |FC_ALL_LS|| / Σ|FP|, **k=1** (columns already
scaled). Free set = `FreeR_flag == 0`, work = any other value. Reproduces 0.0990 /
0.1164 vs logged 0.0992 / 0.1167. Use **`FC_ALL_LS`** (= `FC`; the unweighted total
Fcalc+Fpart with overall scale) — NOT `FC_ALL`, the σA/ML-weighted column, which gives
0.098/0.113 and does not match the log.

**Free-flag convention trap:** phenix `_f_model.mtz` uses a binary `R_FREE_FLAGS`
(free = 1). Refmac's `refme.mtz` uses a 20-bin CCP4 `FreeR_flag` (0–19; free = the
single value **0**). They flag the *same* reflections (100% agreement) but with
opposite encoding. Selecting the wrong one silently computes R-work and calls it
R-free (this produced a spurious "refmac Rfree=0.106" earlier).

**Reflection accounting:** `refme.mtz` is a complete set — 34623 HKL rows to d_min but
only 31001 with measured FP (iotbx drops the MNF rows). Observed = 29398 work + 1603
free = 31001, matching refmac's "29398 used". The log's "1767 free" counts flag==0 in
the full generated set (includes 164 unobserved that never enter R). Refmac rejects no
outliers (uses all 29398/1603). Phenix rejects exactly 14 (its `remove_outliers`) and
cuts low-res at 14.12 Å vs refme's 15.82 Å, giving 30985 — which is why phenix's
FMODEL set and refmac's set differ slightly and each log matches only its own set.

### Output MTZ columns, program by program

MTZ type codes: F=amplitude, Q=sigma, P=phase(°), I=integer, A=Hendrickson-Lattman
coeff, R=real, W=weight, H=Miller index. All relationships below were verified from the
column data, not asserted.

**Phenix** `rerefine_003_f_model.mtz` (the `export_final_f_model` output; 30985
observed post-rejection reflections; contains model + scales, NOT map coefficients):

| Column(s) | Type | Meaning |
|---|---|---|
| H K L | H | Miller indices |
| FOBS / SIGFOBS | F/Q | observed amplitude (French–Wilson) and σ |
| R_FREE_FLAGS | I | **binary: 1=free (5%), 0=work** |
| FMODEL / PHIFMODEL | F/P | total scaled model = `K_ISO·K_ANISO·(FCALC + K_MASK·FMASK)` (verified to 1.2%); phenix's reported R uses this (0.088/0.106 vs FOBS) |
| FCALC / PHIFCALC | F/P | bare atomic structure factor (no scale, no solvent); matches gemmi to R=0.009 |
| FMASK / PHIFMASK | F/P | bulk-solvent mask F — here the **user-supplied solvent map**, unscaled |
| K_ISOTROPIC | F | per-refl isotropic scale (smooth exp × per-bin free): 0.872–1.063, non-monotonic |
| K_ANISOTROPIC | F | per-refl anisotropic scale (6-param tensor): 0.695–1.106 |
| K_MASK | F | per-refl solvent scale on FMASK: 1.010–1.500, rises with resolution |
| SIGK_MASK | Q | uncertainty of K_MASK |
| FOM | F | figure of merit m (0–1) |
| ALPHA / BETA | F/F | ML params: α = D/σA scale on Fc (0.92–0.99); β = variance (23–564, grows with res) |
| CENTRIC_FLAGS | I | 1=centric, 0=acentric |
| HLmodelA–D | A | Hendrickson–Lattman coeffs for the model phase probability |
| HLcombA–D | A | HL coeffs for combined phases (= model here; no experimental phases) |
| RESOLUTION | R | d-spacing (Å) per reflection |

**Refmac** `refme_orig_ncyc0.mtz` (complete set of 34612 reflections for map coeffs;
only 31001 have an observed FP):

| Column(s) | Type | Meaning |
|---|---|---|
| H K L | H | Miller indices |
| FreeR_flag | I | **20-bin CCP4 flag (0–19); free = value 0 (5%)** |
| FP / SIGFP | F/Q | observed amplitude (refmac's copy, ~0.991× input FP) and σ |
| FC / PHIC | F/P | calculated total **including partial structure (Fpart)** with overall LS scale — numerically identical to FC_ALL_LS (diff 0.0000) |
| FC_ALL / PHIC_ALL | F/P | **σA/D-weighted** total (ML); 0.970× FC_ALL_LS; used for ML target + map coeffs |
| FWT / PHWT | F/P | **2mFo−DFc** main-map coefficient (corr 0.994 with 2·FP−FC_ALL) |
| DELFWT / PHDELWT | F/P | **mFo−DFc** difference-map coefficient |
| FOM | W | figure of merit m |
| FC_ALL_LS / PHIC_ALL_LS | F/P | **unweighted** LS total (Fcalc+Fpart, overall scale, no D-weight) = FC; **refmac's reported R uses this** (0.099/0.117 vs FP) |

**Cross-program gotchas:**
- Free flag encoding is opposite: phenix `R_FREE_FLAGS==1` = free; refmac `FreeR_flag==0`
  = free. Same reflections, inverse convention.
- Reported R uses different columns: phenix → FMODEL; refmac → FC_ALL_LS (=FC), NOT the
  σA-weighted FC_ALL (which gives 0.098/0.113). Both at k=1.
- phenix's FMASK is the user solvent map; the Fpart handed to refmac = K_MASK·FMASK exactly.
- phenix exposes K_ISOTROPIC / K_ANISOTROPIC / K_MASK as explicit per-refl columns;
  refmac bakes its overall scale into FC/FC_ALL/FC_ALL_LS and reports only the summary
  "scale, B" in the log.
- phenix file = observed post-rejection set (30985); refmac file = complete set (34612)
  so map coeffs exist for unmeasured reflections, FP present only for the 31001 observed.

### Equations relating refmac columns to phenix columns

Notation: `_P` = phenix column, `_R` = refmac column; per-reflection. All coefficients
below were fitted from the column data, not asserted.

Bare atomic Fcalc (SHANNON grid; protein Fcalcs AGREE between programs):
```
FCALC_R ~= FCALC_P               (shell-mean ratio 1.00 at all resolution; R=0.009
                                  overall = per-reflection ensemble scatter, averages to 1)
```
CORRECTION: an earlier note claimed FCALC_R = 1.045*FCALC_P "form-factor scale". That was
WRONG -- the 1.045 was the b_add_loc GRID artifact of a default-grid MODE SFCALC run. With
SHANNON 4-6 the grid artifact is gone and the two protein Fcalcs are equal. The tables the
programs actually use (phenix n_gaussian, refmac Cromer-Mann/it1992) differ by only 0.08%
(R), <0.05% per shell even at 0.96 A -- form factors are NEGLIGIBLE, not a cause of anything.
Bulk solvent (exact; refmac was fed phenix's already-scaled mask):
```
Fpart_R = K_MASK_P * FMASK_P     (complex, exact ratio 1.0000)
```
Phenix total model (verified to 1.2%):
```
FMODEL_P = K_ISO_P * K_ANISO_P * ( FCALC_P + K_MASK_P * FMASK_P )
```
Refmac total model = the column its reported R uses (FC == FC_ALL_LS):
```
FC_ALL_LS_R = FC_R = s_R * ( FCALC_R + Fpart_R )
                   = s_R * ( FCALC_P + K_MASK_P*FMASK_P )   (FCALC_R = FCALC_P, SHANNON grid)
```
  s_R = refmac's flat overall scale (~1.0); the log's "B=5" is a Wilson diagnostic,
  NOT applied (effective envelope flat to +/-2%).

Bridge between the two totals (same shape, flat scale apart):
```
FMODEL_P ~= 0.933 * FC_ALL_LS_R      (relative B = +0.13 ~ 0)
```
Side by side, this isolates the entire source of the R gap:
```
  phenix:  [ K_ISO*K_ANISO ] * ( FCALC + K_MASK*FMASK )   per-bin env (0.70-1.13)
  refmac:  [   0.933*s_R    ] * ( FCALC + K_MASK*FMASK )   flat env
```
Same FCALC (form factors negligible, 0.08%); the ONLY difference is phenix's per-bin
K_ISO*K_ANISO envelope vs refmac's flat scale. That is the whole ~0.008 R-free gap.
sigmaA/likelihood weighting (refmac's D-factor equals phenix's ALPHA):
```
FC_ALL_R = D_R * FC_ALL_LS_R ,   D_R(bin) ~= ALPHA_P(bin)
  per shell D_R  = 0.99 1.01 1.01 0.99 0.95 0.89
  per shell ALPHA= 0.98 0.99 0.98 0.98 0.95 0.92
```
Figure of merit and map coefficients (FC_ALL is already D-weighted):
```
FOM_R ~= FOM_P = m
FWT_R    = 2m*|FP| - |FC_ALL_R|  (phase PHIC_ALL)   -> 2mFo-DFc  (corr 0.994)
DELFWT_R =  m*|FP| - |FC_ALL_R|  (phase PHIC_ALL)   -> mFo-DFc
```
Observed amplitudes (same intensities, minor per-program handling):
```
FP_R  ~= 0.991 * FP_input ;   FOBS_P = French-Wilson(I)
```

### Data flow through refmac, and the atomic-Fcalc drift (why FC never adds up exactly)

Pipeline (single external partial, no SCPART, NCYC 0), with file:line anchors in
`/programs/ccp4-8.0-src/checkout/refmac5.8/`:

```
HKLIN --MTZRED(oppro_allocate.f:4071)--> FP scaled by scale_sharp*exp(b_sharp) [no-op];
        FPART1/PHIP1 stored RAW, no scale/no B (:4337). NPART:=NSCPART=0, NPARTALL=1 (:4653)
XYZIN --REF_ALL(hkon_secder_tch.f:4690)--> atomic Fcalc, FFT-based:
        (a) if min atomic B too small for the grid, inflate EVERY atom's U by u_add_loc
            (:4808-4826), build gridded Gaussian density (subag_scale_hessian.f ELDEN),
            FFT (subvag1.f:471);
        (b) remove that grid B analytically per reflection: fc *= exp(+rsq*b_add_loc) (:5012)
        --> atomic Fcalc is APPROXIMATE (gridding residual)
RDIND1(subag_scale_hessian.f:1079): atomic Fcalc -> ICPOS=IREF+NPART*NOBS block;
        external partial ADDED AT UNIT SCALE: FC(ICPOS)+=FPP*cos(APP) (:1212) --> EXACT
LS/ML scaling: fit SCALE_LS_OVER, B_LS_OVER, aniso; applied to FOBS (FO/(scale*aniso))
        and to the A_ALL/B_ALL accumulators; overall B is also folded into the per-atom
        exponent for the SF/variance calc and removed analytically
MTZWRITE(oppro_allocate.f:5580):
   oc_FC        = FC(IR+NPART*NOBS,1) = |atomic Fcalc + Fpart|, UNIT scale (:5995)
   oc_FC_ALL_LS = sqrt(a_all_ls^2+b_all_ls^2)  (unscaled sum, ~= oc_FC)
   oc_FC_ALL    = sqrt(A_ALL^2+B_ALL^2)  (sigmaA/D per-bin scaled)
   oc_FWT/DELFWT/FOM = 2mFo-DFc / mFo-DFc / m
```

Empirically confirmed: `oc_FC` is NOT damped by the overall B (model-only FC / phenix
FCALC is flat 1.00 across resolution). The reported "scale/B" go on Fobs, not on oc_FC.

**The drift = imperfect cancellation of an added/removed B in the atomic FFT path.**
refmac adds a B-factor to the atoms for its internal computation (the grid `b_add_loc`
and/or the fitted overall `B_LS_OVER` folded into the per-atom exponent) and removes it
analytically. The two do not cancel exactly (gridded FFT vs analytic exp), leaving a
resolution-graded residual whose size scales with |B correction|. The external partial
has NO such path (added verbatim), so it is exact.

**Root trigger: phenix's very small atomic B-factors force too coarse an FFT grid.**
rerefine_003.pdb has minB=1.23, median 6.3 (rerefine_002.pdb: minB=0.96; e.g. a water O
at B=0.96 -- physically implausible, likely a phenix over-sharpening/b_min-collapse
artifact). Atoms this sharp fall below refmac's default-grid minimum, triggering
`b_add_loc` inflation; the imperfect analytic removal leaves the drift AND shows up as a
spurious ~+5 A^2 overall B. Verified two ways:

Bumping all B (drift tracks |the B refmac must add/remove|):
| all-atom B shift | refmac overall B fit | FC(Fpart=model)/FC(model) |
|---|---|---|
| 0 (as refined)  | +5.2  | 1.91-2.01  (~2.5% drift) |
| +5              | -4.7  | 1.91-2.01  (~2.5%, \|B\|~5) |
| +20             | -20.0 | 1.58-1.80  (~20%, \|B\|~20) |

Finer grid (the real fix -- `SHANNON 3.0`, default 1.5):
| grid | refmac overall B fit | FC(Fpart=model)/FC(model) |
|---|---|---|
| default (SHANNON 1.5) | +5.2 | 1.91-2.01 (~2.5% drift) |
| SHANNON 3.0           | +0.4 | 2.000-2.004 (drift GONE) |

So the ~+5 A^2 overall B was NOT a real sharp-model-vs-data mismatch -- it was refmac's
grid-coarseness (`b_add_loc`) compensation, and a finer grid removes both the overall B
and the Fcalc drift at once. The whole effect is a grid artifact set off by the tiny
B-factors.

**Fixes:**
- refmac: add `SHANNON 3.0` (finer FFT grid; larger factor = finer = smaller spacing).
  Confirmed to cut the drift from ~2.5% to <0.3% and drop the spurious overall B to ~0,
  without altering the model. This is the clean fix.
- phenix: no hard per-atom floor for individual ADP refinement (only
  `group_b_iso.b_iso_min=1.0` for grouped B, plus a global `b_iso_max` cap). Floor the
  model instead with pdbtools:
  `phenix.pdbtools model.pdb modify.selection="bfactor < 5" modify.adp.set_b_iso=5`.

Consequence (default grid): an externally computed Fpart cannot reproduce refmac's FC
below the few-percent gridding floor, and `Fpart = FMODEL - FC_model` constructions fail
/ diverge because they assume the atomic Fcalc cancels exactly between runs -- it does
not, until the grid is refined. NOTE: the user's refmac runs use rerefine_002.pdb; this
repo's analysis used rerefine_003.pdb -- different iterations, ~6% different FC, same
mechanism.

### Making refmac reproduce phenix's R-free exactly (derived Fpart)

Goal: MTZ columns that make refmac's reported R equal phenix's. R = Sum|Fobs-|Ftotal||/
Sum|Fobs|, so two programs match iff (same reflections) they use the same Fobs AND the
same |Ftotal|.

- phenix total:  FMODEL = K_ISO*K_ANISO*(FCALC_P + K_MASK*FMASK)   (verified exact)
- refmac total (SHANNON fine grid, no SCPART, SOLVENT NO):  FC = FCALC_R + Fpart
  (unit-scale add, no envelope on FC; FCALC_R = refmac's own Cromer-Mann atomic Fcalc)

Set FC = FMODEL and solve for the partial column:
```
   Fpart  =  FMODEL - FCALC_R        (complex; PHIFpart = its phase)
          =  (K_ISO*K_ANISO - 1)*FCALC_P  +  K_ISO*K_ANISO*K_MASK*FMASK
             \___ per-bin envelope on protein ___/   \____ scaled solvent ____/
```
The coefficient on FCALC is ~(K_ISO*K_ANISO - 1), NOT a form-factor scale: FCALC_R = FCALC_P
(protein Fcalcs agree once the grid is fixed with SHANNON; form-factor tables differ 0.08%).
So the model-envelope correction is purely phenix's per-bin k_iso*k_aniso deviating from 1
(range 0.70-1.13). The partial is NOT physical solvent -- it packs phenix's per-bin
k_iso*k_aniso and k_mask that refmac's flat scaling cannot reproduce.

Three conditions (all now met):
1. FCALC_R must be reproducible between the model-only run (where it's measured) and the
   Fpart run (where refmac recomputes it) -> requires SHANNON 4-6. Without it the grid
   drift (~2.5%) breaks the cancellation; this is why an earlier `Fpart=FMODEL-FC_model`
   attempt gave 0.115/0.126 and diverged.
2. No SCPART (partial added at unit scale; else refmac re-scales it).
3. refmac's residual overall scale re-fit is flat (k~1, B~0 under SHANNON) -> cancels in
   R's scale-invariance.

Recipe (two refmac passes + one build):
1. refmac SHANNON 6 model-only (no Fpart) -> FCALC_R (its FC column).
2. Fpart = FMODEL - FCALC_R (complex); write MTZ with real FP/SIGFP/FreeR_flag + Fpart.
   (FMODEL from phenix output.export_final_f_model=True)
3. refmac SHANNON 6 with that Fpart:
   LABIN FP=FP SIGFP=SIGFP FREE=FreeR_flag FPART1=Fpart PHIP1=PHIFpart
   SOLVENT NO ; SHANNON 6 ; NCYC 0

CONFIRMED (rerefine_003, refme.mtz):
| run | R-work | R-free | overall B |
|---|---|---|---|
| refmac, derived Fpart, SHAN 6 | 0.0898 | 0.1067 | +0.03 |
| phenix (rerefine_003)         | 0.0880 | 0.1060 | -- |
| same construction WITHOUT shannon | 0.115 | 0.126 (diverged) | +5 |

refmac reproduces phenix to dR-free = 0.0007 (dR-work 0.0018) -- from a 0.010-0.016 gap
down to noise. Residual is refmac's flat scale re-fit + slightly different target.

**This is a fixed-geometry TRANSFER (validator/round-trip), not a refinement setup.** The
moment coordinates move (NCYC>0), FCALC_R changes and Fpart=FMODEL-FCALC_R goes stale. To
refine in refmac with phenix-equivalent scaling you'd need per-cycle Fpart updates or a
refmac per-bin-scaling patch (the refmac analog of the phenix set_user_f_masks patch).

Scorecard for transferring phenix scaling into refmac:
| phenix term | transfer | works? |
|---|---|---|
| k_mask (per-bin solvent) | bake into Fpart as K_MASK*FMASK | yes (ratio 1.0000) |
| grid/b_add_loc artifact | SHANNON 4-6 or Blim 2 | yes (0.122->0.116) |
| k_isotropic (per-bin total) + everything else | Fpart = FMODEL - FCALC_R + SHANNON | yes, exact at fixed geometry |
| Wilsonification: Fobs/(K_ISO*K_ANISO) + physical Fpart + SHANNON | data rescale | YES: R-free exact (0.106); PREFERRED (physical solvent, survives NCYC>0) |
| (early rescale_fobs.py: extrapolated scales, no SHANNON) | -- | failed -- fixed by using exact scales + SHANNON |

### Derived Fpart vs original (physical) Fpart -- what actually differs

Cross-program transfer is a core goal of this project (move solvent/scaling between refmac,
phenix, and others). Two Fpart columns both make refmac work but mean very different things:
- ORIGINAL Fpart = K_MASK*FMASK: the physical bulk solvent (user map, k_mask-scaled).
- DERIVED  Fpart = FMODEL - FCALC_R: solvent PLUS a non-physical correction that packs
  phenix's per-bin k_iso*k_aniso (form factors are negligible, 0.08%).
  (FCALC_R = refmac's own bare atomic Fcalc from a SHANNON model-only run -- no solvent,
  no scaling; the exact thing refmac recomputes and adds Fpart to.)

Amplitudes: derived mean 7.47 vs original 6.71 (1.11x); phases differ 15.8 deg mean.
The difference is RESOLUTION-DEPENDENT (Wilson falloff):

| res (A)   | <\|Fpart\|> orig | <\|Fpart\|> derived | ratio |
|---|---|---|---|
| 14-2.2    | 40.4 | 39.8 | 0.99  (identical: both are the solvent) |
| 1.35-1.25 | 1.90 | 2.67 | 1.41 |
| 1.06-1.01 | 0.24 | 1.51 | 6.2  |
| 1.01-0.97 | 0.15 | 1.41 | 9.5  (derived carries a high-res model tail) |

effective Wilson B:  original 50.4 A^2 (steep -> diffuse low-res solvent);
                     derived  26.7 A^2 (shallow -> carries atomic-model high-res content).

- At LOW resolution the two AGREE -- the derived Fpart's low-res content IS the real
  solvent. So a PHYSICAL solvent transfers cleanly between programs (use K_MASK*FMASK).
- At HIGH resolution they diverge (up to ~10x): the derived tail is NOT solvent (FMASK->0),
  it is ~100% PROTEIN = (K_ISO*K_ANISO - 1)*FCALC_protein -- phenix's per-bin envelope times
  the protein Fcalc. TESTED: not form factors (protein FCALC_P = FCALC_R to shell-mean 1.00;
  n_gaussian vs Cromer-Mann only 0.08%) and not the model Fcalc (the programs agree) -- it
  is purely the per-bin k_iso*k_aniso deviating from 1 (range 0.70-1.13). It is program-pair-
  and geometry-specific.
- Transfer guidance: use ORIGINAL Fpart (K_MASK*FMASK) to move a physical solvent model to
  ANY program; use DERIVED Fpart only to make refmac exactly reproduce phenix's score at
  fixed geometry (validation) -- it is not physical and goes stale when coordinates move.

### Wilsonification: transferring phenix's scaling to refmac via Fobs (THIS WORKS)

Once we established the whole gap is phenix's per-bin K_ISO*K_ANISO envelope (NOT Fcalc,
form factors, or model), the clean transfer is to divide Fobs by that envelope
("Wilsonification" -- it straightens the Wilson plot of Fobs) and hand refmac the PHYSICAL
solvent. CONFIRMED to reproduce phenix:
```
  Fobs' = FOBS / (K_ISO*K_ANISO)          (EXACT per-reflection values from the f_model)
  SIGFP'= SIGFOBS / (K_ISO*K_ANISO)
  Fpart = K_MASK*FMASK                    (the real bulk solvent -- physical, transferable)
  refmac: SOLVENT NO ; SHANNON 6 ; NCYC 0
```
| run | R-work | R-free | refmac fit |
|---|---|---|---|
| refmac Wilsonified + physical Fpart + SHAN6 | 0.0892 | 0.1060 | scale 0.99, B +0.04 |
| phenix                                       | 0.0880 | 0.1060 | -- |

R-free matches to the digit; R-work off 0.0012 (a 1/(K*K) reweighting). refmac's fitted
scale is flat (0.99, B~0) because it inherited phenix's per-bin envelope THROUGH the data.
This BEATS the derived-Fpart trick: the solvent stays physical/transferable, and because it
is a data transformation (not a geometry-specific patch) it should survive NCYC>0 (refresh
K_ISO/K_ANISO periodically as the model moves).

Why it works now and rescale_fobs.py did not: (1) EXACT per-reflection K_ISO*K_ANISO from
the f_model, not extrapolated per-bin k_iso + a fitted aniso tensor; (2) SHANNON 6 so
refmac does not re-impose the b_add_loc grid B that fought the rescaling.

The tool: `wilsonify.py <phenix_f_model.mtz> [out.mtz]` builds the Wilsonified MTZ
(Fobs/(K_ISO*K_ANISO), rescaled sigmas, physical Fpart=K_MASK*FMASK, refmac-convention
free flag) and prints the LABIN line.

HOLDS UNDER REFINEMENT (NCYC>0). refmac 3-cycle refinement (SHANNON 4, Blim 2, wgt 0.5),
Wilsonified vs original data (both with physical Fpart):
| cycle | Wilsonified Rwork/Rfree | Original Rwork/Rfree |
|---|---|---|
| 0 | 0.0892/0.1060 | 0.0925/0.1140 |
| 1 | 0.0887/0.1054 | 0.0923/0.1133 |
| 2 | 0.0889/0.1048 | 0.0924/0.1126 |
| 3 | 0.0896/0.1062 | 0.0929/0.1136 |
Wilsonified stays at phenix-level R-free (~0.105) through refinement; original stays
~0.008 worse (~0.113) -- the per-bin scaling advantage is maintained as coordinates move,
not just at NCYC 0. Geometry improves identically (bonds 0.012->0.007). Both keep flat
overall B (transfer holds). Slight cycle-3 uptick = geometry weight + onset of scale
staleness; for long runs refresh K_ISO/K_ANISO every few cycles (re-export f_model,
re-wilsonify).

Aside -- the per-HKL "Fcalc discrepancy" that earlier made this look impossible was a
PDB-vs-CIF serialization artifact, not real Fcalc error (kept here as a caution):
refmac's Fcalc = gemmi's = EXACT (ratio 1.000-1.001), and PHENIX's FCALC appeared to be the
outlier -- it seemed to deviate from exact bare-atomic by +/-2-6% per reflection, in BOTH
directions, even at low resolution:

| HKL | gemmi(exact) | phenix | refmac | phenix/gemmi | refmac/gemmi |
|---|---|---|---|---|---|
| (0,6,0) 6.78A | 285.75 | 303.75 | 286.06 | 1.063 | 1.001 |
| (1,4,0)       | 529.76 | 515.93 | 529.93 | 0.974 | 1.000 |
| (4,3,0)       | 192.21 | 179.99 | 192.30 | 0.936 | 1.000 |

This is NOT form factors (at low res all tables give ~Z; deviation is bidirectional, not a
smooth scale) and NOT an FFT-accuracy issue (`algorithm=direct` changes nothing -- cctbx
fft = direct = gemmi). CONFIRMED cause: **PDB-vs-CIF precision loss in a 45-conformer
ensemble.** rerefine_003 is a 45-altloc ensemble, EVERY atom partial occupancy (0.00-0.54).
The .pdb, .cif and _f_model.mtz are one job (same timestamp), so it is the SAME model --
but the PDB serialization is lossy: occupancy stored to 2 dp, coords to 3 dp, vs CIF's
3 dp / 5 dp. FCALC(0,6,0) with n_gaussian:
```
  PDB  285.7   CIF  301.5   deposited f_model  303.75    (gemmi from PDB = 285.75)
```
The CIF reproduces the f_model (<1%); the PDB is the outlier. Mechanism: hundreds of
low-occupancy alt-atoms round down (480 atoms hit occ==0.00 in the PDB vs 339 in the CIF)
and the rest redistribute, shifting sensitive low-res reflections ~6%. Scattering table
was n_gaussian (from the .eff), not wk1995.

FIX: feed refmac the CIF, not the PDB (`refmac5 XYZIN rerefine_003.cif`), or use the
f_model directly. Consequence: part of the refmac-vs-phenix Fcalc gap we chased was
self-inflicted -- refmac read the rounded PDB (285.8) while phenix's FMODEL used the
full-precision ensemble (303.75). This is a general gotcha for multi-conformer ensembles:
the PDB format cannot faithfully serialize many low-occupancy conformers.

CORRECTED: that apparent Fcalc "discrepancy" was the PDB-vs-CIF serialization, not real
Fcalc error -- the protein Fcalcs actually AGREE (FCALC_P = FCALC_R). So rescaling Fobs by
K_ISO*K_ANISO is NOT compensating for a Fcalc difference; it is transferring the ONLY real
difference (the per-bin envelope), which is exactly why Wilsonification works (top of this
section). The earlier rescale_fobs.py failed only because it used extrapolated scales and
no SHANNON -- with exact per-reflection K_ISO*K_ANISO + SHANNON 6 it reproduces phenix.

Two working transfers, use whichever fits:
- Wilsonification: Fobs/(K_ISO*K_ANISO) + physical Fpart K_MASK*FMASK + SHANNON.
  Physical, survives refinement (refresh scales as model moves). PREFERRED.
- Derived Fpart: Fpart = FMODEL - FCALC_R + SHANNON. Fixed-geometry validator only.

### Is any of phenix's R-free a gridding artifact? No -- the grid artifact was refmac's

phenix's FFT Fcalc is accurate at its default grid (grid_resolution_factor=1/3):
fft = direct = gemmi, verified per reflection. So phenix's R-free (0.106) is genuine,
NOT inflated by a phenix gridding artifact. (A coarse phenix grid, 1/2, does introduce
~6% Fcalc error, but the default is fine; algorithm=direct changes nothing.)

The gridding artifact was REFMAC's (b_add_loc from the tiny B-factors) and it HURT refmac:
default refmac R-free 0.122 -> 0.116 with SHANNON 4-6 or Blim 2. So part of the *apparent*
phenix-over-refmac advantage was refmac's grid artifact (~0.006), not phenix doing
anything special.

Clean single-model decomposition (rerefine_003; all loose ends now closed):
| refmac setup (NCYC 0, k_mask Fpart) | R-free | adds |
|---|---|---|
| default grid                        | 0.117 | baseline |
| + SHANNON 6                         | 0.114 | grid artifact fixed (~0.003 here) |
| + CIF instead of PDB                | 0.114 | serialization: NO effect on R |
| derived Fpart (FMODEL-FCALC_R)+SHAN6| 0.107 | per-bin k_iso envelope (~0.007) |
| phenix                              | 0.106 | -- |

PDB-vs-CIF serialization is a NON-ISSUE for R: refmac gives identical R (0.1140 PDB vs
0.1139 CIF) and its output FC agrees to 0.4% between the two, because (a) refmac's own
reader/reprocessing lands ~PDB-precision either way (its Fcalc(0,6,0)~287 from both files,
not the CIF's 301.5 that cctbx reads), and (b) the per-reflection Fcalc differences are
bidirectional and wash out in Sum|Fo-Fc|/Sum|Fo|. Serialization affects individual
reflections and maps, not R. So the whole gap = grid artifact (~0.003-0.006, fixed by
SHANNON/Blim) + genuine per-bin scaling (~0.007, closed only by Fpart=FMODEL-FCALC_R).
None of it is a phenix inflation; phenix's number is real.

Free-flag reminder (NOT different free sets): phenix R_FREE_FLAGS is binary (free=1),
refmac FreeR_flag is 20-bin CCP4 (free=value 0). They flag the SAME reflections (100%
agreement) -- opposite encoding only. Select ==1 for phenix, ==0 for refmac when
reconstructing R.

Planned test (user): re-refine phenix with low-occupancy atoms dropped, to check whether
the diffuse low-occ conformers are real signal or overfitting (watch R-free vs R-work),
and to make the model faithfully PDB-representable (those are the atoms the PDB rounds to
occ 0.00).

## Known limitations

- Joint X-ray/neutron refinement path in `validate()` is not wired up
- Twin refinement untested (user f_mask_twin is not set)
- ADP refinement harmful for high-copy ensembles (see above)
