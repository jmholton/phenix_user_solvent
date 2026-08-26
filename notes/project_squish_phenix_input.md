---
name: squish-phenix-input
description: Augmented squish_solvent_runme.com to accept phenix _f_model.mtz input; iterative solvent-map updater
metadata: 
  node_type: memory
  type: project
  originSessionId: 2b5521f9-1ab5-46ea-8597-a7b2c393f71f
  modified: 2026-08-26T15:21:04.164Z
---

squish_solvent_runme.com (`/home/jamesh/projects/squish_solvent/`) is James's "squeeze"-style
bulk-solvent updater: Fo-Fc map -> mask out protein per conformer -> keep significant solvent-region
features -> write `refme_Fpart.mtz` (Fpart/PHIpart = missing-solvent INCREMENT) + `refme_minusol.mtz`.
Built for refmac output. Deps (in /home/jamesh/Develop, on PATH): map_scaleB.com map_func.com
Tc_maskify.com addup_maps_runme.com pick.com. Iterative wrapper: squishNfit_runme.com (loops
squish -> refmac refine scpart -> rescale). See [[User-supplied bulk solvent map feature]].

**Augmentation (2026-08-20): squish now accepts phenix _f_model.mtz directly.**
- `/home/jamesh/projects/phenix_user_solvent/phenix2squish.py` (iotbx, on PATH via Develop symlink):
  converts a phenix _f_model.mtz -> squish-ready mtz. Flips free flag (phenix binary free=1 ->
  CCP4 free=0), keeps phenix column names, adds Fpart=K_MASK*FMASK (current physical solvent).
  Free-flip verified exact (1601 free <-> 1601).
- squish_solvent_runme.com edited (backup `.prePhenix`): auto-detects phenix input (FMODEL present,
  no FC_ALL_LS) -> runs phenix2squish.py -> adopts phenix columns (FP=FOBS, SIGFP=SIGFOBS,
  FC_subtract=FMODEL,PHIFMODEL, FC_addback=FMODEL,PHIFMODEL, FreeR_flag flipped) + guards the
  refmac-only diagnostic FFT blocks (DELFWT/FWT/FC_ALL - dead code anyway). tcsh -n clean.
- Run high-res phenix data (0.96 A) with oversample=5, NOT the default 10 (10 -> 0.45 A grid).
- Validated end-to-end on best_so_far (rerefine_rerun_003) -> produced refme_Fpart.mtz.
  Work dir: `/home/jamesh/projects/phenix_user_solvent/squish_phenix/`.

**float_func stale-binary gotcha (cost real debugging time):** map_func.com self-compiles float_func
from embedded C, but a STALE `~/bin/float_func` (2023) shadows it on PATH and SEGFAULTS on `-func
prefft` (the "pre-fft conditioning of sharp edges" step). The bespoke_amber_restraints/float_func
(2026-07) also segfaults. FIX: recompiled fresh from map_func.com's embedded source (works) and
deployed to `/home/jamesh/Develop/float_func`; runs must have Develop ahead of ~/bin on PATH.
~/bin/float_func should be replaced (not done - it's the user's binary).

**ELECTRON UNITS (verified 2026-08-20):** gemmi map2sf preserves electron units EXACTLY -- CCP4
fft(FMODEL) -> gemmi map2sf round-trip gives sum/median ratio 1.0000 vs FMODEL. (An earlier claim
that gemmi used an "arbitrary scale" was WRONG.) Everything in this pipeline is/should be in electron
units. James's principle: keep everything in electron units so combining is plain addition.
squish's refme_Fpart Fpart looked inflated (max ~2774 vs FOBS ~567) ONLY because of squish's own
fofc_prescale=2 and Wilsonification (per-shell reshaping) -- both are DETECTION aids, not gemmi.

**THE OPEN PROBLEM - one-Fpart hand-back:** squish's Fpart is the INCREMENT; refmac's loop scales it
independently (scpart 2, separate Fpart2) then merges. Patched phenix supports only ONE user Fpart,
so old solvent + increment must be pre-combined into ONE Fpart/PHIpart. Because everything is electron
units, the combine is plain addition: new_solvent = K_MASK*FMASK + increment (electron units), handed
to phenix which re-fits one per-bin k_mask (~1). REQUIREMENT: the increment must be in TRUE electron
units -- generate it with squish fofc_prescale=1 wilsonify=0 (prescale/wilson distort the scale), OR
de-scale the output. NOT yet built; confirm detection-vs-electron-scale tradeoff with James.

**FIRST squish->phenix ITERATION WORKED (2026-08-20).** Full loop: augmented squish accepts phenix
f_model -> guarded Wilsonification -> combine_solvent.py -> phenix.refine. Result vs best_so_far
(same rerefine_002.pdb start, same rerefine_003.eff, seed 2679941, only the solvent map differs):
| solvent | Start Rfree | Final Rwork | Final Rfree |
| old (example/solvent_Fpart) = best_so_far | 0.1082 | 0.0874 | 0.1065 |
| squish-UPDATED (this loop)                | 0.1056 | 0.0819 | 0.1021 |
dRfree = -0.0044, dRwork = -0.0055. Start Rfree already better (0.1082->0.1056) => the updated
solvent MAP itself fits better at cycle 0. This BREAKS past the refmac ceiling (~0.110) that the
1280-cycle refmac longrun/squish25 plateaued at. Output: squish_phenix/refine_updated/updatedsolvent_003*.
NEW BEST -- not yet promoted to best_so_far (ask James).
Confirmed by refmac NCYC0 free LLG (the sensitive metric): scoring both phenix models the same way,
-LLfree 4853.9 (old) -> 4815.0 (updated) = dLLGfree -38.9; refmac Rwork 0.0913->0.0863,
Rfree 0.1124->0.1095, FreeFSC 0.9858->0.9861. All metrics agree. That single squish->phenix cycle
moved free LLG MORE than the whole 1280-cycle refmac longrun (-22 total). Scorer:
squish_phenix/llg_score/ (refmac NCYC0, SHANNON 6, Fpart=K_MASK*FMASK via phenix2squish.py).

Wilsonification GUARD (in repo squish_solvent_runme.com): the low-res line is fit with a SQUARED slope
(loBr, loB=loBr^2) so loB>=0, preventing the negative-loB over-sharpening (was -27.98 on this 0.96A
data -> now 0.77). No-op when loB is naturally positive (e.g. squish25's loB=+25). Root trigger:
anomalously low origin Wilson bin tilts the narrow lofrac=0.1 fit to a rising line.

combine_solvent.py (repo root): new_FMASK = K_MASK*FMASK + k*increment/(K_ISO*K_ANISO), electron units.
k = single divot scale fit to fofc=FOBS-FMODEL over the WORK set (phenix analog of refmac scpart;
got k=0.83). /(K_ISO*K_ANISO) moves the increment from data scale to FMASK's pre-overall scale so
phenix's fresh k_mask starts ~1. Generate the increment with fofc_prescale=1 + the loB guard.

FREE-SET EXCLUSION VERIFIED (result is NOT overfit): free reflections are excluded from squish's
fofc difference map (exclude_freeR=1 -> fft FREE=FreeR_flag; phenix2squish flips R_FREE_FLAGS free=1
-> FreeR_flag free=0, 1601<->1601) AND from combine_solvent.py's k-fit (work-only). GOTCHA: CCP4 fft
prints "No reflections excluded" even when FREE= IS working -- that message is about a different check;
proven by fft-with vs fft-without giving different maps (rms 0.0418 vs 0.0403). So Rfree/free-LLG gains
are genuine cross-validation.

CONVERGED IN ONE CYCLE (2026-08-20). Free-LLG trajectory (refmac NCYC0 scorer, -LLfree):
cycle0/old 4853.9 -> cycle1 4815.0 (-38.9) -> cycle2 4815.7 (+0.7 = flat). phenix Rfree
0.1065 -> 0.1021 -> 0.1024. cycle2 increment was tiny (solvent map changed ~1%, combine k=1.96
but net +1% amplitude), no further gain. best_so_far stays at CYCLE 1 (updatedsolvent_003,
Rfree 0.1021, -LLfree 4815.0); cycle2 outputs in squish_phenix/refine_cycle2/ (not promoted --
edges out by cycle1 on free LLG + phenix Rfree). Promotion done: best_so_far now = cycle1,
old best in best_so_far/prev_rerefine_rerun_003/.

FLOODFILL vs SIGNIFICANCE-MASK A/B (2026-08-20, same cycle-1 input, only mask mode differs):
squish25's floodfill mode (run_long_next.com): exclude_protein=0 floodfill=1 significance=none
floodfill_sigma=1 floodfill_Bsmooth=1 floodfill_maskBsmooth=20 floodfill_rethresh=0.9
floodfill_maxradius=5 floodfill_recycles=3 floodfill_fftmargin=5 (was fftbuffer in old script)
floodfill_highprob=0.6 floodfill_lowprob=0.5 floodfill_retries=10. (squish25 itself mostly used the
SIGNIFICANCE mask exclude_protein=0; floodfill was in run_long_next.com.) Results (refmac NCYC0 -LLfree
/ phenix Rfree): cycle0/old 4853.9/0.1065; significance mask 1cyc 4815.0/0.1021 (-38.9); floodfill 1cyc
4847.0/0.1049 (-6.9). Floodfill = ONE segmented peak per cycle (this divot: solvent mean +0.6%), so one
cycle ~1/6 of the all-features significance-mask gain -- it's built for the ~620-cycle iterated loop.
Per-cycle the significance mask wins; best_so_far stays significance-mask cycle 1 (0.1021). Iterated
floodfill MIGHT reach a cleaner endpoint (squish25's refmac loop hit -LLfree ~4744) but needs many
cycles and wasn't run to convergence here. Outputs: squish_phenix/floodfill/, refine_floodfill/.

FUNDAMENTAL LIMIT CHECK vs SIGF (2026-08-21, James's framing: data should be explainable to within
experimental error; Rfree is nowhere near that). On best_so_far updatedsolvent_003 f_model:
- Rfree 0.1020 vs experimental-error floor Sum(sig)/Sum(F) = 0.0342 -> Rfree is 3.0x the floor
  (a perfect model ~0.8x floor ~0.027). NOT overfitting -- large genuine unexplained signal remains.
- <|Fo-Fmodel|/sigF> free = 4.03 mean, 2.43 median (vs ~0.8 for a model at experimental error).
- Per-shell (free): the misfit is at LOW/MID resolution -- 1.9-2.7A has <|dF|/sig> = 6-7 (data most
  precise there, R_floor ~0.012, yet model 6sig off). High-res R (0.24 at 0.96A) is mostly
  experimental error (weak data, R_floor 0.15, misfit only 1.6sig) -- little room there.
- Implication: the frontier is unmodeled STRUCTURE at 1.4-2.7A (disorder/ADP/anisotropy/ordered
  solvent), not the last milli-R of bulk solvent. squish->phenix (0.1065->0.1021) was a real but
  small step; don't get greedy chasing 0.101 -- the factor-of-3 gap to experimental error is real
  signal concentrated at medium resolution.
Tool: recompute anytime from a phenix _f_model.mtz (FOBS/SIGFOBS/FMODEL/R_FREE_FLAGS).

WILSONIFICATION vs SNR (2026-08-21): James asked if Wilsonification is the way to blur the diff
density and capture the right resolution range. Measured the guarded (loB=0.77) applied scale per
shell (=exp(0.5*(logF2_after-logF2_before))): suppresses low/mid res (0.80 @9.6A, 0.56 @2.9A) and
BOOSTS high res (up to 3.83 @1.26A, 2.4 @0.99A) -- a HIGH-PASS/sharpen, NOT the intended blur. This
is ANTI-correlated with the SNR profile (|dF|/sig): it downweights the 2-2.9A band where SNR peaks
(~6-7sig) and boosts the ~1A end where data is ~noise (1.6sig). Root: straightening normalizes by
the diff map's OWN power (signal+noise), so it amplifies high-res noise. Explains why Wilsonification
still beat RAW fofc (raw is low-res-dominated; Wils de-weights low res but overshoots into high-res
noise; optimum is in between at mid-res). PROPER tool = sigma-weighted Wiener/matched filter:
w(s) proportional to SNR(s) from experimental sigF -- emphasize 1.9-2.7A (6-7sig), roll off ~1A shells
(genuine blur). We have the per-shell |dF|/sig to build it. Proposed: prototype an SNR-weighted diff
map for squish to replace Wilson-plot straightening; targets the 6sig medium-res frontier directly.

WILSONIFICATION FIX + RESULT (2026-08-21): James saw the Wilson output going flat (not preserving the
low-res slope). Root cause: (a) the origin bin (stol2=0, logF2~2.3, a difference-map DC outlier ~4 log
units below neighbors) was in the low-res line fit, flipping the slope (incl origin -> loB=-31.8;
excl -> +6.6); my earlier loB>=0 guard then clamped it flat (~0.77). (b) lofrac=0.1 fit only the lowest
9 noisy bins (d>3A) -> unstable shallow loB. FIXES in repo squish_solvent_runme.com: drop origin bin
(if(bin+0>0) in the Wilson binning END, both occurrences) AND wilsonify_lofrac 0.1 -> 0.25. Result:
loB=19.4 (matches James's eye ~20.3), output Wilson plot drops 10.7 log units (follows B~20), applied
scale is now a genuine BLUR (0.93 low-res -> 0.04 at 0.96A) instead of boosting high-res.
BUT pushed through combine+refine it did NOT beat the flat version: corrected-blur -> phenix Rfree 0.1025
(vs flat 0.1021), refmac Rfree 0.1083 (BEST), free LLG 4835.2 (WORSE than flat 4815.0). Mixed; on the
sensitive metric (free LLG) the monotonic B~20 blur captured only ~half the signal. WHY: the SNR peaks
at 1.9-2.7A; a monotonic blur SUPPRESSES that mid-res band (scale 0.7-0.8) while the flat straightening
incidentally BOOSTED it (1.6-2.9). CONCLUSION: the optimum is a BAND-PASS (sigma-weighted Wiener filter
peaking at the 1.9-2.7A SNR max), not a monotonic blur nor flat straightening. best_so_far stays flat-
Wils significance cycle1 (0.1021/4815). Wilsonification fixes are correct/committed-worthy but the
next real lever is the band-pass extraction. Outputs: squish_phenix/wilsonfix2/, refine_wf2/.

RESOLUTION-FILTER FRAMEWORK + BAND-PASS RESULT (2026-08-21). Added a pluggable resfilter= to squish
(all filters write the same per-refl scale ${t}scales.txt, then existing masking): resfilter=
wilson(default)|bandpass|blur built; wiener|snr PLANNED (need per-shell Pnn=<sigF^2> from $mtzfile
SIGFOBS -> H=Pss/(Pss+Pnn); phenix2squish already carries SIGFOBS into the squish mtz). bandpass params
bandpass_center/bandpass_width (stol^2), blur_B. Band-pass profile verified textbook (0.73 low -> 1.0 @
2.27A -> 0.002 @ 0.97A).
RESULT (refmac NCYC0 -LLfree | phenix Rfree): old 4853.9/0.1065; flat-Wils significance 4815.0/0.1021;
band-pass 4830.0/0.1029; monotonic blur 4835.2/0.1025. ALL resolution-filter variants that suppress
high-res LOSE to the flat straightening (which boosts high-res). Reshaping the solvent-extraction
resolution filter is NOT the lever. This CONFIRMS the SIGF diagnosis: the mid-res 6sig misfit is
unmodeled STRUCTURE (disorder/ADP/anisotropy), NOT bulk solvent -- so no solvent filter (band-pass
correctly targeted 2A and still failed) can capture it. Solvent loop is PLATEAUED at 0.1021; the real
frontier is the protein model at medium resolution. Caveat: significance mask co-varies with the
filter (not a perfectly isolated resolution-weight test), but the ordering (all <= flat) is consistent.
Framework is a keeper capability (resfilter=wilson|bandpass|blur; wiener/snr next) even though it
didn't beat flat. Outputs: squish_phenix/bandpass/, refine_bandpass/.

WIENER/Pnn FILTER BUILT + TESTED (2026-08-22) -- completes the resfilter framework.
sigma_power.py (repo): per-shell Pnn=<sigF^2> binned by stol^2 (free excluded, origin bin dropped),
matched to squish's Wilson bins. squish resfilter=wiener|snr branch: H=1-Pnn/Ptot (wiener) or
Pss/Pnn normalized (snr), Ptot=exp(wilson log<F^2>) per map, Pnn from sigma_power.py. Wiener H came out
a self-tuned band-pass: ~0.97 across 6.8-1.9A (signal band), rolling to 0.16 at 0.97A, no hand-tuning.
FINAL LADDER (refmac NCYC0 -LLfree | phenix Rfree): flat-Wils 4815.0/0.1021 (BEST); WIENER 4825.0/0.1027;
band-pass 4830.0/0.1029; blur 4835.2/0.1025; old 4853.9/0.1065. The Wiener is the BEST filter (closest
to flat, as expected for the optimal filter) but STILL loses to flat. DEFINITIVE: no resolution filter
on the solvent extraction beats flat straightening. Likely flat's edge is a SPATIAL effect (sharper map
-> better-localized significance peaks; mask co-varies with filter), not resolution weighting. CONCLUSION
(rigorously confirmed): the recoverable signal is NOT extractable bulk solvent -- it's unmodeled protein
structure at mid-res (the SIGF 6sig). Solvent loop PLATEAUED at 0.1021 (best_so_far). Next frontier =
protein model (disorder/ADP/anisotropy at ~2A). resfilter framework (wilson|bandpass|blur|wiener|snr) +
sigma_power.py are keeper capabilities; commit-worthy. Outputs: squish_phenix/wiener/, refine_wiener/.

PROTEIN MODEL STATE (2026-08-22, best_so_far updatedsolvent_003, the mid-res frontier target):
1aho scorpion toxin, 64 residues, 10016 atoms, ISOTROPIC ADPs, 0 waters (solvent = bulk map only).
Disorder modeled ENTIRELY by a massive explicit ensemble: ~45 altloc conformers (16 near-complete
copies A-P + minor tail), total occ ~0.90 copies. ADPs are IDLE: median B 6.26, mean 6.70, 83% of atoms
B4-8, meanB ~6.5 in every occ band -- B carries almost no variation, all displacement signal is in the
ensemble geometry. INCOMPLETE: mean per-(atom,residue) total occupancy = 0.903 (CIF-confirmed, not a PDB
rounding artifact), 124 atoms <0.5 occ (real gaps), 34 >1.05 (over-modeled). So ~10% of ordered
scattering is unmodeled. Two structural candidates for the 2A/6sig misfit (NOT solvent): (a) the ~10%
missing occupancy (124 thin atoms -> unmodeled density -> Fo-Fc); (b) pinned-low ADPs can't absorb the
continuous-disorder residual the discrete 45-copy ensemble misses (individual_adp is unstable for high-
copy ensembles per the 1/occ gradient blowup, so B was left low/idle -- may now COST the mid-res fit).
Next: locate WHERE the 124 thin atoms / 6sig free reflns concentrate in real space (which residues).
Reminder: PDB rounds low occ to 0.00 (2252 atoms) -- use the CIF for real occupancies.

Fo-Fc PEAK LOCALIZATION (2026-08-22, pick.com symmetry-aware). Biggest Fo-Fc features of best_so_far
are MODEST (<=4.8sig) and PROTEIN-REGION: 51/52 positive + 26/26 negative peaks within 3.2A of a heavy
atom (most 0.6-2.7A, on/beside atoms), only 1 positive in bulk solvent. Cluster on flexible/charged
surface side chains: +near PRO41,ARG18,HIS54,GLY43,LYS30,ASP9,VAL10; -near ARG62,LYS50,GLU32,PRO52,
ASP8. Across filters the top +peak is 4.7sig in ALL (flat/bandpass/wiener) -- invariant to solvent
filter; flat gives the FEWEST peaks (52/26 vs bandpass 62/35, wiener 55/31) = cleanest residual = best
fit. DEFINITIVE spatial confirmation: residual signal is unmodeled protein (positional/ADP/H fine-
structure at flexible side chains, ~2A), NOT solvent. Next lever = protein model there, not solvent.
pick.com usage (James): pick.com <map> <sigma_bare_number> <model.pdb> <dist>A  (e.g. 0.0001A sets
CLOSE_peaks so ALL peaks kept + annotated with sym-aware dist to nearest atom). GOTCHAS: bare number w/o
A = sigma, number WITH 'A' suffix = CLOSE_peaks; site pdb MUST have CRYST1 (coordconv needs the cell) or
0 sites load; ~5000 sites x symmetry is slow (run backgrounded). Outputs: squish_phenix/fofc_peaks/.

EXCLUDE_PROTEIN=0 GAMUT (2026-08-22, James's test: my solvent cycles all used exclude_protein=1 =
protein masked out, so they couldn't touch the protein-region residual; ep0 lets squish density-modify
the protein region too -- paint observed Fo-Fc protein density into the Fpart). Ran wilson/bandpass/
blur/wiener with exclude_protein=0, same base (prev-best f_model)/eff/seed, free excluded from fofc.
RESULT (phenix Rfree | refmac Rfree | -LLfree): wilson 0.1028|0.1078|4824.0; bandpass 0.1033|0.1093|
4835.3; blur 0.1053|0.1111|4850.3; wiener 0.1023|0.1085|4818.9. NONE beat the ep1 flat best (0.1021 |
0.1095 | 4815.0). Closest: wiener-ep0 (0.1023/4818.9), marginally worse. TELL: wilson-ep0 gave the best
refmac Rfree yet (0.1078) with large k (1.4-2.4 = lots of protein density added) and better R-WORK, but
free LLG did NOT improve -> OVERFITTING the work set; the protein-region residual is distributed fine-
structure (<=4.7sig over many side chains), not coherent missing density, so density-modification (Fpart
paint-in) can't recover it -- needs a real MODEL change (geometry/ADP/conformers). PLATEAU at 0.1021 now
confirmed 3 ways: protein masked (ep1), protein density-modified (ep0), and free-excluded throughout.
best_so_far stays 0.1021. Outputs: squish_phenix/gamut_ep0/ (run_gamut.sh driver + summary.txt).

DIVOT (floodfill) ep0 GAMUT (2026-08-23): James distinguished my filter tests (all WHOLE-MAP
significance masks, floodfill=0) from DIVOTS (floodfill=1, one segmented feature). Ran divot+ep0 gamut
(squish25 floodfill params, free excluded): wilson 0.0868/0.1050/-LLfree4849.5; bandpass 0.0868/0.1054/
4850.2; wiener 0.0867/0.1048/4845.1; blur FLOODFILL FAILED -- actual error "flood-fill hit edge of box"
(NOT "no peak"; corrected). The B~20 blur makes the diff density broad/diffuse (no sharp foothills), so
floodfill's flooded region never closes into a bounded blob -- it spreads to the 5A box wall; retry
escalated floodfill_sigma 1.0->1.9 but still 22 voxels on the box edge -> fail. CORRECTED (James): NOT "incompatible" -- just box-vs-shape. Re-ran blur+floodfill with
floodfill_maxradius=10 (was 5): completes fine, phenix Rfree 0.1052 / -LLfree 4845.3, right in the
divot-ep0 pack. So a diffuse (blurred) feature just needs a box big enough to contain it. Completed
divot-ep0 gamut: wilson 0.1050/4849.5, bandpass 0.1054/4850.2, blur(box10) 0.1052/4845.3, wiener
0.1048/4845.1 -- all cluster ~0.105, all worse than flat ep1 (0.1021); per-cycle divot is too
conservative (solvent delta ~+0.3%), filter choice is second-order. Output: squish_phenix/blur_bigbox/.. All WORSE than whole-map ep0 (4818-4824) and flat ep1 (0.1021/4815). Divot
per-cycle is too conservative: one small feature (map rms 0.0018, k~1.2-1.6) ~ no-op (Rfree 0.1065->
0.1050 like the earlier 1-cycle floodfill). Overfitting check DID come out clean for divots (work-free
gap 0.0182 < flat's 0.0202 -- divots don't overfit like whole-map ep0 did) but they capture too little
per cycle; the divot approach only pays off ITERATED (squish25 ran ~620 cycles). NET: tested every combo
-- whole-map/divot x ep0/ep1 x wilson/bandpass/blur/wiener -- NOTHING beats flat solvent-only 0.1021/
4815. Plateau confirmed from every angle; frontier is a genuine MODEL change at the flexible side chains.
Outputs: squish_phenix/gamut_divot_ep0/.

SCALE-STABILIZATION PROTOCOL (James, 2026-08-23/24) -- a gap in my phenix workflow. squishNfit_runme.com
(lines 183-230) is James's protocol: add a NEW divot as a SEPARATE partial (Fpart2, scpart 2), then
iterate converge_refmac.com + refmac_rescale_runme.com until the divot's scale STABILIZES
(scaling_converged: |scale-1|<0.001 & |B|<0.01), OR it collapses to zero (zeroFparts) and is REJECTED as
not-real. ONLY THEN is the stabilized divot merged into the main Fpart (next cycle: Fpart += Fpart2).
"The scales must stabilize before the divot can be added to the main map." Use converge_refmac.com (NOT
converge_refmac_autopart.com -- less mature). refmac_rescale_runme.com is in ~/projects/squish_solvent
(not on PATH; add it, run don't edit). refmac_opts.txt config (longrun): solvent no; weight matrix 0.3;
scpart; damp 0.001 0.0002 0.0065; NCYC 50; make hout/hydr Y; occupancy refine + ~48 occ groups (one per
altloc) + potential_overrides.txt -- the heavy damping is for the 48-altloc Hessian instability.
MY GAP: combine_solvent.py fits one LS k and merges the divot IMMEDIATELY -- no stabilization, no zero-
collapse rejection. ROOT CAUSE = the single-Fpart limit: patched phenix allows only ONE user Fpart, so I
can't keep the divot as a separate Fpart2 to stabilize independently. CORRECT ADAPTATION: stabilize the
divot's scale in REFMAC (multi-partial, this damped config) -> get stable k* (or reject if it collapses)
-> merge THAT into phenix's single Fpart. Also answers if the ep0 protein divots are real (stabilize) or
noise (collapse).

refmac5-newhess (Hessian-fixed) BROKEN restore (2026-08-26): /data/home/programs/ccp4-8.0/bin/refmac5-
newhess and /programs/ccp4-8.0/bin/refmac5-newhess are byte-identical and BOTH SIGSEGV at startup
(refmac_xml.f MAIN) with every libmmdb2 tried (ccp4-8.0 and ccp4-9). No rpath (needs explicit
LD_LIBRARY_PATH=/programs/ccp4-8.0/lib, which ccp4.setup-sh does NOT set) AND ABI-mismatched with the
available libmmdb2 (its matching build-libs weren't restored; /data/home/programs/ccp4-8.0/lib lacks
libmmdb2 entirely). Stock /programs/ccp4-8.0/bin/refmac5 (5.8.0425) works fine. So the Hessian fix is
UNUSABLE until its libs/rebuild are restored -- told James. Fallback: stock refmac5 (ccp4-9) + light
damping (James: "shouldn't need all that damping with this model").

METRIC: watch FREE LLG (refmac_Rplot.txt col6 = -LL free), NOT just Rfree -- refmac was damped hard for
the 48-altloc Hessian instability so Rfree is coarse/frozen (longrun 0.1099 vs squish25 0.1089 differ
by a whisker in Rfree but 4760 vs 4744 in free LLG). Score phenix models by refmac NCYC0 free LLG.

**DO NOT EDIT ~/Develop or ~/bin (James's instruction, 2026-08-20).** Those hold mature scripts; James
copies things in once established. Make own copies in THIS repo (phenix_user_solvent). Working copies:
`phenix_user_solvent/squish_solvent_runme.com` (augmented), `phenix_user_solvent/bin/{float_func,
addup_maps_runme.com}`, `phenix_user_solvent/phenix2squish.py`. For runs prepend
`$REPO:$REPO/bin:~/Develop` so my float_func beats stale ~/bin and map_func/Tc_maskify/pick come from
~/Develop (read-only use OK). If ~/Develop or ~/bin look stale, TELL James. ~/bin/float_func (2023) is
stale (a restoration artifact, per James); ~/Develop/map_func.com is current (self-compiles a working
float_func from embedded source).
