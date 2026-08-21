# Restoration instructions — phenix_user_solvent + 1aho refinement data

Written 2026-08-05 after a stale backup was restored over live work. This file tells the
restoration agent exactly what was lost, what is safe, and how to recover each piece.
**Read all of it before restoring anything — do NOT re-apply the same stale backup.**

## What happened
A system-level restore (files now owned `root:root`, dated `Aug 5 11:09`) overwrote the
working directories with a backup that PREDATES the current work. It clobbered:
- this git repo (`/home/jamesh/projects/phenix_user_solvent`) — `.git` is now corrupt
  (has `objects/`+`refs/` but no `HEAD` or `config`); `CLAUDE.md` and several scripts gone.
- the crystallographic data tree `/home/jamesh/projects/amber.old/1aho_refine/best_Rfree/`
  — restored to a PRE-`rerefine_003` state (only old files like `scales.mtz`,
  `kmask_vs_stol.txt` remain; the actual model/data are gone).

The live work spanned **2026-06-21 to 2026-07-21**. Any backup used to restore the data
MUST be from **2026-07-10 or later** (ideally 2026-07-21+). The stale one predates 06-21.

## SAFE — already recoverable, no action from a backup needed
1. **Repo code + documentation** → GitHub `github.com/jmholton/phenix_user_solvent`
   - branch **`refmac-phenix-transfer`** @ commit **`f43742e`** = the current work
     (`CLAUDE.md` with the full refmac<->phenix study, plus `wilsonify.py`,
     `rescale_fobs.py`, `extrapolate_kmask.py`, `package_for_refmac.py`,
     `wilson_b_fmask.py`, `slide.html`).
   - branch `master` @ `9975c8c` = older, pre-session.
   Recover by re-initialising git from the remote and checking out `refmac-phenix-transfer`
   (do NOT overwrite with the stale working tree). The corrupt local `.git` should be
   replaced from the remote.
2. **Deconform / realdata line** (survived on disk, do NOT touch/overwrite):
   `/scratch/jamesh/projects/amber/1aho_refine/deconform_under20/` — contains
   `realdata.mtz`, `realdata_0003_phenix/weightscan_005_001.pdb`, `phenix_001..004.{pdb,mtz}`,
   `best_for_000*.pdb`.

## LOST — must be restored from a NON-STALE backup (NOT in git, no copy on disk)
Confirmed absent everywhere on `/home` and `/scratch`. Restore these from a backup dated
>= 2026-07-10:

`/home/jamesh/projects/amber.old/1aho_refine/best_Rfree/phenix_usersolvent/`
  - `rerefine_003.pdb`          (45-conformer ensemble model)
  - `rerefine_003.cif`          (SAME model, FULL precision — needed; the PDB is lossy)
  - `rerefine_003_f_model.mtz`  (CRITICAL: columns FOBS, R_FREE_FLAGS, FMODEL, FCALC,
                                 FMASK, K_ISOTROPIC, K_ANISOTROPIC, K_MASK — everything
                                 the Wilsonification/transfer tooling needs)
  - `rerefine_003.mtz`, `rerefine_003.eff`, `rerefine_003.log`
  - `rerefine_002.pdb`

`/home/jamesh/projects/amber.old/1aho_refine/best_Rfree/refmac_phenixsolvent/`
  - `refme.mtz`  (symlink -> `refme_with_fpart.mtz`; Fpart = K_MASK*FMASK, the
                  k_mask-scaled physical solvent for refmac)
  - `ncyc0/`     (`notzero.pdb`, `refmacout.mtz`, `converge0.log`, `converge_shan*.log`,
                  `converge_Blim*.log`)

Timestamps to sanity-check a candidate backup (these are the versions we need):
  `rerefine_003.*` = 2026-06-21 ~08:02;  `refme_with_fpart.mtz` = 2026-06-29;
  `ncyc0/*` (converge_shan6.log, notzero.pdb) = 2026-06-30 .. 2026-07-10.
A backup whose `best_Rfree/phenix_usersolvent/` contains only `scales.mtz` /
`kmask_vs_stol.txt` and NO `rerefine_003*` is the STALE one — reject it.

## If the data backup cannot be found
The `best_Rfree/rerefine_003` line is then unrecoverable and must be re-refined from
scratch. The deconform line (survived) can continue instead: any phenix `_f_model.mtz`
(none survived anywhere) is regenerable with
`phenix.refine <model> <data> output.export_final_f_model=True`
(use `/programs/phenix-2.1rc2-6037/bin/phenix.refine`), after which `wilsonify.py` works.

## Post-restore verification
- `ls best_Rfree/phenix_usersolvent/rerefine_003_f_model.mtz` exists; its columns include
  `K_ISOTROPIC K_ANISOTROPIC K_MASK FMASK` (check with `phenix.mtz.dump` or iotbx).
- `CLAUDE.md` present and contains the string "Wilsonif" (confirms non-stale repo).
- Fix ownership: restored dirs are `root:root`; chown back to `jamesh:831_staff`.
- Do not overwrite the surviving `/scratch/.../deconform_under20/` tree.
