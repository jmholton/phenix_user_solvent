---
name: feedback-outputs-in-project-not-scratch
description: "Put refinement results/outputs in a subdir of the project, not the session scratch dir"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 2b5521f9-1ab5-46ea-8597-a7b2c393f71f
  modified: 2026-08-17T23:01:39.880Z
---

When a run produces results worth keeping (refinement outputs: pdb/cif/mtz/log/f_model, etc.),
write them into a subdir under the current project (e.g. `/home/jamesh/projects/phenix_user_solvent/<rundir>/`),
NOT the session scratch directory. The user asked for this explicitly and with some frustration
("why do you always use scratch dirs? just make a subdir under THIS project").

**Why:** the scratch dir is session-temporary and hard for the user to find/keep; real results
need to persist in-project alongside the rest of the work.

**How to apply:** default run output dirs to a named subdir of the project. Scratch is fine only
for genuinely throwaway temp files (intermediate scripts, parsing scratch). If a run must go to
scratch first (e.g. to avoid clobbering), copy the results into a project subdir afterward.
Still honor "don't overwrite" — use a new subdir/prefix rather than writing over existing outputs.
See [[User-supplied bulk solvent map feature]].
