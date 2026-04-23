Example: user-supplied bulk solvent map (1AHO)
==============================================

The solvent model in solvent_Fpart.mtz was computed by Refmac using a
48-copy ensemble model of the same crystal. It contains Fpart (amplitudes)
and PHIpart (phases in degrees), equivalent to Refmac's FPART/PHIPART.

These structure factors replace the internally computed flat 0/1 mask in
phenix.refine's bulk solvent model:

  F_bulk = k_sol * exp(-B_sol * s^2/4) * F_user

k_sol and B_sol continue to refine against the user-supplied map.

Command:

  phenix.refine \
    1aho.pdb \
    1aho.mtz \
    refinement.input.bulk_solvent_map.file_name=solvent_Fpart.mtz \
    refinement.input.bulk_solvent_map.amplitudes_label=Fpart \
    refinement.input.bulk_solvent_map.phases_label=PHIpart

Result (1 macro cycle, no RSR):
  Input:  R-work = 0.165  R-free = 0.160
  Output: R-work = 0.151  R-free = 0.147
