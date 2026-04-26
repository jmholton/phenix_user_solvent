Example: user-supplied bulk solvent map (multi-conf structure)
=============================================================

Files:
  evenmoreconf.pdb  - multi-conformer model
  refme.mtz         - experimental data (Fobs, FreeR_flag)
  solvent_Fpart.mtz - ensemble-derived bulk solvent (Fpart, PHIpart)

The solvent model in solvent_Fpart.mtz was computed by Refmac from an
ensemble model of the same crystal. It contains Fpart (amplitudes) and
PHIpart (phases in degrees), equivalent to Refmac's FPART/PHIPART.

These structure factors replace the internally computed flat 0/1 mask in
phenix.refine's bulk solvent model:

  F_model = k_total * (F_calc + k_mask * F_mask)

k_mask continues to be determined analytically per resolution bin.

Command:

  phenix.refine \
    evenmoreconf.pdb \
    refme.mtz \
    refinement.input.bulk_solvent_map.file_name=solvent_Fpart.mtz \
    refinement.input.bulk_solvent_map.amplitudes_label=Fpart \
    refinement.input.bulk_solvent_map.phases_label=PHIpart \
    refinement.main.number_of_macro_cycles=1 \
    "refinement.refine.strategy=individual_sites individual_adp occupancies" \
    output.prefix=test_usermap --overwrite

Result (1 macro cycle, no RSR):
  Default bulk solvent:   R-work = 0.177  R-free = 0.180
  Ensemble bulk solvent:  R-work = 0.108  R-free = 0.119
