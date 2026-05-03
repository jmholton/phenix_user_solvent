# Example: user-supplied bulk solvent map

Canonical comparison of user-supplied vs default bulk solvent for a
multi-conformer structure of 1AHO refined with `phenix.refine`.

## Files

| File | Description |
|---|---|
| `under20.pdb` | Multi-conformer model (conformers with occ ≥ 0.20 retained) |
| `refme.mtz` | Experimental data (Fobs, FreeR_flag) |
| `solvent_Fpart.mtz` | Ensemble bulk solvent (Fpart amplitudes, PHIpart phases) |
| `under20_occ_groups.eff` | Per-residue constrained occupancy groups for phenix |
| `under20_usersolvent_xyzoB_001.log` | Reference log — user solvent run |
| `under20_defaultsolvent_xyzoB_001.log` | Reference log — default flat-mask run |
| `fill_missing_atoms.py` | Fill missing atoms in chain-as-conformer PDB format |
| `make_occ_groups.py` | Generate per-residue constrained occupancy group .eff file |
| `reoccupy.awk` | Convert chain-as-conformer to altloc-as-conformer PDB format |

`solvent_Fpart.mtz` was computed from an AMBER MD trajectory of the same
crystal, followed by density editing in which strong difference features from
Fo-Fc maps (computed without free-flagged reflections) were added back.
The amplitudes and phases replace the internally computed flat
0/1 mask in phenix.refine's bulk solvent model:

```
F_model = k_total * (F_calc + k_mask * F_mask)
```

`k_mask` continues to be fit analytically per resolution bin.

## Commands

User-supplied bulk solvent (10 macro cycles, xyz + ADP + occ):

```
/programs/phenix-2.1rc2-6037/bin/phenix.refine \
  under20.pdb refme.mtz \
  refinement.input.bulk_solvent_map.file_name=solvent_Fpart.mtz \
  refinement.input.bulk_solvent_map.amplitudes_label=Fpart \
  refinement.input.bulk_solvent_map.phases_label=PHIpart \
  under20_occ_groups.eff \
  refinement.main.number_of_macro_cycles=10 \
  "refinement.refine.strategy=individual_sites individual_adp occupancies" \
  output.prefix=under20_usersolvent_xyzoB --overwrite
```

Default flat-mask bulk solvent:

```
/programs/phenix-2.1rc2-6037/bin/phenix.refine \
  under20.pdb refme.mtz \
  under20_occ_groups.eff \
  refinement.main.number_of_macro_cycles=10 \
  "refinement.refine.strategy=individual_sites individual_adp occupancies" \
  output.prefix=under20_defaultsolvent_xyzoB --overwrite
```

## Results

| Bulk solvent | R-work | R-free |
|---|---|---|
| Default (flat mask) | 0.129 | 0.144 |
| Ensemble Fpart/PHIpart | **0.086** | **0.104** |
