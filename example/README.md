# Example: user-supplied bulk solvent map

Canonical comparison of user-supplied vs default bulk solvent for a
multi-conformer structure of 1AHO refined with `phenix.refine`.

## Files

| File | Description |
|---|---|
| `starthere.pdb` | Starting multi-conformer model |
| `refme.mtz` | Experimental data (Fobs, FreeR_flag) |
| `solvent_Fpart.mtz` | Ensemble bulk solvent (Fpart amplitudes, PHIpart phases) |
| `usersolvent_001.*` | Full phenix job — user-supplied bulk solvent |
| `usersolvent_data.mtz` | Data MTZ written by phenix for usersolvent run |
| `defaultsolvent_001.*` | Full phenix job — default flat-mask bulk solvent |
| `fill_missing_atoms.py` | Fill missing atoms in chain-as-conformer PDB format |
| `make_occ_groups.py` | Generate per-residue constrained occupancy group .eff file |
| `reoccupy.awk` | Convert chain-as-conformer to altloc-as-conformer PDB format |

`solvent_Fpart.mtz` was computed from an AMBER MD trajectory of the same
crystal, followed by density editing in which strong difference features from
Fo-Fc maps (computed without free-flagged reflections) were added back.
The amplitudes and phases replace the internally computed flat 0/1 mask in
phenix.refine's bulk solvent model:

```
F_model = k_total * (F_calc + k_mask * F_mask)
```

`k_mask` continues to be fit analytically per resolution bin.

## Commands

User-supplied bulk solvent (3 macro cycles, xyz + ADP):

```
/programs/phenix-2.1rc2-6037/bin/phenix.refine \
  starthere.pdb refme.mtz \
  refinement.input.bulk_solvent_map.file_name=solvent_Fpart.mtz \
  refinement.input.bulk_solvent_map.amplitudes_label=Fpart \
  refinement.input.bulk_solvent_map.phases_label=PHIpart \
  refinement.main.number_of_macro_cycles=3 \
  "refinement.refine.strategy=individual_sites individual_adp" \
  output.prefix=usersolvent --overwrite
```

Default flat-mask bulk solvent:

```
/programs/phenix-2.1rc2-6037/bin/phenix.refine \
  starthere.pdb refme.mtz \
  refinement.main.number_of_macro_cycles=3 \
  "refinement.refine.strategy=individual_sites individual_adp" \
  output.prefix=defaultsolvent --overwrite
```

## Results

| Bulk solvent | R-work | R-free |
|---|---|---|
| Default (flat mask) | 0.128 | 0.147 |
| Ensemble Fpart/PHIpart | **0.093** | **0.109** |
