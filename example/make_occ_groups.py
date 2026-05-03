"""
Generate a phenix .eff file with explicit per-residue occupancy constrained
groups for a multi-conformer PDB.

Handles two formats:
  - chain-as-conformer (minconf1.pdb): each chain letter = one conformer,
    altloc = chain. Selection uses chain.
  - altloc-as-conformer (onechain.pdb): single chain, altloc field = conformer.
    Selection uses altloc + chain.

Usage: python make_occ_groups.py model.pdb > occ_groups.eff
"""
import sys
from collections import OrderedDict

WATER_RESNAMES = {"HOH", "WAT", "DOD", "H2O", "SOL"}

pdb_file = sys.argv[1]

# conformer key: (chain, altloc_stripped)
# residue key: (resseq_str, icode)
# map: residue_key -> list of conformer keys (in order of first appearance)
residues = OrderedDict()

chains_seen = set()
altlocs_seen = set()

with open(pdb_file) as f:
    for line in f:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resname = line[17:20].strip()
        if resname in WATER_RESNAMES:
            continue
        altloc = line[16].strip()
        chain  = line[21]
        resseq = line[22:26]
        icode  = line[26]
        chains_seen.add(chain)
        if altloc:
            altlocs_seen.add(altloc)
        key = (resseq, icode)
        conf = (chain, altloc)
        if key not in residues:
            residues[key] = []
        if conf not in residues[key]:
            residues[key].append(conf)

def make_selection(chain, altloc, resid):
    if altloc:
        return "chain %s and altloc %s and resid %s" % (chain, altloc, resid)
    else:
        return "chain %s and resid %s" % (chain, resid)

print("refinement {")
print("  refine {")
print("    occupancies {")

n_groups = 0
for (resseq, icode), confs in residues.items():
    if len(confs) < 2:
        continue
    resid = resseq.strip() + (icode.strip() or "")
    print("      constrained_group {")
    for (chain, altloc) in confs:
        print("        selection = %s" % make_selection(chain, altloc, resid))
    print("      }")
    n_groups += 1

print("      individual = water")
print("    }")
print("  }")
print("}")
print("# %d constrained groups written" % n_groups, file=sys.stderr)
