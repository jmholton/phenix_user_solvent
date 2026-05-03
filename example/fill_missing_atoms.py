"""
fill_missing_atoms.py — chain-as-conformer PDB format

For each protein residue, every atom that appears in *any* conformer (chain)
is added to *all* conformers that are missing it.  Coordinates are copied from
the first chain that has the atom; occupancy is taken from another atom already
present in the target chain (so the filled-in atom gets the correct per-chain occ).

Usage: python3 fill_missing_atoms.py model.pdb > model_filled.pdb
"""
import sys
from collections import defaultdict, OrderedDict

lines = open(sys.argv[1]).readlines()

# ── 1. Parse all ATOM/HETATM records ──────────────────────────────────────
atom_rec   = {}           # (chain, resseq, icode, resname, atomname) -> line (str, no newline)
res_order  = []           # residue keys in file order, unique
res_chains = OrderedDict()  # rkey -> list of chains in file order (unique)
res_atoms  = OrderedDict()  # rkey -> list of atomnames in file order (unique)
chain_occ  = {}           # (chain, resseq, icode, resname) -> occupancy string (6 chars)

def _add_unique(lst, val):
    if val not in lst:
        lst.append(val)

for line in lines:
    rec = line.rstrip('\n')
    if not rec.startswith(("ATOM", "HETATM")):
        continue
    atomname = rec[12:16]
    resname  = rec[17:20]
    chain    = rec[21]
    resseq   = rec[22:26]
    icode    = rec[26]
    occ_str  = rec[54:60]
    rkey = (resseq, icode, resname)
    akey = (chain, resseq, icode, resname, atomname)

    if rkey not in res_chains:
        res_order.append(rkey)
        res_chains[rkey] = []
        res_atoms[rkey]  = []
    _add_unique(res_chains[rkey], chain)
    _add_unique(res_atoms[rkey],  atomname)
    atom_rec[akey] = rec
    # store per-chain occ (any atom in this chain/residue will do)
    ckey = (chain, resseq, icode, resname)
    if ckey not in chain_occ:
        chain_occ[ckey] = occ_str

# ── 2. Fill missing atoms ──────────────────────────────────────────────────
n_filled = 0
for rkey in res_order:
    resseq, icode, resname = rkey
    chains = res_chains[rkey]
    if len(chains) < 2:
        continue
    for atomname in res_atoms[rkey]:
        # find the first chain that has this atom (source of coordinates)
        source = None
        for c in chains:
            sk = (c, resseq, icode, resname, atomname)
            if sk in atom_rec:
                source = atom_rec[sk]
                break
        if source is None:
            continue
        for chain in chains:
            akey = (chain, resseq, icode, resname, atomname)
            if akey in atom_rec:
                continue
            # use target chain's occupancy; fall back to source occ
            ckey = (chain, resseq, icode, resname)
            occ_str = chain_occ.get(ckey, source[54:60])
            # build new record: source coords, target chain letter + occupancy
            new_rec = source[:21] + chain + source[22:54] + occ_str + source[60:]
            atom_rec[akey] = new_rec
            chain_occ.setdefault(ckey, occ_str)
            n_filled += 1

print(f"Filled {n_filled} missing atoms", file=sys.stderr)

# ── 3. Write output ────────────────────────────────────────────────────────
# Header / REMARK lines
for line in lines:
    rec = line.rstrip('\n')
    if rec.startswith(("ATOM", "HETATM", "TER", "END", "MASTER")):
        break
    print(rec)

# Atoms: residue file order → chain file order → atom file order
atnum = 0
for rkey in res_order:
    resseq, icode, resname = rkey
    for chain in res_chains[rkey]:
        for atomname in res_atoms[rkey]:
            akey = (chain, resseq, icode, resname, atomname)
            if akey not in atom_rec:
                continue
            rec = atom_rec[akey]
            atnum += 1
            print(rec[:6] + f"{atnum:5d}" + rec[11:])

# Trailer
seen_end = False
for line in lines:
    rec = line.rstrip('\n')
    if rec.startswith(("TER", "MASTER", "CONECT")):
        print(rec)
    elif rec.startswith("END") and not seen_end:
        seen_end = True
        print("END")
if not seen_end:
    print("END")
