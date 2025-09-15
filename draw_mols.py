#!/usr/bin/env python3
# draw_mols.py — read .smi and output a grid image of structures
from rdkit import Chem
from rdkit.Chem import Draw
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("--smi", required=True, help="Input .smi (SMILES [tab] Name)")
ap.add_argument("--out", default="molecule_grid.png", help="Output PNG")
ap.add_argument("--per-row", type=int, default=4, help="Molecules per row")
ap.add_argument("--size", type=int, default=300, help="Size of each cell (px)")
args = ap.parse_args()

smiles, names, mols = [], [], []
with open(args.smi, "r", encoding="utf-8") as f:
    for line in f:
        if not line.strip(): continue
        parts = line.strip().split()
        smi = parts[0]
        name = parts[1] if len(parts) > 1 else "mol"
        m = Chem.MolFromSmiles(smi)
        if m:
            mols.append(m)
            names.append(name)

if not mols:
    raise SystemExit("No molecules parsed.")

img = Draw.MolsToGridImage(
    mols, molsPerRow=args.per_row, subImgSize=(args.size, args.size),
    legends=names, useSVG=False
)
img.save(args.out)
print(f"Saved {args.out} ({len(mols)} molecules)")
