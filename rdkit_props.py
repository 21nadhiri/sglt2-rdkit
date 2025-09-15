#!/usr/bin/env python3
# rdkit_props.py
import argparse, os, sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen, Lipinski, QED
from rdkit.Chem.MolStandardize import rdMolStandardize

def standardize(m: Chem.Mol) -> Chem.Mol | None:
    if m is None: return None
    m = rdMolStandardize.Cleanup(m)
    m = rdMolStandardize.LargestFragmentChooser().choose(m)
    m = rdMolStandardize.Normalize(m)
    m = rdMolStandardize.Uncharger().uncharge(m)
    Chem.SanitizeMol(m)
    return m

def calc_props(m: Chem.Mol) -> dict:
    return dict(
        InChIKey=Chem.inchi.MolToInchiKey(m),
        MolWt=Descriptors.MolWt(m),
        ExactMolWt=Descriptors.ExactMolWt(m),
        cLogP=Crippen.MolLogP(m),
        TPSA=rdMolDescriptors.CalcTPSA(m),
        HBD=Lipinski.NumHDonors(m),
        HBA=Lipinski.NumHAcceptors(m),
        RotBonds=Lipinski.NumRotatableBonds(m),
        RingCount=rdMolDescriptors.CalcNumRings(m),
        FractionCSP3=rdMolDescriptors.CalcFractionCSP3(m),
        HeavyAtomCount=Descriptors.HeavyAtomCount(m),
        FormalCharge=Chem.GetFormalCharge(m),
        QED=QED.qed(m),
    )

def lipinski_ok(r) -> bool:
    return (r["MolWt"] <= 500 and r["cLogP"] <= 5 and r["HBD"] <= 5 and r["HBA"] <= 10)

def veber_ok(r) -> bool:
    return (r["TPSA"] <= 140 and r["RotBonds"] <= 10)

def load_smi(path: str):
    # Expects: SMILES [tab/space] NAME (NAME optional)
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith("#"): continue
            parts = line.split()
            smi = parts[0]
            name = parts[1] if len(parts) > 1 else os.path.basename(path)
            m = Chem.MolFromSmiles(smi)
            if m: yield name, m
            else: sys.stderr.write(f"[WARN] bad SMILES: {line}\n")

def load_sdf(path: str):
    supp = Chem.SDMolSupplier(path, removeHs=False)
    for m in supp:
        if m is None: continue
        name = m.GetProp("_Name") if m.HasProp("_Name") else os.path.basename(path)
        yield name, m

def main():
    ap = argparse.ArgumentParser(description="RDKit standardize + property calculator")
    ap.add_argument("--input", "-i", required=True, nargs="+", help="Input .smi/.txt or .sdf")
    ap.add_argument("--out", "-o", default="properties.csv", help="Output CSV")
    ap.add_argument("--lipinski", action="store_true", help="Add Lipinski_OK column")
    ap.add_argument("--veber", action="store_true", help="Add Veber_OK column")
    args = ap.parse_args()

    rows, seen = [], set()
    for path in args.input:
        if not os.path.isfile(path):
            sys.stderr.write(f"[WARN] not found: {path}\n")
            continue
        ext = os.path.splitext(path)[1].lower()
        loader = load_sdf if ext == ".sdf" else load_smi
        for name, m in loader(path):
            try:
                m = standardize(m)
            except Exception as e:
                sys.stderr.write(f"[WARN] standardize failed for {name}: {e}\n")
                continue
            if m is None: continue
            ik = Chem.inchi.MolToInchiKey(m)
            if ik in seen:  # de-duplicate by structure
                continue
            seen.add(ik)
            row = {"Name": name} | calc_props(m)
            rows.append(row)

    if not rows:
        print("No valid molecules processed.")
        sys.exit(2)

    df = pd.DataFrame(rows)
    if args.lipinski:
        df["Lipinski_OK"] = df.apply(lipinski_ok, axis=1)
    if args.veber:
        df["Veber_OK"] = df.apply(veber_ok, axis=1)

    df.to_csv(args.out, index=False)
    print(f"Saved {len(df)} rows -> {args.out}")

if __name__ == "__main__":
    main()
