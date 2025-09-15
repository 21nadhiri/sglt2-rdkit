#!/usr/bin/env python3
# chembl_to_smi.py
import argparse, sys, time
from pathlib import Path
import requests

API = "https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"

def get_smiles(chembl_id: str) -> str | None:
    url = API.format(chembl_id=chembl_id)
    r = requests.get(url, timeout=20)
    if r.status_code != 200:
        sys.stderr.write(f"[WARN] {chembl_id}: HTTP {r.status_code}\n")
        return None
    data = r.json()
    ms = data.get("molecule_structures") or {}
    smi = ms.get("canonical_smiles")
    if not smi:
        sys.stderr.write(f"[WARN] {chembl_id}: no canonical_smiles\n")
    return smi

def main():
    ap = argparse.ArgumentParser(description="Fetch canonical SMILES from ChEMBL IDs")
    ap.add_argument("--ids", nargs="*", help="ChEMBL IDs (space-separated)")
    ap.add_argument("--ids-file", help="Path to a file with one ChEMBL ID per line")
    ap.add_argument("--out", default="chembl.smi", help="Output .smi (SMILES [tab] NAME)")
    ap.add_argument("--sleep", type=float, default=0.2, help="Sleep seconds between requests")
    args = ap.parse_args()

    ids = []
    if args.ids_file:
        ids += [line.strip() for line in Path(args.ids_file).read_text().splitlines() if line.strip()]
    if args.ids:
        ids += args.ids
    if not ids:
        print("Provide --ids or --ids-file")
        sys.exit(1)

    with open(args.out, "w", encoding="utf-8") as f:
        for cid in ids:
            smi = get_smiles(cid)
            if smi:
                # Write SMILES + name (here we use ChEMBL ID as the name)
                f.write(f"{smi}\t{cid}\n")
            time.sleep(args.sleep)

    print(f"Wrote {args.out} with {len(ids)} entries (check warnings for any missing).")

if __name__ == "__main__":
    main()

