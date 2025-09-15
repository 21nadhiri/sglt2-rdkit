#!/usr/bin/env python3
# plot_props.py — make simple, clean charts from sglt_props.csv
import argparse
import pandas as pd
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--csv", default="sglt_props.csv", help="Properties CSV")
args = ap.parse_args()

df = pd.read_csv(args.csv)

# A: Histograms (MolWt, cLogP, TPSA, RotBonds)
plots = [
    ("MolWt", "Molecular Weight"),
    ("cLogP", "cLogP (Crippen)"),
    ("TPSA", "Topological Polar Surface Area"),
    ("RotBonds", "Rotatable Bonds"),
]
for col, label in plots:
    plt.figure()
    df[col].hist(bins=10, edgecolor="black")
    plt.title(label)
    plt.xlabel(label)
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(f"{col}_hist.png", dpi=200)
    plt.close()

# B: Scatter charts (cLogP vs TPSA; MolWt vs TPSA)
pairs = [
    ("cLogP", "TPSA", "cLogP_vs_TPSA.png"),
    ("MolWt", "TPSA", "MolWt_vs_TPSA.png"),
]
for x, y, fn in pairs:
    plt.figure()
    plt.scatter(df[x], df[y])
    for i, row in df.iterrows():
        plt.annotate(str(row["Name"]), (row[x], row[y]), fontsize=8, xytext=(3,3), textcoords="offset points")
    plt.xlabel(x)
    plt.ylabel(y)
    plt.title(f"{x} vs {y}")
    plt.tight_layout()
    plt.savefig(fn, dpi=200)
    plt.close()

print("Saved charts: MolWt_hist.png, cLogP_hist.png, TPSA_hist.png, RotBonds_hist.png, cLogP_vs_TPSA.png, MolWt_vs_TPSA.png")
