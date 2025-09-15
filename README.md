# sglt2-rdkit
RDKit Analysis of SGLT2 Inhibitors

📊 RDKit Analysis of SGLT2 Inhibitors
🔎 Project Overview

This project uses RDKit to compute molecular properties and drug-likeness of SGLT2 inhibitors (sotagliflozin, dapagliflozin, empagliflozin, and canagliflozin).
The pipeline automatically:

Fetches canonical SMILES from ChEMBL API using ChEMBL IDs.

Standardizes molecules with RDKit (cleanup, largest fragment, normalize, uncharge).

Calculates a set of medicinal chemistry descriptors.

Applies Lipinski and Veber filters for drug-likeness.

Outputs structured results (.smi, .csv) and visualization plots (molecule grid, histograms, scatter plots).

Provides interpretation of results for scientific presentation.

⚙️ Installation
1. Clone this repo
git clone https://github.com/yourusername/sglt2-rdkit.git
cd sglt2-rdkit

2. Create environment with Conda (recommended)
conda create -n rdkit -c conda-forge rdkit pandas matplotlib pillow python-pptx requests -y
conda activate rdkit

This is four marketed SGLT2 inhibitors:

Sotagliflozin → CHEMBL3039507

Dapagliflozin → CHEMBL429910

Empagliflozin → CHEMBL2107830

Canagliflozin → CHEMBL4594217

Step 1  Create ids.txt
CHEMBL3039507
CHEMBL429910
CHEMBL2107830
CHEMBL4594217

Step 2  Run the pipeline
python sglt_pipeline.py --ids-file ids.txt --out sglt_props.csv --save-smi sglt.smi --lipinski --veber


Outputs:

sglt.smi → canonical SMILES + ChEMBL IDs

sglt_props.csv → computed descriptors + Lipinski/Veber flags

Step 3 Visualize molecules
python draw_mols.py --smi sglt.smi --out sglt_grid.png --per-row 4 --size 320

Step 4  Generate property plots
python plot_props.py --csv sglt_props.csv


Creates:

MolWt_hist.png

cLogP_hist.png

TPSA_hist.png

RotBonds_hist.png

cLogP_vs_TPSA.png

MolWt_vs_TPSA.png

Step 5 — (Optional) Export to PowerPoint
python make_ppt.py --title "SGLT2 Inhibitors — RDKit Analysis" --grid sglt_grid.png --out Report.pptx

📂 Pipeline Description

Input: ChEMBL IDs or a .smi file.

Core Script: sglt_pipeline.py

Additional Scripts:

draw_mols.py → molecule grid PNG

plot_props.py → histograms & scatter plots

make_ppt.py → slide deck

Descriptors calculated: MolWt, ExactMolWt, cLogP, TPSA, HBD, HBA, RotBonds, RingCount, FractionCSP3, HeavyAtomCount, FormalCharge, QED.

📊 Results & Interpretation

Molecular Weight: 408–451 Da → within Lipinski’s cutoff (<500).

cLogP: 1.6–3.2 → moderate lipophilicity (good balance).

TPSA: 90–109 Å² → <140 Å² (favorable absorption).

Rotatable Bonds: 5–6 → <10 (moderate flexibility).

Lipinski & Veber Rules: All compounds pass.

QED: 0.49–0.66 (good drug-likeness).

Visualizations:

Histograms show clustering in drug-like ranges.

Scatter plots (cLogP vs TPSA; MolWt vs TPSA) confirm clustering in oral drug “sweet spot.”

Molecule grid highlights structural similarity (glucoside core) with diverse substituents.

🚀 Implications

Confirms why all four are approved oral drugs: they occupy favorable chemical space.

Pipeline can be applied to other ChEMBL/ZINC datasets for:

hit-to-lead profiling,

QSAR modeling,

scaffold analysis,

automated reporting (CSV + PPT).

📌 References

ChEMBL Database

ZINC20 Database

RDKit Documentation

