# sglt2-rdkit
RDKit Analysis of SGLT2 Inhibitors
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
📂 Pipeline Description
A. Input

ChEMBL IDs of molecules (e.g., CHEMBL3039507 for sotagliflozin).

You can provide IDs in a text file (ids.txt) or directly via command line.

B. Scripts

sglt_pipeline.py → main pipeline script.

draw_mols.py → generates molecule grid image.

plot_props.py → generates histograms & scatter plots.

make_ppt.py → optional PowerPoint export.

C. Descriptors Calculated

Molecular Weight (MolWt)

Exact Molecular Weight (ExactMolWt)

LogP (Crippen cLogP)

Topological Polar Surface Area (TPSA)

H-bond donors (HBD) / acceptors (HBA)

Rotatable bonds (RotBonds)

Ring count (RingCount)

Fraction sp³ carbons (FractionCSP3)

Heavy atom count

Formal charge

QED (Quantitative Estimate of Drug-likeness)

▶️ Usage
1. Fetch SMILES + compute properties
# ids.txt contains one ChEMBL ID per line
python sglt_pipeline.py --ids-file ids.txt --out sglt_props.csv --save-smi sglt.smi --lipinski --veber


Outputs:

sglt.smi → SMILES + ChEMBL ID

sglt_props.csv → property table with Lipinski/Veber flags

2. Visualize molecules
python draw_mols.py --smi sglt.smi --out sglt_grid.png --per-row 4 --size 320

3. Generate property plots
python plot_props.py --csv sglt_props.csv


Creates:

MolWt_hist.png

cLogP_hist.png

TPSA_hist.png

RotBonds_hist.png

cLogP_vs_TPSA.png

MolWt_vs_TPSA.png

4. (Optional) Export to PowerPoint
python make_ppt.py --title "SGLT2 Inhibitors — RDKit Analysis" --grid sglt_grid.png --out Report.pptx

📊 Results & Interpretation
Structural Insights

All four inhibitors share a glucoside scaffold with variations in the aglycone region.

Substituents (ethoxy, cyclic ether, thiophene, fluorophenyl, sulfur groups) modulate potency and PK properties.

Physicochemical Properties

Molecular Weight: 408–451 Da (all <500; Lipinski compliant).

cLogP: 1.6–3.2 (moderate lipophilicity; favorable for oral drugs).

TPSA: 90–109 Å² (well below 140 Å²; supports good absorption).

Rotatable Bonds: 5–6 (below Veber’s cutoff of 10).

HBD/HBA: Balanced (3–4 donors, 6–7 acceptors).

QED: 0.49–0.66 (medium–high drug-likeness).

Visualizations

Histograms show clustering in drug-like ranges (MW, LogP, TPSA, flexibility).

Scatter plots confirm clustering in the “sweet spot” of oral bioavailability (LogP ~2, TPSA ~100).

Molecule grid illustrates scaffold conservation and chemical diversity.

Drug-likeness Evaluation

Lipinski’s Rule of Five: All pass.

Veber’s Rule: All pass.

Conclusion: Consistent with their clinical success as marketed oral drugs.

🚀 Implications

Confirms that RDKit pipelines can reproduce known drug-likeness profiles of approved drugs.

Provides a template workflow for analyzing other classes of molecules (from ChEMBL, ZINC20, etc.).

The clustering of all four drugs in a favorable physicochemical space highlights rational design principles for oral SGLT2 inhibitors.

This pipeline can be extended to:

Larger inhibitor libraries (e.g., full SGLT2 screening sets from ChEMBL/ZINC).

QSAR modeling.

Scaffold analysis and visualization.

Automated reporting (CSV + PPT).

📌 References

ChEMBL Database

ZINC20 Database

RDKit: Open-source cheminformatics, http://www.rdkit.org

👉 Would you like me to also prepare a short “Getting Started” demo (with the exact 4 SGLT2 inhibitors) inside the README, so visitors can run the pipeline without extra setup?

You said:
YES
ChatGPT said:

Perfect 🙌 — here’s the updated README.md draft with a clear Getting Started demo.
You can paste this directly into your GitHub repo.

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

▶️ Getting Started Demo

This demo reproduces results for four marketed SGLT2 inhibitors:

Sotagliflozin → CHEMBL3039507

Dapagliflozin → CHEMBL429910

Empagliflozin → CHEMBL2107830

Canagliflozin → CHEMBL4594217

Step 1 — Create ids.txt
CHEMBL3039507
CHEMBL429910
CHEMBL2107830
CHEMBL4594217

Step 2 — Run the pipeline
python sglt_pipeline.py --ids-file ids.txt --out sglt_props.csv --save-smi sglt.smi --lipinski --veber


Outputs:

sglt.smi → canonical SMILES + ChEMBL IDs

sglt_props.csv → computed descriptors + Lipinski/Veber flags

Step 3 — Visualize molecules
python draw_mols.py --smi sglt.smi --out sglt_grid.png --per-row 4 --size 320

Step 4 — Generate property plots
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
