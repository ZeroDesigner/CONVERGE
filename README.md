# CONVERGE

**Version:** v1.0.0  

## Introduction

This directory implements a fragment-level virtual screening pipeline for LolDF ATPase inhibitors. It combines RDKit fragment extraction, a Chemprop fragment classifier, and a QUBO-based diversity selection:

- `01_mol2frag.py` converts a molecule library into fragments. 
- `02_chemprop_predict.sh` scores candidate-library fragments with the trained
  model.
- `03_score_select.py` aggregates fragment probabilities into a molecule-level score (`avg_frag_score`, AFS), drops candidates with AFS <= 0.6, applies physicochemical filters (MW, LogP, HBA, HBD, QED), keeps the top 100 by AFS after scaffold deduplication, and selects 20 diverse molecules with QUBO.
- `utils.py` contains the shared RDKit fragment and scaffold helpers.


## Directory Tree


```
├── 01_mol2frag.py
├── 02_chemprop_predict.sh
├── 03_score_select.py
├── utils.py
├── CONVERGE.yml
├── chemprop.yml
├── Material/
├── OpenAnti/
└── fda_fragment_Chemprop/  
```

## Environments

Create the environments from the provided YAML files by [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html):

```bash
# Create and activate the CONVERGE environment
mamba env create -f CONVERGE.yml
mamba activate CONVERGE
# Create and activate the Chemprop environment
mamba env create -f chemprop.yml
mamba activate chemprop
```

## Tutorial: Screen a New Molecule Library

### 1. Prepare the input

Create a CSV with at least a `SMILES` column. A `Compound_id` column is recommended because it is carried through to the final table.

Example:

```csv
Compound_id,SMILES
ID_0001,CCOc1ccccc1
ID_0002,CC(=O)Nc1ccc(O)cc1
```

### 2: Molecules split into fragments

```bash
mamba activate CONVERGE
python 01_mol2frag.py \
  --input       ./Material/candidate_mols.csv \
  --output      ./Material/candidate_frags.csv
```

The output has one column, `frag_SMILES`, containing the unique ring-centred fragments.

### 3: Predict the fragments

```bash
mamba activate chemprop
./02_chemprop_predict.sh \
  --input_frags  ./Material/candidate_frags.csv \
  --input_models ./fda_fragment_Chemprop \
  --output ./Material/candidate_frags_pred.csv
```

The output has `frag_SMILES,label`, where `label` is the predicted probability.

### 4: Selection

```bash
mamba activate CONVERGE
python 03_score_select.py \
  --candidate_mols           ./Material/candidate_mols.csv \
  --candidate_frags          ./Material/candidate_frags.csv \
  --fragment_pred            ./Material/candidate_frags_pred.csv \
  --reference_openanti_dir   ./OpenAnti \
  --output                   ./Material/candidate_select.csv 
```

### 5. Inspect the result

Open `candidate_select.csv`. Each row is a selected molecule; `avg_frag_score` is the mean predicted probability of its fragments. 

#### Expected Output

Running the bundled sample data (`Material/candidate_mols.csv` → `candidate_select.csv`) is a self-check of the installation. The expected output is a table of **20 selected molecules** that  include:

| Compound_id | ID in Paper | SMILES                         |
| ----------- | -------------- | ------------------------------ |
| CMPD_0238   | LL236         | `OC1(C#Cc2ccccc2)CCCCCCCCCCC1` |

> Note: `Material/candidate_mols.csv` is a **sample test set** (1000 molecules) provided for verification. It is **not** training data and contains no experimental activity values.

## License

This project is licensed under the MIT License. See [MIT_license](./LICENSE) for the full license text.

## Data Sources and Acknowledgements

The reference library in `OpenAnti/` was compiled from publicly available data sources in September 2024, including [ChEMBL](https://www.ebi.ac.uk/chembl/) and [PubChem](https://pubchem.ncbi.nlm.nih.gov/).

We acknowledge the original data providers for making these resources available to the scientific community. 
