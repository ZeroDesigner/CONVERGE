# CONVERGE

## Introduction

This directory implements a fragment-level virtual screening pipeline for LolDF ATPase inhibitors. It combines RDKit fragment extraction, a Chemprop fragment classifier, and a QUBO-based diversity selection:

- `01_mol2frag.py` converts a molecule library into fragments. 
- `02_chemprop_predict.sh` scores candidate-library fragments with the trained
  model.
- `03_score_select.py` aggregates fragment probabilities into a molecule-level score (`avg_frag_score`, AFS), applies physicochemical and scaffold filters, and selects a diverse set of molecules with QUBO.
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

Create a CSV with at least a `SMILES` column. A `Compound_id` column is
recommended because it is carried through to the final table.

Example:

```csv
Compound_id,SMILES
CMPD_0001,CCOc1ccccc1
CMPD_0002,CC(=O)Nc1ccc(O)cc1
```

### 2: molecules split into fragments

```bash
mamba activate CONVERGE
python 01_mol2frag.py \
  --input       ./Material/candidate_mols.csv \
  --output      ./Material/candidate_frags.csv

```

The output has one column, `frag_SMILES`, containing the unique ring-centred fragments.

### 3: predict the fragments

```bash
mamba activate chemprop
./02_chemprop_predict.sh \
  --input_frags  ./Material/candidate_frags.csv \
  --input_models ./fda_fragment_Chemprop \
  --output ./Material/candidate_frags_pred.csv

```

The output has `frag_SMILES,label`, where `label` is the predicted probability.

### 4: selection

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

## License

This project is licensed under the MIT License. See [MIT_license](./MIT_license) for the full license text.



The other materials are coming soon.
