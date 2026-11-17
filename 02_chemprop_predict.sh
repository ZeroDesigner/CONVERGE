#!/usr/bin/env bash
set -euo pipefail

INPUT_FRAGS=""
INPUT_MODELS=""
OUTPUT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input_frags)  INPUT_FRAGS="$2"; shift 2 ;;
    --input_models) INPUT_MODELS="$2"; shift 2 ;;
    --output)       OUTPUT="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

if [[ -z "$INPUT_FRAGS" || -z "$INPUT_MODELS" || -z "$OUTPUT" ]]; then
  echo "Usage: $0 --input_frags FRAGS.csv --input_models MODEL_DIR --output PRED.csv" >&2
  exit 2
fi

chemprop_predict \
    --test_path "$INPUT_FRAGS" \
    --checkpoint_dir "$INPUT_MODELS" \
    --preds_path "$OUTPUT" \
    --smiles_columns frag_SMILES \
    --features_generator rdkit_2d_normalized \
    --no_features_scaling
