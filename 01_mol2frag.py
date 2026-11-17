#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stage 1: convert a molecule library to fragments.
"""

import argparse
from pathlib import Path

import pandas as pd

from utils import get_scaffold, longest_component, ring_khop_fragments


def build_parser():
    parser = argparse.ArgumentParser(
        description="Stage 1: molecule library -> fragment table."
    )
    parser.add_argument("--input", required=True, help="Input molecule CSV.")
    parser.add_argument("--output", required=True, help="Output fragment CSV.")
    parser.add_argument("--smiles_col", "--smiles-col", dest="smiles_col", default="SMILES")
    parser.add_argument("--scaffold_col", "--scaffold-col", dest="scaffold_col", default="Scaffold")
    parser.add_argument("--k", type=int, default=2, help="Ring-centred hop count.")
    return parser


def main():
    args = build_parser().parse_args()

    df = pd.read_csv(args.input)
    if args.smiles_col not in df.columns:
        raise SystemExit(f"Missing SMILES column: {args.smiles_col}")

    df = df[df[args.smiles_col].notna()].copy()
    df[args.smiles_col] = df[args.smiles_col].astype(str)

    if args.scaffold_col in df.columns:
        scaffolds = df[args.scaffold_col].tolist()
    else:
        scaffolds = [None] * len(df)

    fragment_lists = []
    for smiles, scaffold in zip(df[args.smiles_col], scaffolds):
        if not isinstance(scaffold, str) or not scaffold.strip():
            scaffold = get_scaffold(longest_component(smiles))
        fragment_lists.append(
            ring_khop_fragments(scaffold, args.k)
            if isinstance(scaffold, str) and scaffold
            else []
        )

    fragments = sorted(
        {fragment for row in fragment_lists for fragment in row}
    )
    output = pd.DataFrame({"frag_SMILES": fragments})

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)
    print(f"Wrote {args.output}: {len(output)} rows")


if __name__ == "__main__":
    main()
