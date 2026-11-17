#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stage 3: score molecules from fragment predictions and select ~K diverse
molecules with a QUBO.
"""

import argparse
import multiprocessing as mp
import random
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import openjij as oj
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors, QED

from utils import get_scaffold, longest_component, ring_khop_fragments


_PARALLEL_Q = None

def _sample_seed(seed):
    response = oj.SASampler().sample_qubo(
        _PARALLEL_Q, num_reads=500, seed=seed
    )
    sample = response.first
    return float(sample.energy), dict(sample.sample)


def make_sample_seeds(n=200):
    return random.sample(range(1, 100001), n)


# ==============================================================================
# Fragment and property helpers
# ==============================================================================

def calc_props(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None,) * 5
    try:
        return (
            Descriptors.MolWt(mol),
            Descriptors.MolLogP(mol),
            Descriptors.NumHAcceptors(mol),
            Descriptors.NumHDonors(mol),
            QED.qed(mol),
        )
    except Exception:
        return (None,) * 5


def get_frag_fingerprint(fragments):
    if not fragments:
        return None
    fingerprints = []
    for smiles in fragments:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fingerprints.append(
                AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            )
    return fingerprints or None


def open_anti(opendir):
    opendir = Path(opendir)
    keep_cols = ["ChEMBL ID", "Smiles"]
    df_c1 = pd.read_csv(
        opendir / "ChemBL_Antibacterial.csv", sep=";", usecols=keep_cols
    )
    df_c2 = pd.read_csv(
        opendir / "ChemBL_Antibiotic.csv", sep=";", usecols=keep_cols
    )
    df_p = pd.read_csv(
        opendir / "PubChem_compound_Antibiotic.csv",
        usecols=["Compound_CID", "SMILES"],
    )

    df_chembl = pd.concat([df_c1, df_c2], ignore_index=True)
    df_chembl["Library"] = "ChemBL"
    df_chembl["ID"] = df_chembl["ChEMBL ID"]
    df_chembl["SMILES"] = df_chembl["Smiles"]
    df_p["Library"] = "PubChem"
    df_p["ID"] = df_p["Compound_CID"]

    df_all = pd.concat(
        [
            df_chembl[["Library", "ID", "SMILES"]],
            df_p[["Library", "ID", "SMILES"]],
        ],
        ignore_index=True,
    )
    df_all["SMILES"] = df_all["SMILES"].astype(str)
    df_all["Scaffold"] = df_all["SMILES"].apply(get_scaffold)
    df_all = df_all[df_all["Scaffold"].map(lambda x: bool(x))].reset_index(drop=True)
    df_all["Lib"] = "OpenAnti"
    df_all["Compound_id"] = df_all["ID"]
    return df_all[["Lib", "Compound_id", "SMILES", "Scaffold"]]


# ==============================================================================
# QUBO selection
# ==============================================================================

def qubo_select_molecules_by_frags(
    df_all,
    df_align_frag,
    col_name="avg_frag_score",
    K=20,
    alpha=0.5,
    beta=0.5,
    lam1=1.0,  # activity
    lam2=1.0,  # structural redundancy
    lam3=1.0,  # reference alignment
    lam4=2.0,  # cardinality
    jobs=10,
):
    df_all = df_all.copy().reset_index(drop=True)
    N = len(df_all)
    if N == 0:
        return df_all, [], np.zeros(0)

    df_all["frag_fps"] = df_all["fragments"].apply(get_frag_fingerprint)
    align_fps = []
    for smiles in df_align_frag["frag_SMILES"]:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            align_fps.append(
                AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            )

    S = np.zeros((N, N))
    for i in range(N):
        fp_list_i = df_all.at[i, "frag_fps"]
        if fp_list_i is None:
            continue
        for j in range(i + 1, N):
            fp_list_j = df_all.at[j, "frag_fps"]
            if fp_list_j is None:
                continue
            max_sim = 0.0
            for fi in fp_list_i:
                for fj in fp_list_j:
                    sim = DataStructs.TverskySimilarity(fi, fj, alpha, beta)
                    if sim > max_sim:
                        max_sim = sim
            S[i, j] = S[j, i] = max_sim

    S_align = np.zeros(N)
    for i in range(N):
        fp_list_i = df_all.at[i, "frag_fps"]
        if fp_list_i is None or not align_fps:
            continue
        max_sim_ext = 0.0
        for fi in fp_list_i:
            sims = DataStructs.BulkTverskySimilarity(fi, align_fps, alpha, beta)
            local_max = max(sims)
            if local_max > max_sim_ext:
                max_sim_ext = local_max
        S_align[i] = max_sim_ext

    cat_score = df_all[col_name].fillna(0).values
    if cat_score.max() > cat_score.min():
        cat_score = (cat_score - cat_score.min()) / (
            cat_score.max() - cat_score.min()
        )

    Q = {}
    for i in range(N):
        Q[(i, i)] = (
            lam3 * S_align[i] + lam4 * (1 - 2 * K) - lam1 * cat_score[i]
        )
        for j in range(i + 1, N):
            value = lam2 * S[i, j] + 2 * lam4
            if value != 0:
                Q[(i, j)] = value

    sample_seeds = make_sample_seeds()
    # print(f"Sample seeds: {sample_seeds}")
    global _PARALLEL_Q
    _PARALLEL_Q = Q
    with ProcessPoolExecutor(
        max_workers=min(jobs, len(sample_seeds)),
        mp_context=mp.get_context("fork"),
    ) as executor:
        raw = list(executor.map(_sample_seed, sample_seeds, chunksize=1))
    records = []
    for current_seed, (energy, sample_dict) in zip(sample_seeds, raw):
        idx = [i for i, value in sample_dict.items() if value == 1]
        records.append({
            "seed": current_seed,
            "energy": float(energy),
            "n": len(idx),
            "idx": idx,
        })
        # print(
        #     f"  seed={current_seed} n={len(idx)} "
        #     f"energy={energy:.3f}"
        # )
    candidates = [r for r in records if r["n"] == K]
    #print(f"Seeds sampled: {len(records)}, n=={K}: {len(candidates)}")
    if not candidates:
        print(
            f"  WARNING: no solution with exactly {K} molecules "
            "among the sampled seeds"
        )
        return df_all.iloc[[]].copy(), [], S_align
    best = min(candidates, key=lambda r: r["energy"])
    selected_idx = list(best["idx"])
    # print(
    #     f"  final: seed={best['seed']} n={len(selected_idx)} "
    #     f"energy={best['energy']:.3f}"
    # )
    return df_all.iloc[selected_idx].copy(), selected_idx, S_align


# ==============================================================================
# CLI
# ==============================================================================

def build_parser():
    parser = argparse.ArgumentParser(
        description="Stage 3: fragment scores -> AFS scores -> QUBO selection."
    )
    parser.add_argument(
        "--candidate_mols", "--candidate-mols", dest="candidate_mols", required=True
    )
    parser.add_argument(
        "--candidate_frags", "--candidate-frags", dest="candidate_frags", default=None
    )
    parser.add_argument(
        "--fragment_pred", "--fragment-pred", dest="fragment_pred", default=None
    )
    parser.add_argument(
        "--reference_training", "--reference-training", dest="reference_training", default=None
    )
    parser.add_argument(
        "--reference_openanti_dir",
        "--reference-openanti-dir",
        dest="reference_openanti_dir",
        default=None,
    )
    parser.add_argument("--output", required=True)
    parser.add_argument("--smiles_col", "--smiles-col", dest="smiles_col", default="SMILES")
    parser.add_argument("--scaffold_col", "--scaffold-col", dest="scaffold_col", default="Scaffold")
    parser.add_argument("--score_col", "--score-col", dest="score_col", default="label")
    parser.add_argument(
        "--reference_label_col",
        "--reference-label-col",
        dest="reference_label_col",
        default="Label",
    )
    parser.add_argument("--k", type=int, default=2)
    parser.add_argument("--top_n", "--top-n", dest="top_n", type=int, default=100)
    parser.add_argument("--min_mw", "--min-mw", dest="min_mw", type=float, default=200.0)
    parser.add_argument("--min_logp", "--min-logp", dest="min_logp", type=float, default=3.0)
    parser.add_argument("--min_hba", "--min-hba", dest="min_hba", type=float, default=1.0)
    parser.add_argument("--min_hbd", "--min-hbd", dest="min_hbd", type=float, default=1.0)
    parser.add_argument("--min_qed", "--min-qed", dest="min_qed", type=float, default=0.5)
    parser.add_argument("--select_k", "--select-k", dest="select_k", type=int, default=20)
    parser.add_argument(
        "--lam1", type=float, default=1.0,
        help="weight for AI-predicted activity",
    )
    parser.add_argument(
        "--lam2", type=float, default=1.0,
        help="weight for pairwise structural redundancy",
    )
    parser.add_argument(
        "--lam3", type=float, default=1.0,
        help="weight for reference-space alignment",
    )
    parser.add_argument(
        "--lam4", type=float, default=2.0,
        help="weight for the cardinality penalty",
    )
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--beta", type=float, default=0.5)
    parser.add_argument(
        "--jobs", type=int, default=10,
        help="worker processes for seed sampling (default: 10)",
    )
    return parser


def reference_fragments(df, smiles_col, scaffold_col, label_col, positive, k):
    if label_col in df.columns:
        df = df[pd.to_numeric(df[label_col], errors="coerce") == positive]

    fragments = set()
    for _, row in df.iterrows():
        scaffold = row[scaffold_col] if scaffold_col in df.columns else None
        if not isinstance(scaffold, str) or not scaffold.strip():
            smiles = row[smiles_col] if smiles_col in df.columns else None
            scaffold = get_scaffold(smiles)
        if isinstance(scaffold, str) and scaffold:
            fragments.update(ring_khop_fragments(scaffold, k))
    return fragments


def main():
    args = build_parser().parse_args()

    fragment_pred = args.fragment_pred
    if not fragment_pred:
        if not args.candidate_frags:
            raise SystemExit(
                "--fragment_pred is required (or provide --candidate_frags)"
            )
        path = Path(args.candidate_frags)
        fragment_pred = str(path.with_name(path.stem + "_pred.csv"))

    predictions = pd.read_csv(fragment_pred)
    if "frag_SMILES" not in predictions.columns:
        raise SystemExit(f"Missing frag_SMILES column in {fragment_pred}")
    if args.score_col not in predictions.columns:
        raise SystemExit(f"Missing {args.score_col} column in {fragment_pred}")

    pred_scores = dict(
        zip(
            predictions["frag_SMILES"].astype(str),
            pd.to_numeric(
                predictions[args.score_col], errors="coerce"
            ).fillna(0.0),
        )
    )

    if args.candidate_frags:
        candidate_frags = pd.read_csv(args.candidate_frags)
        if "frag_SMILES" in candidate_frags.columns:
            candidate_set = set(candidate_frags["frag_SMILES"].astype(str))
            covered = len(candidate_set & set(pred_scores))
            print(
                f"Candidate fragments: {len(candidate_set)}, "
                f"with predictions: {covered}"
            )

    df = pd.read_csv(args.candidate_mols)
    if args.smiles_col not in df.columns:
        raise SystemExit(f"Missing SMILES column: {args.smiles_col}")
    df = df[df[args.smiles_col].notna()].copy()
    df[args.smiles_col] = df[args.smiles_col].astype(str)

    scaffolds = (
        df[args.scaffold_col].tolist()
        if args.scaffold_col in df.columns
        else [None] * len(df)
    )

    resolved_scaffolds = []
    fragment_lists = []
    for smiles, scaffold in zip(df[args.smiles_col], scaffolds):
        if not isinstance(scaffold, str) or not scaffold.strip():
            scaffold = get_scaffold(longest_component(smiles))
        resolved_scaffolds.append(scaffold)
        fragment_lists.append(
            ring_khop_fragments(scaffold, args.k)
            if isinstance(scaffold, str) and scaffold
            else []
        )

    df["Scaffold"] = resolved_scaffolds
    df["fragments"] = fragment_lists

    properties = [calc_props(smiles) for smiles in df[args.smiles_col]]
    df["MolWt"] = [value[0] for value in properties]
    df["LogP"] = [value[1] for value in properties]
    df["HBA"] = [value[2] for value in properties]
    df["HBD"] = [value[3] for value in properties]
    df["QED"] = [value[4] for value in properties]

    df["fragment_scores"] = df["fragments"].apply(
        lambda fragments: [pred_scores.get(fragment, 0.0) for fragment in fragments]
    )
    df["avg_frag_score"] = df["fragment_scores"].apply(
        lambda values: float(np.mean(values)) if values else 0.0
    )

    print(f"Before AFS cutoff: {len(df)}")
    df = df[df["avg_frag_score"] > 0.6].copy()
    print(f"After AFS > 0.6: {len(df)}")

    df = df[
        (df["MolWt"] >= args.min_mw)
        & (df["LogP"] >= args.min_logp)
        & (df["HBA"] >= args.min_hba)
        & (df["HBD"] >= args.min_hbd)
        & (df["QED"] >= args.min_qed)
    ].copy()
    print(f"After property filter: {len(df)}")

    df = (
        df.sort_values("avg_frag_score", ascending=False)
        .drop_duplicates(subset=["Scaffold"], keep="first")
        .head(args.top_n)
        .reset_index(drop=True)
    )
    print(f"After scaffold dedupe/top{args.top_n}: {len(df)}")

    align_frags = set()
    if args.reference_training:
        reference = pd.read_csv(args.reference_training)
        align_frags.update(
            reference_fragments(
                reference,
                args.smiles_col,
                args.scaffold_col,
                args.reference_label_col,
                1,
                args.k,
            )
        )
    if args.reference_openanti_dir:
        openanti = open_anti(args.reference_openanti_dir)
        for scaffold in openanti["Scaffold"]:
            align_frags.update(ring_khop_fragments(scaffold, args.k))
    df_align_frag = pd.DataFrame({"frag_SMILES": sorted(align_frags)})

    selected, _, _ = qubo_select_molecules_by_frags(
        df_all=df,
        df_align_frag=df_align_frag,
        col_name="avg_frag_score",
        K=args.select_k,
        alpha=args.alpha,
        beta=args.beta,
        lam1=args.lam1,  # activity
        lam2=args.lam2,  # structural redundancy
        lam3=args.lam3,  # reference alignment
        lam4=args.lam4,  # cardinality
        jobs=args.jobs,
    )
    selected = selected.drop(columns=["frag_fps"], errors="ignore")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.output, index=False)
    print(
        f"Candidates: {len(df)}, selected: {len(selected)}, "
        f"reference fragments: {len(align_frags)}"
    )
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
