#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Shared RDKit fragment helpers for the S3 pipeline."""

from collections import deque

from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

def build_adj(mol):
    adj = {i: [] for i in range(mol.GetNumAtoms())}
    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        adj[a1].append(a2)
        adj[a2].append(a1)
    return adj


def k_hop(adj, seeds, k):
    visited = set(seeds)
    queue = deque((atom, 0) for atom in seeds)
    while queue:
        atom, depth = queue.popleft()
        if depth >= k:
            continue
        for neighbour in adj[atom]:
            if neighbour not in visited:
                visited.add(neighbour)
                queue.append((neighbour, depth + 1))
    return visited


def extract_fragments(mol, atom_set):
    editable = Chem.RWMol(mol)
    for idx in sorted(range(mol.GetNumAtoms()), reverse=True):
        if idx not in atom_set:
            editable.RemoveAtom(idx)
    fragment = editable.GetMol()
    try:
        Chem.SanitizeMol(fragment)
    except Exception:
        return None
    return Chem.MolToSmiles(fragment).split(".")


def ring_khop_fragments(smiles, k=2):
    if not isinstance(smiles, str) or not smiles.strip():
        return []
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    adj = build_adj(mol)
    fragments = []
    for ring in mol.GetRingInfo().AtomRings():
        parts = extract_fragments(mol, k_hop(adj, ring, k))
        if not parts:
            continue
        for part in parts:
            if Chem.MolFromSmiles(part):
                fragments.append(part)
    return list(set(fragments))


def longest_component(smiles):
    """Return the longest dot-separated component of a SMILES string."""
    if not isinstance(smiles, str) or not smiles.strip():
        return smiles
    parts = smiles.split(".")
    return max(parts, key=len) if len(parts) > 1 else smiles


def get_scaffold(smiles):
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
