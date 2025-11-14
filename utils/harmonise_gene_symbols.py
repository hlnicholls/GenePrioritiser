#!/usr/bin/env python

import pandas as pd
from pathlib import Path

# Adjust to your repo layout: utils/hgnc_complete_set.txt
REPO_ROOT = Path(__file__).resolve().parents[2]
HGNC_FILE = REPO_ROOT / "utils" / "hgnc_complete_set.txt"


def build_symbol_map(hgnc_path: Path = HGNC_FILE) -> dict:
    """
    Build a dict mapping any HGNC symbol / alias / previous symbol
    to the current approved HGNC symbol.

    - Keys: approved symbol, alias_symbol entries, prev_symbol entries
    - Value: approved symbol
    """
    hgnc = pd.read_csv(hgnc_path, sep="\t", dtype=str, low_memory=False)

    syn_to_approved: dict[str, str] = {}

    def add_mapping(raw: str | None, approved: str) -> None:
        if raw is None:
            return
        raw = str(raw).strip()
        if not raw:
            return
        syn_to_approved[raw] = approved

    for _, row in hgnc.iterrows():
        approved = row.get("symbol")
        if pd.isna(approved):
            continue
        approved = str(approved).strip()
        if not approved:
            continue

        # approved symbol -> itself
        add_mapping(approved, approved)

        # previous symbols
        prev = row.get("prev_symbol")
        if isinstance(prev, str):
            for x in prev.split("|"):
                add_mapping(x, approved)

        # alias symbols
        alias = row.get("alias_symbol")
        if isinstance(alias, str):
            for x in alias.split("|"):
                add_mapping(x, approved)

    return syn_to_approved


def harmonise_gene_series(genes: pd.Series, symbol_map: dict) -> pd.Series:
    """
    Map each gene name to approved HGNC symbol if possible.
    Unmapped names are left untouched.
    """
    def _map_one(g):
        if pd.isna(g):
            return g
        g = str(g).strip()
        return symbol_map.get(g, g)

    # dict map is vectorised here; no Python-level for-loop over rows.
    return genes.map(symbol_map).fillna(genes)


def harmonise_gene_column(df: pd.DataFrame, col: str = "Gene",
                          symbol_map: dict | None = None) -> pd.DataFrame:
    """
    Return a copy of df with df[col] harmonised to HGNC-approved symbols.
    """
    if symbol_map is None:
        symbol_map = build_symbol_map()
    if col not in df.columns:
        raise KeyError(f"Column {col!r} not found in DataFrame")
    out = df.copy()
    out[col] = harmonise_gene_series(out[col], symbol_map)
    return out


if __name__ == "__main__":
    # tiny demo
    smap = build_symbol_map()
    test = pd.Series(["WAPAL", "WAPL", "AC104794.4"])
    print("raw   :", list(test))
    print("mapped:", list(harmonise_gene_series(test, smap)))
