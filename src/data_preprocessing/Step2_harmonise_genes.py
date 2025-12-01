#!/usr/bin/env python

import os
import sys
from pathlib import Path
from typing import Optional, Dict

import pandas as pd


def _read_table_with_sep(path: Path, **kwargs) -> pd.DataFrame:
    """Read a CSV/TSV file and try to auto-detect the delimiter (tab or comma).

    Falls back to pandas default if detection fails.
    """
    # try to sniff first line
    with open(path, 'r', encoding=kwargs.pop('encoding', 'utf-8'), errors='replace') as fh:
        first = fh.readline()
    sep = '\t' if '\t' in first else ',' if ',' in first else None
    if sep is None:
        return pd.read_csv(path, engine='python', **kwargs)
    return pd.read_csv(path, sep=sep, engine='python', **kwargs)


def _find_gene_column(df: pd.DataFrame) -> Optional[str]:
    """Return the actual column name for 'Gene' in a case-insensitive way, or None."""
    for c in df.columns:
        if str(c).strip().lower() == 'gene':
            return c
    return None


def get_repo_root() -> Path:
    # This file lives in src/data_preprocessing -> src -> repo root
    return Path(__file__).resolve().parents[2]


def load_config():
    """
    Import config/config.py with minimal path hacking.
    """
    repo_root = get_repo_root()
    config_dir = repo_root / "config"
    sys.path.insert(0, str(config_dir))
    import config  # type: ignore  # noqa: E402

    return config


def build_symbol_map(hgnc_path: Path) -> Dict[str, str]:
    """
    Build a dict mapping any HGNC symbol / alias / previous symbol
    to the current approved HGNC symbol.

    Keys:
      - approved symbol
      - each 'prev_symbol' entry (split by '|')
      - each 'alias_symbol' entry (split by '|')

    Values:
      - approved HGNC symbol (row['symbol'])
    """
    if not hgnc_path.exists():
        raise FileNotFoundError(f"HGNC file not found at {hgnc_path}")

    hgnc = pd.read_csv(hgnc_path, sep="\t", dtype=str, low_memory=False)

    syn_to_approved: Dict[str, str] = {}

    def add_mapping(raw: Optional[str], approved: str) -> None:
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

        # approved symbol maps to itself
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


def harmonise_gene_column(df: pd.DataFrame, symbol_map: Dict[str, str], col: str = "Gene") -> pd.DataFrame:
    """
    Return a copy of df with df[col] harmonised to HGNC-approved symbols.

    Unmapped names are left as-is.
    """
    if col not in df.columns:
        raise KeyError(f"Column {col!r} not found in DataFrame")

    out = df.copy()

    # Vectorised map via dict; no explicit Python loop over rows
    def _map_one(g):
        if pd.isna(g):
            return g
        g = str(g).strip()
        return symbol_map.get(g, g)

    out[col] = out[col].map(_map_one)
    return out


def main():
    repo_root = get_repo_root()
    utils_dir = repo_root / "utils"
    hgnc_file = utils_dir / "hgnc_complete_set.txt"

    config = load_config()
    variant_dir = Path(getattr(config, "variant_output_directory"))

    if not hgnc_file.exists():
        raise SystemExit(
            f"ERROR: HGNC file not found at {hgnc_file}\n"
            "Download it first, e.g. from:\n"
            "  https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
        )

    if not variant_dir.exists():
        raise SystemExit(f"ERROR: variant_output_directory does not exist: {variant_dir}")

    print(f"Repo root         : {repo_root}")
    print(f"Variant directory : {variant_dir}")
    print(f"HGNC file         : {hgnc_file}")

    print("Building HGNC synonym map...")
    symbol_map = build_symbol_map(hgnc_file)
    print(f"Loaded {len(symbol_map)} synonym mappings")

    # Find annotated GWAS and variant_data files
    annotated_files = sorted(variant_dir.glob("Annotated_GWAS_*.csv"))
    variant_files = sorted(variant_dir.glob("variant_data_*.csv"))

    if not annotated_files and not variant_files:
        print("No Annotated_GWAS_*.csv or variant_data_*.csv found. Nothing to do.")
        return

    print("\n=== Harmonising Annotated_GWAS_* files ===")
    for path in annotated_files:
        print(f"  Processing {path.name} ...", end="", flush=True)
        try:
            df = _read_table_with_sep(path, dtype=str, keep_default_na=False)
        except Exception as e:
            print(f" [ERROR reading file: {e}]")
            continue

        gene_col = _find_gene_column(df)
        if gene_col is None:
            print(" [SKIP: no 'Gene' column]")
            continue

        before = df[gene_col].astype(str).copy()
        df_h = harmonise_gene_column(df, symbol_map, col=gene_col)
        after = df_h[gene_col].astype(str)

        n_changed = (before != after).sum()
        df_h.to_csv(path, index=False)
        print(f" [OK, {n_changed} gene names changed]")

    print("\n=== Harmonising variant_data_* files ===")
    for path in variant_files:
        print(f"  Processing {path.name} ...", end="", flush=True)
        try:
            df = _read_table_with_sep(path)
        except Exception as e:
            print(f" [ERROR reading file: {e}]")
            continue

        gene_col = _find_gene_column(df)
        if gene_col is None:
            print(" [SKIP: no 'Gene' column]")
            continue

        before = df[gene_col].astype(str).copy()
        df_h = harmonise_gene_column(df, symbol_map, col=gene_col)
        after = df_h[gene_col].astype(str)

        n_changed = (before != after).sum()
        df_h.to_csv(path, index=False)
        print(f" [OK, {n_changed} gene names changed]")

    print("\nDone. All applicable files overwritten with HGNC-harmonised Gene symbols.")


if __name__ == "__main__":
    main()
