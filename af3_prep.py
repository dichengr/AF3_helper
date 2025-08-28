#!/usr/bin/env python3
# af3_prep.py
# Create AlphaFold3 input JSONs from an input table (supports PTMs, mutations, and homomeric copies).

from pathlib import Path
import argparse
import json
import sys
import re
from typing import List, Tuple, Optional

import pandas as pd
from tqdm import tqdm

from prep_utils import (
    fetch_uniprot_sequence,
    parse_ptms,
    parse_mutations,
    apply_mutations,
    format_complex_name,
    build_alphafold_json,
)

def _parse_model_seeds(seed_spec: str) -> List[int]:
    seeds: List[int] = []
    for part in seed_spec.split(','):
        part = part.strip()
        if not part:
            continue
        if '-' in part:
            start, end = (int(x) for x in part.split('-', 1))
            seeds.extend(range(start, end + 1))
        else:
            seeds.append(int(part))
    return seeds

def _read_input_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        sys.exit(f"Input file '{path}' not found")

    if path.suffix.lower() == '.xlsx':
        df = pd.read_excel(path, engine='openpyxl')
    else:
        df = pd.read_csv(path, sep=None, engine='python')

    df.columns = [str(c).strip().lower() for c in df.columns]
    return df

def _has_long_format(df: pd.DataFrame) -> bool:
    cols = set(df.columns)
    return {'complex', 'name', 'uniprot'} <= cols

def _collect_wide_columns(df: pd.DataFrame) -> List[Tuple[str, str, Optional[str]]]:
    triples: List[Tuple[str, str, Optional[str]]] = []
    for i in range(1, 100):
        name_col = f'name{i}'
        acc_col = f'acc{i}' if f'acc{i}' in df.columns else f'uniprot{i}'
        ptm_col = f'ptms{i}' if f'ptms{i}' in df.columns else (f'ptm{i}' if f'ptm{i}' in df.columns else None)
        if name_col in df.columns and acc_col in df.columns:
            triples.append((name_col, acc_col, ptm_col))
        else:
            break
    if not triples:
        sys.exit('Wide format detected but could not find columns like name1/acc1/name2/acc2')
    return triples

def _fetch_sequences(accs: List[str]) -> List[str]:
    seqs: List[str] = []
    for acc in accs:
        seq = fetch_uniprot_sequence(str(acc))
        seqs.append(seq or '')
    return seqs

def _group_long(df: pd.DataFrame):
    for complex_id, subdf in df.groupby('complex', sort=False):
        yield complex_id, subdf.reset_index(drop=True)

def _excel_col(idx: int) -> str:
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    idx += 1
    s = ''
    while idx > 0:
        idx, rem = divmod(idx - 1, 26)
        s = letters[rem] + s
    return s

def _assign_chain_ids_by_copies(copies: List[int]) -> List[List[str]]:
    total = sum(max(1, int(c or 1)) for c in copies)
    flat = [_excel_col(i) for i in range(total)]
    out: List[List[str]] = []
    k = 0
    for c in copies:
        c = max(1, int(c or 1))
        out.append(flat[k:k+c])
        k += c
    return out

def _render_filename(base_dir: Path, payload_name: str) -> Path:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", payload_name)
    return base_dir / f"{safe}.json"

def main() -> None:
    parser = argparse.ArgumentParser(
        description='Create AlphaFold3 input JSON from an input table (supports PTMs, mutations, and copies).',
    )
    parser.add_argument('-i', '--input', required=True,
                        help='Input table (.xlsx/.csv/.tsv). Long format requires: complex,name,uniprot[,ptms,muts,copies].')
    parser.add_argument('-o', '--output', help='Output directory for JSON files (default: current dir)')
    parser.add_argument('-s', '--seeds', default='1-20',
                        help='Comma-separated list or ranges of model seeds (e.g. 1-4,6,8)')

    args = parser.parse_args()

    df = _read_input_table(Path(args.input))
    model_seeds = _parse_model_seeds(args.seeds)

    out_dir = Path(args.output) if args.output else Path('.')
    out_dir.mkdir(parents=True, exist_ok=True)

    if _has_long_format(df):
        total_complexes = df['complex'].nunique()
        print(f'Processing {total_complexes} complexes in long format...')

        for complex_id, sub in tqdm(_group_long(df), total=total_complexes, desc='Complexes'):
            names = [str(x) for x in sub['name'].tolist()]
            accs = [str(x) for x in sub['uniprot'].tolist()]

            copies = [int(x) if str(x).strip().isdigit() else 1 for x in (sub['copies'].tolist() if 'copies' in sub.columns else [1]*len(names))]
            chain_ids = _assign_chain_ids_by_copies(copies)

            ptm_fields = sub['ptms'].astype(str).tolist() if 'ptms' in sub.columns else [''] * len(names)

            sequences = _fetch_sequences(accs)

            mut_fields = sub['muts'].astype(str).tolist() if 'muts' in sub.columns else [''] * len(names)
            mutation_tags_per_chain: List[List[str]] = []
            all_mut_warnings: List[str] = []
            mutated_sequences: List[str] = []
            for mf, seq in zip(mut_fields, sequences):
                edits, mut_tags, mut_warns = parse_mutations(mf, seq)
                mutation_tags_per_chain.append(mut_tags)
                all_mut_warnings.extend([f'[{complex_id}] {w}' for w in mut_warns])
                mutated_sequences.append(apply_mutations(seq, edits))
            sequences = mutated_sequences

            modifications: List[List[dict]] = []
            all_warnings: List[str] = []
            for ptm_field, seq in zip(ptm_fields, sequences):
                mods, warns = parse_ptms(ptm_field, seq)
                modifications.append(mods)
                all_warnings.extend([f'[{complex_id}] {w}' for w in warns])

            pretty_name = format_complex_name(names, sequences, modifications, mutation_tags_per_chain)

            payload = build_alphafold_json(names, sequences, model_seeds, modifications, chain_ids)
            payload['name'] = pretty_name

            out_path = out_dir / f"{payload['name']}.json"
            with out_path.open("w", encoding="utf-8") as fh:
                json.dump(payload, fh, indent=2)

            for w in all_mut_warnings + all_warnings:
                print(f'WARNING: {w}', file=sys.stderr)

    else:
        triples = _collect_wide_columns(df)
        total_rows = len(df)
        print(f'Processing {total_rows} rows in wide format...')

        for row in tqdm(df.itertuples(index=False), total=total_rows, desc='Rows'):
            names: List[str] = []
            accs: List[str] = []
            ptm_fields: List[str] = []
            mut_fields: List[str] = []
            copies: List[int] = []

            for (name_col, acc_col, ptm_col) in triples:
                names.append(str(getattr(row, name_col)))
                accs.append(str(getattr(row, acc_col)))
                ptm_fields.append(str(getattr(row, ptm_col)) if ptm_col else '')

                idx_str = name_col.replace('name','')
                try:
                    mut_val = getattr(row, f'muts{idx_str}')
                except Exception:
                    mut_val = ''
                mut_fields.append(str(mut_val) if mut_val is not None else '')

                try:
                    cval = getattr(row, f'copies{idx_str}')
                    copies.append(int(cval))
                except Exception:
                    copies.append(1)

            sequences = _fetch_sequences(accs)
            chain_ids = _assign_chain_ids_by_copies(copies)

            mutation_tags_per_chain: List[List[str]] = []
            all_mut_warnings: List[str] = []
            mutated_sequences: List[str] = []
            for mf, seq in zip(mut_fields, sequences):
                edits, mut_tags, mut_warns = parse_mutations(mf, seq)
                mutation_tags_per_chain.append(mut_tags)
                all_mut_warnings.extend(mut_warns)
                mutated_sequences.append(apply_mutations(seq, edits))
            sequences = mutated_sequences

            modifications: List[List[dict]] = []
            all_warnings: List[str] = []
            for ptm_field, seq in zip(ptm_fields, sequences):
                mods, warns = parse_ptms(ptm_field, seq)
                modifications.append(mods)
                all_warnings.extend(warns)

            pretty_name = format_complex_name(names, sequences, modifications, mutation_tags_per_chain)

            payload = build_alphafold_json(names, sequences, model_seeds, modifications, chain_ids=chain_ids)
            payload['name'] = pretty_name

            out_path = out_dir / f"{payload['name']}.json"
            with out_path.open('w', encoding='utf-8') as fh:
                json.dump(payload, fh, indent=2)

            for w in all_mut_warnings + all_warnings:
                print(f'WARNING: {w}', file=sys.stderr)

    print('Done.')

if __name__ == '__main__':
    main()
