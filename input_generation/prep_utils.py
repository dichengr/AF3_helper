import re
import time
from typing import Optional, List, Dict, Tuple, Union
import requests

# --------------------------------------------------------------------------------------
# UniProt
# --------------------------------------------------------------------------------------

def fetch_uniprot_sequence(
    accession: str,
    retries: int = 3,
    backoff: int = 3,
    timeout: int = 30,
) -> Optional[str]:
    """
    Fetch the amino-acid sequence for a UniProt protein accession.

    Parameters
    ----------
    accession : str
        UniProt accession (e.g. "P68871").
    retries : int, default 3
        Number of retry attempts in case of server errors or rate limits.
    backoff : int, default 3
        Seconds to wait between retries (exponential back-off).
    timeout : int, default 30
        Seconds before giving up on a single request.

    Returns
    -------
    str | None
        The raw sequence (no FASTA header) or None if not found.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"

    for attempt in range(retries):
        try:
            response = requests.get(url, timeout=timeout)

            # Success
            if response.status_code == 200:
                fasta_text = response.text
                sequence = "".join(
                    line.strip() for line in fasta_text.splitlines() if not line.startswith(">")
                )
                return sequence

            # Not found – stop early
            if response.status_code == 404:
                return None

            # Rate-limited or server error – wait and retry
            time.sleep(backoff * (attempt + 1))

        except requests.RequestException:
            time.sleep(backoff * (attempt + 1))

    return None


# --------------------------------------------------------------------------------------
# PTM parsing & validation
# --------------------------------------------------------------------------------------

# Simple expectation map for common residue replacements used as PTMs (PDB-style 3-letter codes)
_EXPECTED_RESIDUE: Dict[str, str] = {
    # phosphorylation
    "SEP": "S",  # phospho-Ser
    "TPO": "T",  # phospho-Thr
    "PTR": "Y",  # phospho-Tyr
    # lysine acyl/methyl variants (best-effort; user can still supply any 3-letter code)
    "ALY": "K",  # N6-acetyl-L-lysine
    "MLY": "K",  # N6-methyl-L-lysine
    "M3L": "K",  # N6,N6,N6-trimethyl-L-lysine
    # others
    "MSE": "M",  # selenomethionine
}

def _tokenize_field(field: Optional[str]) -> List[str]:
    """
    Split a free-text field into atomic tokens using commas, semicolons, or whitespace.
    """
    if field is None:
        return []
    s = str(field).strip()
    if not s or s.lower() in {"nan", "none", "-", ""}:
        return []
    return [t for t in re.split(r"[;,\s]+", s) if t]


def parse_ptms(ptm_field: Optional[str], sequence: str) -> Tuple[List[Dict[str, Union[int, str]]], List[str]]:
    """
    Parse a PTM field into AlphaFold3 'modifications' entries.

    Supported syntaxes (case-insensitive):
      • Explicit code@pos:   SEP@560, TPO@597, PTR@123
      • Shorthand phospho:   pS560, pT597, pY123
      • Alternative order:   S560p, T597p, Y123p

    Returns
    -------
    (mods, warnings)
      mods:      list of {"ptmType": <3-letter>, "ptmPosition": <1-based int>}
      warnings:  human-readable validation warnings (residue mismatches / out-of-range)
    """
    tokens = _tokenize_field(ptm_field)
    mods: List[Dict[str, Union[int, str]]] = []
    warnings: List[str] = []

    for tok in tokens:
        t = tok.strip()

        m = re.fullmatch(r"([A-Za-z]{3})@(\d+)", t)
        if m:
            code = m.group(1).upper()
            pos = int(m.group(2))
            mods.append({"ptmType": code, "ptmPosition": pos})
            continue

        m = re.fullmatch(r"p([STY])(\d+)", t, flags=re.IGNORECASE)
        if m:
            aa = m.group(1).upper()
            pos = int(m.group(2))
            code = {"S": "SEP", "T": "TPO", "Y": "PTR"}[aa]
            mods.append({"ptmType": code, "ptmPosition": pos})
            continue

        m = re.fullmatch(r"([STY])(\d+)p", t, flags=re.IGNORECASE)
        if m:
            aa = m.group(1).upper()
            pos = int(m.group(2))
            code = {"S": "SEP", "T": "TPO", "Y": "PTR"}[aa]
            mods.append({"ptmType": code, "ptmPosition": pos})
            continue

        # Forgiving variant: S560-SEP or S560_SEP
        m = re.fullmatch(r"([A-Za-z])(\d+)[\-\_]([A-Za-z]{3})", t)
        if m:
            pos = int(m.group(2))
            code = m.group(3).upper()
            mods.append({"ptmType": code, "ptmPosition": pos})
            continue

        warnings.append(f"Unrecognized PTM token '{tok}' – skipped. Use e.g. 'SEP@560' or 'pS560'.")

    # Validate residue matches and positions
    for mod in mods:
        pos = int(mod["ptmPosition"])
        code = str(mod["ptmType"]).upper()
        if pos < 1 or pos > len(sequence):
            warnings.append(f"PTM {code}@{pos} is out of range (sequence length {len(sequence)}).")
            continue
        expected = _EXPECTED_RESIDUE.get(code)
        if expected:
            actual = sequence[pos - 1]
            if actual != expected:
                warnings.append(
                    f"PTM {code}@{pos} expects residue '{expected}', but sequence has '{actual}'."
                )

    return mods, warnings


def _site_tags_for_chain(mods: List[Dict[str, Union[int, str]]], sequence: str) -> List[str]:
    """Return compact tags like 'S560' that match the base residue in the original sequence."""
    tags: List[str] = []
    for m in mods or []:
        pos = int(m["ptmPosition"])
        if 1 <= pos <= len(sequence):
            res = sequence[pos - 1]
            tags.append(f"{res}{pos}")
        else:
            tags.append(str(pos))
    return tags


# --------------------------------------------------------------------------------------
# Mutations parsing & application
# --------------------------------------------------------------------------------------

_THREE_TO_ONE: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
_ONE_TO_THREE: Dict[str, str] = {v: k for k, v in _THREE_TO_ONE.items()}

def _aa_to_one(aa: str) -> Optional[str]:
    aa = aa.strip().upper()
    if len(aa) == 1:
        return aa
    return _THREE_TO_ONE.get(aa)

def _parse_pos_list(s: str) -> List[int]:
    # supports "1,2,5-9"
    out: List[int] = []
    for chunk in re.split(r"[,\s]+", s.strip()):
        if not chunk:
            continue
        if "-" in chunk:
            a, b = chunk.split("-", 1)
            try:
                out.extend(list(range(int(a), int(b) + 1)))
            except ValueError:
                continue
        else:
            try:
                out.append(int(chunk))
            except ValueError:
                continue
    return out

def parse_mutations(muts_field: Optional[str], sequence: str) -> Tuple[List[Tuple[int, str, Optional[str]]], List[str], List[str]]:
    """
    Parse mutation expressions into edits on the sequence.

    Returns
    -------
    (edits, tags, warnings)
      edits:    list of (position 1-based, new_aa ONE-LETTER, expected_old ONE-LETTER or None)
      tags:     compact tags for naming (e.g., ['A272G','toA:560','A->G*'])
      warnings: human-readable validation warnings
    """
    tokens = _tokenize_field(muts_field)
    edits: List[Tuple[int, str, Optional[str]]] = []
    tags: List[str] = []
    warnings: List[str] = []

    if not tokens:
        return edits, tags, warnings

    for tok in tokens:
        t = tok.strip()

        # Pattern 1: direct substitution like A272G or Ala272Gly
        m = re.fullmatch(r"([A-Za-z]{1,3})(\d+)([A-Za-z]{1,3})", t)
        if m:
            old = _aa_to_one(m.group(1))
            pos = int(m.group(2))
            new = _aa_to_one(m.group(3))
            if old and new:
                edits.append((pos, new, old))
                tags.append(f"{old}{pos}{new}")
                continue

        # Pattern 2: position to target: 272->A or 272->Ala
        m = re.fullmatch(r"(\d+)\s*->\s*([A-Za-z]{1,3})", t)
        if m:
            pos = int(m.group(1))
            new = _aa_to_one(m.group(2))
            if new:
                edits.append((pos, new, None))
                tags.append(f"{pos}->{new}")
                continue

        # Pattern 3: batch to target: toA:560,594 or ->Gly:560-565
        m = re.fullmatch(r"(?:to|->)\s*([A-Za-z]{1,3})\s*:\s*([0-9,\s\-]+)", t, flags=re.IGNORECASE)
        if m:
            new = _aa_to_one(m.group(1))
            pos_list = _parse_pos_list(m.group(2))
            if new and pos_list:
                for pos in pos_list:
                    edits.append((pos, new, None))
                tags.append(f"to{new}:{','.join(str(p) for p in pos_list[:3])}{'...' if len(pos_list)>3 else ''}")
                continue

        # Pattern 4: source->target with list: A->G:272,300 or Ala->Gly:15-20
        m = re.fullmatch(r"([A-Za-z]{1,3})\s*->\s*([A-Za-z]{1,3})\s*:\s*([0-9,\s\-]+)", t)
        if m:
            old = _aa_to_one(m.group(1))
            new = _aa_to_one(m.group(2))
            pos_list = _parse_pos_list(m.group(3))
            if old and new and pos_list:
                for pos in pos_list:
                    edits.append((pos, new, old))
                tags.append(f"{old}->{new}:{','.join(str(p) for p in pos_list[:3])}{'...' if len(pos_list)>3 else ''}")
                continue

        # Pattern 5: global replacement: A->G:*
        m = re.fullmatch(r"([A-Za-z]{1,3})\s*->\s*([A-Za-z]{1,3})\s*:\s*\*", t)
        if m:
            old = _aa_to_one(m.group(1))
            new = _aa_to_one(m.group(2))
            if old and new:
                for i, ch in enumerate(sequence, start=1):
                    if ch == old:
                        edits.append((i, new, old))
                tags.append(f"{old}->{new}*")
                continue

        warnings.append(f"Unrecognized mutation token '{tok}' – skipped.")

    # Validate positions and expected residue matches
    for (pos, new, old) in edits:
        if pos < 1 or pos > len(sequence):
            warnings.append(f"Mutation {old or '?'}{pos}{new} is out of range (len={len(sequence)}).")
            continue
        if old is not None:
            actual = sequence[pos - 1]
            if actual != old:
                warnings.append(f"Mutation {old}{pos}{new} expects '{old}' but sequence has '{actual}'. Applying anyway.")

    return edits, tags, warnings


def apply_mutations(sequence: str, edits: List[Tuple[int, str, Optional[str]]]) -> str:
    if not edits:
        return sequence
    seq_list = list(sequence)
    for (pos, new, _old) in edits:
        if 1 <= pos <= len(seq_list):
            seq_list[pos - 1] = new
    return "".join(seq_list)


def _mut_site_tags_for_chain(tags: List[str]) -> List[str]:
    return tags or []


# --------------------------------------------------------------------------------------
# Naming and JSON assembly
# --------------------------------------------------------------------------------------

def format_complex_name(
    names: List[str],
    sequences: List[str],
    modifications: List[List[Dict[str, Union[int, str]]]],
    mutation_tags: Optional[List[List[str]]] = None,
) -> str:
    """
    Apply the user rule for naming with PTMs, then append concise mutation tags.

    PTM Rules (interpreted at complex level):
      • If total PTM sites across the complex == 1  → keep names unchanged.
      • If total PTM sites < 3 (i.e., 2)            → add exact site tags to chains with PTMs.
      • If total PTM sites >= 3                     → add literal '(ptms)' to chains with ≥1 PTM.

    Mutations:
      • If a chain has mutation tags, append [mut:...].
        If tags per chain > 3, collapse to [mut].
    """
    total_sites = sum(len(m or []) for m in (modifications or []))
    base_parts: List[str] = []
    for name, seq, mods in zip(names, sequences, modifications or [[]] * len(names)):
        # PTM naming
        if total_sites <= 0:
            part = name
        elif total_sites < 3:
            if mods:
                site_tags = _site_tags_for_chain(mods, seq)
                part = f"{name}({','.join(site_tags)})"
            else:
                part = name
        else:
            part = f"{name}(ptms)" if mods else name
        base_parts.append(part)

    # Append mutation tags per-chain (if provided)
    if mutation_tags:
        with_mut: List[str] = []
        for part, tags in zip(base_parts, mutation_tags):
            if tags:
                if len(tags) <= 3:
                    with_mut.append(f"{part}[_{','.join(tags)}]")
                else:
                    with_mut.append(f"{part}[_muts]")
            else:
                with_mut.append(part)
        base_parts = with_mut

    return "_".join(base_parts)


def build_alphafold_json(
    names: List[str],
    sequences: List[str],
    model_seeds: List[int],
    modifications: Optional[List[List[Dict[str, Union[int, str]]]]] = None,
    chain_ids: Optional[List[Union[str, List[str]]]] = None,
) -> Dict[str, object]:
    """
    Assemble an AlphaFold3 input JSON payload.

    Parameters
    ----------
    names : list[str]
        Protein names in the same order as *sequences*.
    sequences : list[str]
        Corresponding raw amino acid sequences.
    model_seeds : list[int]
        Model seed integers as required by AlphaFold3.
    modifications : list[list[dict]], optional
        For each sequence, a list of {"ptmType": CODE, "ptmPosition": POS}.
    chain_ids : list[Union[str, list[str]]], optional
        If provided, element i may be a string ID or a list of IDs (homomeric copies).
    """
    if len(names) != len(sequences):
        raise ValueError("names and sequences must be the same length")
    if modifications is not None and len(modifications) != len(sequences):
        raise ValueError("modifications must match the number of sequences")
    if chain_ids is not None and len(chain_ids) != len(sequences):
        raise ValueError("chain_ids must match the number of sequences")

    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    seq_entries: List[Dict[str, object]] = []
    for i, seq in enumerate(sequences):
        if chain_ids:
            ids = chain_ids[i]
            if isinstance(ids, list):
                pid = [str(x) for x in ids]
            else:
                pid = [str(ids)]
        else:
            pid = [letters[i % len(letters)]]
        entry: Dict[str, object] = {"protein": {"id": pid, "sequence": seq}}
        if modifications and modifications[i]:
            entry["protein"]["modifications"] = modifications[i]
        seq_entries.append(entry)

    return {
        "name": "_".join(names),
        "sequences": seq_entries,
        "modelSeeds": model_seeds,
        "dialect": "alphafold3",
        "version": 1,
    }