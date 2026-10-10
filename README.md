# AF3 helper

Helper scripts around an AlphaFold3 run: build the input JSONs from a table, then
collect and plot the confidence metrics from the inference output.

| Folder | Scripts | Purpose |
| --- | --- | --- |
| `input_generation/` | `af3_prep.py`, `prep_utils.py` | Table → AlphaFold3 input JSON (UniProt sequences, PTMs, mutations, homomer copies) |
| `analysis/` | `af3_results_analyze.py`, `af3_utils.py` | Inference output → per-run metrics CSVs and per-dataset plots |
| `third_party/` | `ipsae.py` | ipSAE / pDockQ / pDockQ2 / LIS scoring, called by the analysis script |
| `tools/` | `seq_differ.py` | Position-by-position comparison of two sequences |

## Install

```bash
pip install -r requirements.txt
```

---

# Part 1: input generation

Turn a simple table into **AlphaFold3 input JSON** files. The tool:

* fetches **UniProt** sequences,
* applies **mutations** to sequences,
* adds **PTMs** as AF3 `modifications`.

## 1) Input formats

The script auto-detects **long** (recommended) vs **wide** input.

### A) Long / normalized (recommended)

Each row = **one chain**.

**Required columns**

* `complex`: complex/group ID for one output JSON (e.g., `Clust1`)
* `name`: chain label to appear in the JSON `"name"`
* `uniprot`: UniProt accession (e.g., `Q9UMX1`)

**Optional columns**

* `ptms`: PTM tokens (see §2)
* `muts`: mutation tokens (see §3)
* `copies`: homomer copy count (default 1)

**Example**

| complex | name     | copies | uniprot | ptms        | muts        |
| ------: | -------- | :----: | :-----: | ----------- | ----------- |
|  Clust1 | SUFU     |    2   |  Q9UMX1 | pS560;pS594 | A272G       |
|  Clust1 | GLI3     |    1   |  P10071 |             |             |
|  Clust1 | GLI3CL   |    1   |  P10071 | pS52        |             |
|  Clust2 | SMO      |    1   |  Q99835 | SEP\@560    | toA:560,594 |
|  Clust2 | NBEAcore |    1   |  Q8N163 |             | A->G:\*     |

### B) Wide format (2+ chains per row)

Numbered columns per chain:

* `name1, acc1 (or uniprot1), ptms1, muts1, copies1`
* `name2, acc2 (or uniprot2), ptms2, muts2, copies2`
* …and so on

## 2) PTM syntax (column `ptms` / `ptmsN`)

Separate multiple tokens with `;`, `,`, or spaces.

**Accepted forms**

* **Explicit code @ position** (3-letter code): `SEP@560`, `TPO@597`, `PTR@123`
* **Shorthand phosphorylation**: `pS560`, `pT597`, `pY123` (also `S560p`, `T597p`, `Y123p`)
* **Forgiving**: `S560-SEP`, `S560_SEP`

**Residue expectations (for warnings)**

* `SEP` → S
* `TPO` → T
* `PTR` → Y
* `ALY`, `MLY`, `M3L` → K
* `MSE` → M

> You can use any PDB 3-letter component for `ptmType`; the script only warns on suspicious residue/site mismatches.

## 3) Mutation syntax (column `muts` / `mutsN`)

Separate multiple tokens with `;`, `,`, or spaces. **1-based positions** (UniProt numbering).

**Accepted forms**

1. **Direct substitution**: `A272G` or `Ala272Gly`
2. **Position → target residue**: `272->A` or `272->Ala`
3. **Batch positions → target**: `toA:560,594` or `->Gly:560-565`
4. **Source→target with positions**: `A->G:272,300` or `Ala->Gly:15-20`
5. **Global replacement (danger!)**: `A->G:*` (mutate **all** alanines to glycine)

**Behavior**

* Mutations are **applied to the sequence before PTM parsing**.
* If a mutation causes a PTM residue mismatch or a site is out-of-range, the script prints a **warning** and continues.

## 4) Homomeric copies & chain IDs

* Set `copies` (or `copiesN`) per chain; default `1`.
* IDs are auto-assigned **across the complex/row** as **A, B, …, Z, AA, AB, …**.
* AF3 JSON uses **list IDs** for homomers:

```json
"protein": {
  "id": ["A","B"],
  "sequence": "..."
}
```

## 5) Naming & filenames

**PTM naming (complex-level rule)**

* If **total PTM sites < 3** (1 or 2): chains with PTMs get **exact site tags**, e.g. `SUFU(S560)` or `SUFU(S560,S594)`.
* If **total PTM sites ≥ 3**: chains with any PTM get **`(ptms)`**, e.g. `SUFU(ptms)_GLI3_GLI3CL(ptms)`.

**Mutation tags**

* A chain with up to 3 mutation tokens gets them appended as `[_...]`, e.g. `SUFU[_A272G]`.
* A chain with more than 3 collapses to `[_muts]`.

**Filename rule**

* The output filename is exactly the JSON `"name"` plus `.json`.
  Example: `"name": "SUFU(S560)[_A272G]_GLI3"` → file `SUFU(S560)[_A272G]_GLI3.json`.

## 6) How to run

```bash
cd input_generation
python af3_prep.py -i my_complexes.xlsx -o out -s 1-20
```

**Args**

* `-i / --input`: `.xlsx`, `.csv`, or `.tsv` (auto-detected)
* `-o / --output`: output folder (default: current directory)
* `-s / --seeds`: model seeds as comma/ranges (e.g., `1-4,6,8`; default `1-20`)

**Output**

* One JSON per complex (long) or per row (wide), named by the final `"name"`.

## 7) Example output (snippet)

The `Clust1` rows above have three PTM sites in total, so the name uses `(ptms)`:

```json
{
  "name": "SUFU(ptms)[_A272G]_GLI3_GLI3CL(ptms)",
  "sequences": [
    {
      "protein": {
        "id": ["A","B"],
        "sequence": "M....",
        "modifications": [
          {"ptmType": "SEP", "ptmPosition": 560},
          {"ptmType": "SEP", "ptmPosition": 594}
        ]
      }
    },
    {
      "protein": {"id": ["C"], "sequence": "M...."}
    },
    {
      "protein": {
        "id": ["D"],
        "sequence": "M....",
        "modifications": [{"ptmType": "SEP", "ptmPosition": 52}]
      }
    }
  ],
  "modelSeeds": [1,2,3,4,5],
  "dialect": "alphafold3",
  "version": 1
}
```

## 8) Troubleshooting

* **Empty sequence / UniProt fetch failed**: check the accession and your network.
* **PTM residue mismatch**: e.g., `SEP@560` but residue at 560 isn’t `S`. Verify positions and mutations.
* **Position out of range**: site beyond sequence length; confirm your UniProt ID and numbering.
* **Many copies**: IDs continue past `Z` as `AA, AB, …`.
* **File cannot be written on Windows**: mutation forms 2–5 put `:`, `>` or `*` into the name, which Windows does not allow in filenames. Use the direct form (`A272G`) or run on Linux.

---

# Part 2: results analysis

Walks an AlphaFold3 inference output folder, collects the confidence metrics of every
seed/sample, scores each model with `ipsae.py`, and plots the metrics per dataset.

## 1) Expected folder layout

```
base_dir/
├── Clust1/                         first level under base_dir = dataset label
│   └── .../seed-1_sample-0/
│       ├── *summary_confidences.json
│       ├── *confidences.json
│       └── *.cif
└── Clust2/
    └── ...
```

A run is any `*summary_confidences.json` whose folder name or filename contains
`seed-N ... sample-M`. The top-level "best model" summary has no seed/sample and is skipped.

## 2) How to run

```bash
cd analysis
python af3_results_analyze.py --base_dir /path/to/inference_output --out_dir results
```

**Args**

* `--base_dir`: root folder containing the dataset subfolders (required)
* `--out_dir`: where CSVs and plots are written (required)
* `--plot`: `dots` (default), `box`, or `violin`
* `--ipsae_py`: path to `ipsae.py` (default: the bundled `third_party/ipsae.py`); pass `''` to skip ipSAE scoring
* `--ipsae_pae`: ipSAE PAE cutoff (default 10)
* `--ipsae_dist`: ipSAE distance cutoff (default 10)

Each plot is also shown on screen. On a machine without a display, set `MPLBACKEND=Agg`.

## 3) Output

| File | Content |
| --- | --- |
| `all_interface_metrics.csv` | One row per seed/sample: `iptm`, `ptm`, `ranking_score`, chain-pair matrices, file paths, plus `iface_pTM_max` and `iface_minPAE` (best off-diagonal chain-pair ipTM and minimum PAE) |
| `ipsae_chainpair_scores.csv` | `ipsae.py` chain-pair scores for every run (ipSAE, pDockQ, pDockQ2, LIS, …) |
| `ipsae_residue_scores.csv` | `ipsae.py` per-residue scores for every run |
| `<plot>_iptm_by_dataset.png`, `<plot>_iface_minPAE_by_dataset.png` | Per-dataset distributions, ordered by median |
| `<plot>_ipSAE_by_dataset.png`, `..._pDockQ_...`, `..._pDockQ2_...`, `..._LIS_...` | Same, for the `Type=max` chain-pair rows |

`ipsae.py` also writes its own `<model>_<pae>_<dist>.txt`, `_byres.txt` and `.pml` files
next to each CIF.

---

# Sequence comparison

`tools/seq_differ.py` lists the positions where two sequences differ. Paste the two
sequences into the `seq1` / `seq2` variables at the bottom of the file and run it:

```bash
python tools/seq_differ.py
```

# Credits

`third_party/ipsae.py` is by Roland Dunbrack (Fox Chase Cancer Center), distributed under
the MIT license. If you use the ipSAE scores, cite
[Dunbrack 2025](https://www.biorxiv.org/content/10.1101/2025.02.10.637595v1).
