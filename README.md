# AF3 Prep — AlphaFold3 JSON Builder (PTMs • Mutations • Homomeric Copies)

Turn a simple table into **AlphaFold3 input JSON** files. The tool:

* fetches **UniProt** sequences,
* applies **mutations** to sequences,
* adds **PTMs** as AF3 `modifications`,

Scripts:

* `af3_prep.py` (CLI)
* `prep_utils.py` (helpers)


## 1) Input formats

The script auto-detects **long** (recommended) vs **wide** input.

### A) Long / normalized (recommended)

Each row = **one chain**.

**Required columns**

* `complex` — complex/group ID for one output JSON (e.g., `Clust1`)
* `name` — chain label to appear in the JSON `"name"`
* `uniprot` — UniProt accession (e.g., `Q9UMX1`)

**Optional columns**

* `ptms` — PTM tokens (see §3)
* `muts` — mutation tokens (see §4)
* `copies` — homomer copy count (default 1)

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

---

## 3) PTM syntax (column `ptms` / `ptmsN`)

Separate multiple tokens with `;`, `,`, or spaces.

**Accepted forms**

* **Explicit code @ position** (3-letter code):

  * `SEP@560`, `TPO@597`, `PTR@123`
* **Shorthand phosphorylation**:

  * `pS560`, `pT597`, `pY123` (also `S560p`, `T597p`, `Y123p`)
* **Forgiving**:

  * `S560-SEP`, `S560_SEP`

**Residue expectations (for warnings)**

* `SEP` → S
* `TPO` → T
* `PTR` → Y
* `ALY`, `MLY`, `M3L` → K
* `MSE` → M

> You can use any PDB 3-letter component for `ptmType`; the script only warns on suspicious residue/site mismatches.

---

## 4) Mutation syntax (column `muts` / `mutsN`)

Separate multiple tokens with `;`, `,`, or spaces. **1-based positions** (UniProt numbering).

**Accepted forms**

1. **Direct substitution**
   `A272G` or `Ala272Gly`
2. **Position → target residue**
   `272->A` or `272->Ala`
3. **Batch positions → target**
   `toA:560,594` or `->Gly:560-565`
4. **Source→target with positions**
   `A->G:272,300` or `Ala->Gly:15-20`
5. **Global replacement (danger!)**
   `A->G:*` (mutate **all** alanines to glycine)

**Behavior**

* Mutations are **applied to the sequence before PTM parsing**.
* If a mutation causes a PTM residue mismatch or a site is out-of-range, the script prints a **warning** and continues.

---

## 5) Homomeric copies & chain IDs

* Set `copies` (or `copiesN`) per chain; default `1`.
* IDs are auto-assigned **across the complex/row** as **A, B, …, Z, AA, AB, …**.
* AF3 JSON uses **list IDs** for homomers:

```json
"protein": {
  "id": ["A","B"],
  "sequence": "..."
}
```

---

## 6) Naming & filenames

**PTM naming (complex-level rule):**

* If **total PTM sites < 3** (i.e., **1 or 2**): add **exact site tags** to chains with PTMs, e.g. `SUFU(S560)` or `SUFU(S560,S594)`.
* If **total PTM sites ≥ 3**: chains with any PTM get **`(ptms)`**, e.g. `SUFU(ptms)_GLI3_GLl3CL(ptms)`.

**Mutation tags**

* If a chain has mutations, append `[mu/...]`.
* If a chain has >3 mutation tags, it collapses to `[muts]`.

**Filename rule**

* The **output filename is exactly equal to the JSON `"name"`**, with `.json` appended.
  Example: `"name": "SUFU(S560)[mut/A272G]_GLI3"` → file `SUFU(S560)[mut/A272G]_GLI3.json`.

---

## 7) How to run

**Long format (recommended)**

```bash
python af3_prep.py -i my_complexes.xlsx -o out -s 1-20
```

**Wide format**

```bash
python af3_prep.py -i my_wide.csv -o out -s 1-10
```

**Args**

* `-i / --input` — `.xlsx`, `.csv`, or `.tsv` (auto-detected)
* `-o / --output` — output folder (default: current directory)
* `-s / --seeds` — model seeds as comma/ranges (e.g., `1-4,6,8`)

**Output**

* One JSON per complex (long) or per row (wide), named by the final `"name"`.

---

## 8) Example output (snippet)

```json
{
  "name": "SUFU(S560,S594)[mut:A272G]_GLI3_GLl3CL(S52)",
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

---

## 9) Troubleshooting

* **Empty sequence / UniProt fetch failed**: check the accession and your network.
* **PTM residue mismatch**: e.g., `SEP@560` but residue at 560 isn’t `S`. Verify positions and mutations.
* **Position out of range**: site beyond sequence length; confirm your UniProt ID and numbering.
* **Many copies**: IDs continue past `Z` as `AA, AB, …`.

---

## 10) Tips

* Prefer **long format** for flexibility and QC.
* Use `muts` for **sequence edits**; `ptms` for **chemical modifications**.
* Keep tokens simple; whitespace and separators are forgiving.
* Consider pinning seeds (e.g., `-s 1-20`) to stabilize your batch workflow.
