# Reproducing the published assemblies

The assemblies in the InstaNexus paper ([doi:10.1016/j.mcpro.2026.101547](https://doi.org/10.1016/j.mcpro.2026.101547))
are deposited on Zenodo ([10.5281/zenodo.18414969](https://doi.org/10.5281/zenodo.18414969)), together with the
predictions they were assembled from:

```
supplementary_data/
├── predictions/<SAMPLE>.csv              # e.g. BSA.csv, NB13.csv, MA1.csv
└── assembly/<key>/<config>/              # e.g. bsa/greedy_fdr0.1, ma1/dbg_weighted_heavy_fdr0.1_ks7
    ├── cleaned.csv
    └── scaffolds/
```

Rerunning them needs three steps that the deposit does not record. **Each one changes the result silently when it is
missed**: the run completes, but with different input or different scaffolds.

## 1. Convert the prediction files

The deposited prediction files cannot be read as they are:

- the peptide is stored only as a Python list of tokens in `prediction`, e.g. `"['[UNIMOD:385]', 'A', 'A', ...]"`,
  with no untokenised peptide column;
- there is no `experiment_name` column, which InstaNexus reads the protease from.

Read as deposited, InstaNexus keeps the brackets, quotes and commas as part of the peptide, and almost every PSM is
filtered out: 34 of BSA's 680 PSMs survive at `--fdr 0.1`.

Add the two columns before running InstaNexus. The peptide is the tokens joined, and the run is `spectrum_id` up to
its last `:`:

```python
import ast

import pandas as pd

df = pd.read_csv("supplementary_data/predictions/BSA.csv", index_col=0)
df["prediction_untokenised"] = df["prediction"].map(lambda s: "".join(ast.literal_eval(s)))
df["experiment_name"] = df["spectrum_id"].str.rsplit(":", n=1).str[0]
df.to_csv("inputs/bsa.csv", index=False)  # name the file after the metadata key, see step 3
```

## 2. Use `--refine` for the `dbg_weighted` assemblies

The folder names give the assembly mode, chain, FDR threshold and k-mer size (`dbg_weighted_heavy_fdr0.1_ks7`), but
not `--refine`. Every published `dbg_weighted` assembly used it. Without it, the scaffolds differ: for BSA,
`dbg_weighted` gives 33 scaffolds without `--refine` and 26 with it, and the 26 are identical to the published
`scaffolds.fasta`. The `greedy` assemblies do not use `--refine`.

## 3. Map the Zenodo sample names to the metadata keys

InstaNexus looks a sample up in `json/sample_metadata.json` by the stem of the input file name, and the
`assembly/` folders use those keys. The prediction files use different names:

| `predictions/` file | `json/sample_metadata.json` key and `assembly/` folder |
|---|---|
| `BSA` | `bsa` |
| `NB1`–`NB6` | `nb1`–`nb6` |
| `NB8` | `nb7` |
| `NB10` | `nb8` |
| `NB12` | `nb9` |
| `NB13` | `nb10` |
| `BIND15`, `BIND16`, `BIND17` | `bind1`, `bind2`, `bind3` |
| `MA1`, `MA2`, `MA3` | `ma1`, `ma2`, `ma3` |
| `PA` (the oligoclonal mix) | no single key (`ma1` + `ma2` + `ma3`) |

In particular, `predictions/NB10.csv` is `nb8`, and sample Nb10 of the published assemblies corresponds to
`predictions/NB13.csv` in the deposit.

## Example

BSA, `dbg_weighted` at FDR 0.1 with k = 7, after converting `BSA.csv` to `inputs/bsa.csv` (step 1):

```bash
instanexus \
    --input-csv inputs/bsa.csv \
    --folder-outputs outputs \
    --metadata-json-path json/sample_metadata.json \
    --contaminants-fasta-path fasta/contaminants.fasta \
    --reference \
    --assembly-mode dbg_weighted \
    --fdr 0.1 \
    --kmer-size 7 \
    --refine
```

For the antibodies, add `--chain light` or `--chain heavy`, as in the folder name.

With these three steps:

- the filtered PSMs in `cleaned.csv` are identical to the published ones (checked for BSA, mAb1 and the samples
  deposited as `NB6`, `NB8`, `NB10` and `NB13`);
- assembling each published `cleaned.csv` gives a `scaffolds.fasta` identical to the published one, for all 40
  configurations: the 20 `dbg_weighted` assemblies with `--refine` and the 20 `greedy` assemblies without it
  (`--min-overlap 3 --size-threshold 10`, the defaults).

`MA1_high.csv` and `MA1_low.csv` are in InstaNovo's own output format (`preds`, `preds_tokenised`), which InstaNexus
reads as is. They have no published assembly folder.

## Expected differences

- **Assembly statistics.** Since [#61](https://github.com/Multiomics-Analytics-Group/InstaNexus/issues/61) the
  reference statistics use the mapping coordinates correctly, so `coverage` is lower and lengths, N50 and N90 are
  usually 1 shorter than in the deposited `statistics/*.json`. See the CHANGELOG.
- **Isoleucine.** The published scaffolds report every I/L position as L. Versions that restore isoleucine
  ([#66](https://github.com/Multiomics-Analytics-Group/InstaNexus/pull/66)) skip the restoration when the reads
  contain only one of the two letters. Every deposited prediction file contains I and no L, so the restoration is
  skipped and the scaffolds are unchanged; `--no-isoleucine-restoration` makes that explicit.
