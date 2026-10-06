# SNP enrichment in gene regions using regioneR

R scripts and example data for the positional enrichment analysis reported in the accompanying
manuscript. Given a set of trait-associated SNPs, the scripts test — gene by gene — whether the SNPs
overlap a gene region more often than expected by chance, using circular permutation (`regioneR`).

---

## Files

### Analysis scripts

| File | Gene region tested |
|---|---|
| `GB_regioneR.R` | **Gene body (GB).** SNPs inside the gene body are counted. |
| `extract_outGB_regioneR.R` | **Flanking only (outGB).** The gene body is subtracted from the expanded region, so only SNPs in the flanking sequence are counted. |

Notes

- The two scripts share the same algorithm, environment and command-line interface; they differ only in
  how the tested region is defined.
- For `extract_outGB_regioneR.R`, the subtraction of the gene body is applied **to the observed SNPs and
  to every permuted replicate**, so a flanking-only definition tests the flanking sequence alone.
- Both are intended to be run from the command line.
- With expanded regions, a loaded region may extend past the end of a chromosome and raise a warning.
  The script trims the overhang, so this is not an error.

### Gene region files (BED)

One BED per region definition. Column 4 is `ENSEMBL_ID|GENE_SYMBOL`, which the scripts use to label the
output (no separate ID-to-symbol mapping file is needed).

| File | Region |
|---|---|
| `protein_genes_up0_down0.bed` | Gene body (`bed_tag = up0_down0`) |
| `protein_genes_up5_down1p5.bed` | 5 kb upstream / 1.5 kb downstream (`up5_down1p5`) |
| `protein_genes_up10_down10.bed` | 10 / 10 kb (`up10_down10`) |
| `protein_genes_up20_down20.bed` | 20 / 20 kb (`up20_down20`) |
| `protein_genes_up35_down10.bed` | 35 / 10 kb (`up35_down10`) |
| `protein_genes_up50_down50.bed` | 50 / 50 kb (`up50_down50`) |

`extract_outGB_regioneR.R` always needs `protein_genes_up0_down0.bed` in addition to the expanded BED,
because it subtracts the gene body.

### Example input (for a trial run)

| File | Trait | SNPs |
|---|---|---|
| `example_inputSNP_Endurance_EUR_LD_r0.2_MUO2_0004887.tsv` | MUO2 (maximal oxygen uptake) | 10 |
| `example_inputSNP_Power_EUR_LD_r0.2_GPSM_0006941.tsv` | GPSM (grip strength) | 203 |

Both are at the baseline condition of the manuscript (LD pruning in EUR at r² = 0.2, 50 kb window).
The SNP file is a **headerless, tab-separated** file with three columns:

```
chromosome<TAB>position<TAB>SNP_ID
1	1430364	rs880315
```

- `position` is the 1-based VCF `POS`.
- **Each SNP must carry a unique ID**; duplicated IDs cause an error.
- Only autosomes are used; other chromosomes are dropped inside the script.

### Environment records

| File | Content |
|---|---|
| `ENV_sessionInfo.txt` | R version, platform and the attached packages |
| `ENV_bioconductor_version.txt` | Bioconductor release |
| `ENV_installed_packages.csv` | All installed packages and versions |
| `renv.lock` | `renv` lockfile for restoring the package environment |

---

## Requirements

- R (see `ENV_sessionInfo.txt` for the version used)
- **macOS or Linux** — the scripts use `parallel::mclapply()`, which does not fork on Windows
- Bioconductor packages: `regioneR`, `GenomicRanges`, `rtracklayer`, `GenomeInfoDb`, `S4Vectors`,
  `IRanges`, `BSgenome.Hsapiens.UCSC.hg19`, `BSgenome.Hsapiens.UCSC.hg19.masked`

The analysis was run on a MacBook Pro (M1 Max, 64 GB). Runtime scales with the number of permutations,
the number of input SNPs and the number of cores. Measured on two cores: about 30 seconds per run at
`K = 200`, and on the order of hours at the default `K = 50000`. See "A faster trial run" below.

## Setup

Restore the package environment with **renv**:

```bash
Rscript -e 'install.packages("renv", repos="https://cloud.r-project.org"); renv::restore(prompt = FALSE)'
```

`renv/library` is not included in this repository. To set the environment up by hand, use the versions
recorded in `ENV_installed_packages.csv`.

`renv.lock` and the `ENV_*` records describe the same environment (R 4.5.1, Bioconductor 3.21,
`regioneR` 1.40.1, `rtracklayer` 1.68.0, `GenomicRanges` 1.60.0,
`BSgenome.Hsapiens.UCSC.hg19.masked` 1.3.993).

---

## Running

Both scripts take the same arguments.

```
Rscript <script> <trait_category> <trait_name> <trait_id> <ld_r2> <ld_pop> <bed_tag> <cores> \
        [out_dir] [input_tsv] [sig_fdr] [bed_dir] [K]
```

### Trial run with the example files

Run these from inside this directory; the BED files are found automatically.

```bash
# Gene body
Rscript GB_regioneR.R Endurance MUO2 0004887 0.2 EUR up0_down0 2 \
        ./out ./example_inputSNP_Endurance_EUR_LD_r0.2_MUO2_0004887.tsv 0.05

# Flanking only (35 kb upstream / 10 kb downstream)
Rscript extract_outGB_regioneR.R Power GPSM 0006941 0.2 EUR up35_down10 2 \
        ./out ./example_inputSNP_Power_EUR_LD_r0.2_GPSM_0006941.tsv 0.05
```

### Arguments

| # | Argument | Meaning |
|---|---|---|
| 1 | `trait_category` | Trait category. **Used only to build output file names.** |
| 2 | `trait_name` | Trait abbreviation. Used in output file names. |
| 3 | `trait_id` | Trait ID (EFO). Used in output file names. |
| 4 | `ld_r2` | r² threshold used for LD pruning of the input. **Recorded, not applied here.** |
| 5 | `ld_pop` | Reference population used for LD pruning. **Recorded, not applied here.** |
| 6 | `bed_tag` | Which BED to use, e.g. `up0_down0`, `up35_down10`. The script reads `protein_genes_<bed_tag>.bed`. |
| 7 | `cores` | Number of cores for `mclapply()`. |
| 8 | `out_dir` | Output directory (default: current directory). Created if absent. |
| 9 | `input_tsv` | Path to the SNP TSV. If omitted, the script looks for `vcf_<ld_pop>_LD_r<ld_r2>_<trait_name>_<trait_id>_withID.tsv` under `$SNP_TSV_DIR` (default: current directory). |
| 10 | `sig_fdr` | FDR threshold for the "significant" output (default `0.05`). |
| 11 | `bed_dir` | Directory holding the BED files. If omitted, `$GENE_BED_DIR` is used; if that is unset, **the directory containing the script** — so the example below works with no configuration. |
| 12 | `K` | Number of permutations (default `50000`, the value used in the manuscript). **Lower it only for a quick trial**; see below. |

> Arguments 1–5 do not change the analysis; they are metadata used to name the inputs and outputs.
> If your files follow a different naming convention, pass the path explicitly with argument 9.

### A faster trial run

With the default `K = 50000` a single run takes on the order of hours. To check that the scripts work,
lower `K` with argument 12:

```bash
# finishes in about 30 seconds on two cores
Rscript GB_regioneR.R Endurance MUO2 0004887 0.2 EUR up0_down0 2 \
        ./out ./example_inputSNP_Endurance_EUR_LD_r0.2_MUO2_0004887.tsv 0.05 "" 200
```

The script prints a notice whenever `K` differs from 50000, and **the value of K appears in the output
file name**, so a trial result cannot be mistaken for a result reported in the manuscript.
Use `K = 50000` for any analysis you intend to interpret.

### Fixed analysis settings

These are set inside the script rather than exposed as arguments, because they define the test:

| Setting | Value |
|---|---|
| Permutations | `K = 50000` (argument 12 overrides it, for trial runs only) |
| Randomisation | `circularRandomizeRegions`, per chromosome |
| Genome / mask | `BSgenome.Hsapiens.UCSC.hg19.masked` |
| Chromosomes | autosomes (chr1–chr22) only |
| Multiple testing | Benjamini–Hochberg across the tested genes |
| Significance | `FDR_BH ≤ sig_fdr` **and** at least one observed overlap |

---

## Output

Two CSV files are written to `out_dir`:

| Script | File |
|---|---|
| `GB_regioneR.R` | `all_gene_<bed_tag>_result_rep_<K>_<pop>_LD_<r2>_<trait>_<id>.csv`<br>`significant_gene_<bed_tag>_result_rep_<K>_<pop>_LD_<r2>_<trait>_<id>.csv` |
| `extract_outGB_regioneR.R` | `all_gene_exp_outGB_<bed_tag>_result_rep_<K>_<pop>_LD_<r2>_<trait>_<id>.csv`<br>`significant_gene_exp_outGB_<bed_tag>_result_rep_<K>_<pop>_LD_<r2>_<trait>_<id>.csv` |

The `all_gene_` file lists every tested gene; the `significant_gene_` file lists only those passing the
FDR threshold.

Main columns: gene ID and symbol, `obs` (observed overlaps), the permutation mean and SD, `z`,
`p_emp` (empirical p) and `FDR_BH`.
