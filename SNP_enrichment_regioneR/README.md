# oku_taiga_paper_2026_1

Analysis scripts for the manuscript

> **Re-detection of previously reported athletic-performance genes by integrating GWAS-associated SNPs
> across related traits**

The study collects GWAS-associated SNPs for 16 physiological traits related to power and endurance,
and asks — gene by gene, trait by trait — whether those SNPs are positionally enriched in gene regions.
The main outcome is how many genes previously reported by candidate-gene studies are re-detected by
this framework, and through which traits.

## Paper information

- Authors / journal / DOI: *to be added on acceptance*
- Correspondence: see the manuscript

## Directory structure

```
.
├── SNP_enrichment_regioneR/   core analysis: positional enrichment of SNPs in gene regions (regioneR)
├── LICENSE
└── README.md
```

`SNP_enrichment_regioneR/` contains the two analysis scripts, the gene-region BED files, example input
for a trial run, and records of the software environment. See the README inside that directory for the
command-line interface and the output format.

### Scope of this repository

This repository publishes the **core analysis**: the positional enrichment test that produces the gene
lists on which every downstream result rests.

Upstream preparation (collecting SNPs from the GWAS Catalog, mapping them to genomic positions against
1000 Genomes, LD pruning with PLINK2) and downstream steps (re-detection tables, GO enrichment with
Metascape, figures) are described in full in the Methods section of the manuscript, with the parameters
needed to reproduce them. They are not included here because they depend on large external datasets
rather than on code specific to this study.

> **Note on an earlier version of this repository.** A previous release also contained
> `Calc_Zscore_regioneReloaded/`, scripts for a local Z-score analysis with `regioneReloaded`.
> That analysis was removed during revision and is **not part of the current manuscript**, so the
> directory has been deleted.

## Environment

The environment used for the analysis is recorded in `SNP_enrichment_regioneR/`
(`ENV_sessionInfo.txt`, `ENV_bioconductor_version.txt`, `ENV_installed_packages.csv`) and can be
restored with the `renv.lock` file in the same directory. The scripts require macOS or Linux because
they parallelise with `parallel::mclapply()`.

## Trial run

A small example is included so that the scripts can be run without preparing any data:

```bash
cd SNP_enrichment_regioneR
Rscript GB_regioneR.R Endurance MUO2 0004887 0.2 EUR up0_down0 2 \
        ./out ./example_inputSNP_Endurance_EUR_LD_r0.2_MUO2_0004887.tsv 0.05
```

See `SNP_enrichment_regioneR/README.md` for the full argument list and for the flanking-region variant.

## Contact and bug reports

Please open an issue on GitHub. The code is provided for research purposes and is supported on a
best-effort basis.

## License

MIT License (see `LICENSE`).
