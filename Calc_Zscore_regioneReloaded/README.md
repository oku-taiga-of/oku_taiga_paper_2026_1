# Analysis script for local Z score in significant genes based on positional information using regioneReloaded
This directory provides R scripts and example data used for calculationg local Z score per significant genes by using regioneReloaded in the accompanying manuscript.

- **Version**: v1.0  
- **Last updated**: 2026-01-19

This script used **regioneReloaded**, the relevant Bioconductor page and papers are listed below.

[Bioconductor page](https://www.bioconductor.org/packages/release/bioc/html/regioneReloaded.html)

[Roberto Malinverni, David Corujo, Bernat Gel, Marcus Buschbeck, regioneReloaded: evaluating the association of multiple genomic region sets, Bioinformatics, Volume 39, Issue 11, November 2023, btad704, https://doi.org/10.1093/bioinformatics/btad704](https://academic.oup.com/bioinformatics/article/39/11/btad704/7439591)

We thank the developers of **regioneReloaded** for providing and maintaining a robust and well-documented framework for genomic region enrichment analysis.

## Files
### Analysis R scripts
`GB_calcZscore_regioneReloaded.R`: When referring to the gene body as the gene region, use this one.

`outGB_calcZscore_regioneReloaded.R`: When the gene region is the expansion region, use this one.

Note:
- The two scripts differ only in analysis mode; they share the same software environment.
- These scripts are intended to be executed via the CLI.
- The analysis algorithms for the two scripts are identical, but when using the extended region, overlaps are counted as the difference between the extended region and the gene body region (see the manuscript's method section for details).
- Note that the command arguments at runtime differ between the two scripts (the actual commands are described later).
- When using extended regions, the loaded extended region information may exceed the endpoints of the chromosome in the referenced genome, potentially causing a warning. Within the script, the exceeding portion is trimmed, so it does not directly cause a problem.


### Example input file
`example_inputSNP_Endurance_nonAFR_LD_r2_MUO2_0004887.tsv` : Inputting SNP data

`example_inputSNP_Power_nonAFR_LD_r2_GPSM_0006941.tsv` : Inputting SNP data

`gene_id_name_map.csv` : File for referencing gene symbol names using gene IDs as keys when outputting results

`gene_region_ORF.bed` : `.bed` Files for each gene region (gene body data)

`gene_region_u5.0k_d1.5k.bed` : `.bed` Files for each gene region (data for extended regions), File for the region expanded 5.0k bp upstream and 1.5k bp downstream.

`significant_gene_GB_result_rep_50000_nonAFR_LD_r2_MUO2_0004887_.csv` : Data of significant genes (output file from SNP_enrichment_regioneR)

`significant_gene_outGB_u5.0k_d1.5k_result_rep_50000_nonAFR_LD_r2_GPSM_0006941_.csv` : Data of significant genes (output file from SNP_enrichment_regioneR)

Note:
- These are example input files, so they are the minimum files required for demonstration and reproducibility checks.
- Since the script references columns, when running it in a different file, ensure the file structure matches the one containing the columns.
- Regarding SNP data, to prevent errors, each SNP must be assigned a unique SNP ID.

### Environment records
`renv.lock`: R package environment lockfile

`ENV_sessionInfo.txt`: Information about the R 

`ENV_bioconductor_version.txt`: Version of bioconductor

`ENV_installed_packages.csv`: Information about installed packages


## Requirements
- R version: see `ENV_sessionInfo.txt`
- macOS or Linux (the scripts use `parallel::mclapply()`)
- The analysis was performed using a MacBook Pro (M1 Max) with 64GB of memory. Depending on the number of cores used and the SNP data, the analysis completed within a few minutes to a few hours.

## Setup
**Assuming R is already installed**
Restore the R package environment using **renv**:

1) Ensure `renv` auto-activation is safe (recommended)

If `.Rprofile` contains `source("renv/activate.R")`, R may fail to start in a fresh directory
where `renv/activate.R` does not exist yet. Replace it with a guarded version:
```bash
printf 'if (file.exists("renv/activate.R")) source("renv/activate.R")\n' > .Rprofile
```
2) Initialize the renv project
```bash
Rscript -e 'install.packages("renv", repos="https://cloud.r-project.org"); renv::init(bare = TRUE)'
```
3) Restore the exact package environment from renv.lock
```bash
Rscript -e 'renv::restore(prompt = FALSE)'
```
4) Verify installation

Confirm that R is using the project library:
```bash
Rscript -e 'print(.libPaths())'
```
Check that the project is synchronized with the lockfile:
```bash
Rscript -e 'renv::status()'
```
At this point, the environment setup is complete. The analysis script will now run.

When published on GitHub, `renv/library` is not included.
* When manually setting up the environment, refer to the information contained in the Environment records.

## Run CLI commands (example)
### When the gene body is defined as the gene region
```bash
Rscript GB_calcZscore_regioneReloaded.R <trait_category> <trait_name> <trait_ID> <LD_r> <LD_pop> <bed_tag> <CORES>
```
Command examples corresponding to the sample files
```bash
Rscript GB_calcZscore_regioneReloaded.R Endurance MUO2 0004887 r2 nonAFR ORF 4
```
### When the expansion region is defined as the gene region
```bash
Rscript outGB_calcZscore_regioneReloaded.R <trait_category> <trait_name> <trait_ID> <LD_r> <LD_pop> <bed_tag_upstream> <bed_tag_downstream> <CORES>
```
Command examples corresponding to the sample files
```bash
Rscript outGB_calcZscore_regioneReloaded.R Power GPSM 0006941 r2 nonAFR 5.0 1.5 4
```
Both analysis scripts will return an error message and terminate without generating an output file if no analyzable genes are found.

### The command arguments
- `<trait_category>, <trait_name>, <trait_ID>, <LD_r>, <LD_pop>` : For each SNP dataset loaded: - Trait category - Trait abbreviation - Trait ID - r^2 threshold for LD pruning - Population used for LD pruning in the genomic dataset. Arguments for creating paths for SNP data to be loaded and output files.
- `<bed_tag>` : An argument used only when using the gene body, to be declared explicitly in the command.
- `<bed_tag_upstream>` <bed_tag_downstream> : Arguments used only when using the expansion region, used to create paths for the expansion　region file to be loaded and the output file. The former is the length of the upstream expanded region, and the latter is the length of the downstream expanded region.
- `<CORES>` : The analysis script can split for loop processing into chunks and run them in parallel across multiple cores to reduce run time. This argument specifying the number of cores to use for execution.

Note:
- Arguments other than the core count are used to create input and output file paths, so they do not impact analysis except for the presence or absence of input files.
- When inputting data that follows naming rules different from the example provided, the script must be modified to match those naming rules.

## Output file
This R script outputs two files: 
- The first table combines an overview of z-score calculations with regioneR results for significant genes. Note that p-values, average and standard deviation of overlaps, and z-scores were calculated separately using regioneR and regioneReloaded; therefore, values may not match due to differences in tools or the number of replicates.
- The second table outputs the raw data for z-score calculations. It lists the window offset distance and the corresponding z-score for each gene.


