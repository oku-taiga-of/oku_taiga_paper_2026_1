# Analysis script for SNP enrichment in gene regions based on positional information using regioneR
This directory provides R scripts and example data used for analysis of genomic region enrichment by using regioneR in the accompanying manuscript.

## Files
### Analysis R scripts
- `GB_regioneR.R`: When referring to the gene body as the gene region, use this one.
- `outGB_regioneR.R`: When the gene region is the expansion region, use this one.
Note:
- The two scripts differ only in analysis mode; they share the same software environment.
- These scripts are intended to be executed via the CLI.
- The analysis algorithms for the two scripts are identical, but when using the extended region, overlaps are counted as the difference between the extended region and the gene body region (see the manuscript's method section for details).
- Note that the command arguments at runtime differ between the two scripts (the actual commands are described later).
- When using extended regions, the loaded extended region information may exceed the endpoints of the chromosome in the referenced genome, potentially causing a warning. Within the script, the exceeding portion is trimmed, so it does not directly cause a problem.

### Example input file
- `example_inputSNP_Endurance_nonAFR_LD_r2_MUO2_0004887.tsv` : Inputting SNP data
- `gene_id_name_map.csv` : File for referencing gene symbol names using gene IDs as keys when outputting results
- `gene_region_ORF.bed` : `.bed` files for each gene region (gene body data)
- `gene_region_u5.0k_d1.5k.bed` : `.bed` files for each gene region (data for extended regions), File for the region expanded 5.0k bp upstream and 1.5k bp downstream.
Note:
- These are example input files, so they are the minimum files required for demonstration and reproducibility checks.
- Since the script references columns, when running it in a different file, ensure the file structure matches the one containing the columns.
- Regarding SNP data, to prevent errors, each SNP must be assigned a unique SNP ID.

### Environment records
- `renv.lock`: R package environment lockfile
- `ENV_sessionInfo.txt`
- `ENV_bioconductor_version.txt`
- `ENV_installed_packages.csv`

## Requirements
- R version: see `ENV_sessionInfo.txt`
- macOS or Linux (the scripts use `parallel::mclapply()`)
The analysis was performed using a MacBook Pro (M1 Max) with 64GB of memory. Depending on the number of cores used and the SNP data, the analysis completed within a few minutes to a few hours.

## Setup
Restore the R package environment using **renv**:
```bash
Rscript -e 'install.packages("renv", repos="https://cloud.r-project.org"); renv::restore(prompt = FALSE)'
```
When published on GitHub, `renv/library` is not included.
* When manually setting up the environment, refer to the information contained in the Environment records.

## Run CLI commands (example)
### When the gene body is defined as the gene region
`Rscript GB_regioneR.R <trait_category> <trait_name> <trait_ID> <LD_r> <LD_pop> <bed_tag> <CORES>`
Command examples corresponding to the sample files
`Rscript GB_regioneR.R Endurance MUO2 0004887 r2 nonAFR ORF 4`
### When the expansion region is defined as the gene region
`Rscript outGB_regioneR.R <trait_category> <trait_name> <trait_ID> <LD_r> <LD_pop> <bed_tag_upstream> <bed_tag_downstream> <CORES>`
Command examples corresponding to the sample files
`Rscript outGB_regioneR.R Endurance MUO2 0004887 r2 nonAFR 2.0 0.2 4`
### The command arguments
- <trait_category>, <trait_name>, <trait_ID>, <LD_r>, <LD_pop> : For each SNP dataset loaded: - Trait category - Trait abbreviation - Trait ID - r^2 threshold for LD pruning - Population used for LD pruning in the genomic dataset. Arguments for creating paths for SNP data to be loaded and output files.
- <bed_tag> : An argument used only when using the gene body, to be declared explicitly in the command.
- <bed_tag_upstream> <bed_tag_downstream> : Arguments used only when using the expansion region, used to create paths for the expansion　region file to be loaded and the output file. The former is the length of the upstream expanded region, and the latter is the length of the downstream expanded region.
- <CORES> : The analysis script can split for loop processing into chunks and run them in parallel across multiple cores to reduce run time. This argument specifying the number of cores to use for execution.
Note:
- Arguments other than the core count are used to create input and output file paths, so they do not impact analysis except for the presence or absence of input files.
- When inputting data that follows naming rules different from the example provided, the script must be modified to match those naming rules.

## Output file
This R script outputs two files: one containing results for all genes used in the analysis, and another containing results only for significant genes.Both files are saved in the same directory. To change this, modify the output file path in the script.


