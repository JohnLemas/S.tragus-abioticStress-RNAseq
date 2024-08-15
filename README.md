# S.tragus-abioticStress-RNAseq

## Overview

This project was performed to understand abiotic stress response in the allotetraploid _Salsola tragus_. A reference genome and annotation file for haplome 1 of _S. tragus_ was obtained from [WeedPedia](https://weedpedia.weedgenomics.org/). 

## Differential expression pipeline

Pipeline scripts were adapted from those by [Dr. Erin Nishimura](github.com/erinosb/2023_RNAseq_pairedend). Adaptations for this project from Erin's pipeline include substituting STAR for Hisat2 during alignment. This pipline was run at [CU Boulder's Research Computing](https://www.colorado.edu/rc/). 

## DESeq2 Scripts 

The DESeq2 scripts for this analysis were adapted from the RNA sequencing course offered at Colorado State University. The original author is [Dr. Erin Nishimura](github.com/erinosb/2023_RNAseq_pairedend).

Erin's script served as a basis for the DESeq2 design. Several versions were created to accomidate the specific experimental design for this study.

# The Experimental Design

This experiment was performed to understand abiotic stress response in the tetraploid tumbleweed _Salsola tragus_. The analysis included a transcriptome analysis of both root and shoot tissue under salinity stress and temperature stress over a 24 hour period. Each sampling time included three biological replicates. Three root and three shoot tissues were sampled at 0 hours after treatment to serve as a baseline control. The salinity stress treatment included three biological reps at 3, 8, and 24 hours after treatment. The cold stress treatment included three biological reps at 6 and 24 hours after treatment. RNA was extracted from these samples and sequenced on Illumina resulting in 72 paired-end reads. These reads wre run through the differential expression pipeline to obtain a counts file from featureCounts for use in RStudio DESeq2. 

## RNA read QC, alignment, and counts data table generation

The original `analyze.sh` script was obtained from [Dr. Erin Nishimura](github.com/erinosb/2023_RNAseq_pairedend). This original script was worked through several versions resulting in `analyze_v9.sh`. [CURC](https://www.colorado.edu/rc/) required SLURM commands for batch scripting so an additional `execute.sbatch` script was utilized to run the `analyze.sh` script. `execute.sbatch` followed version creation for `analyze.sh` resulting in `execute_v9.sbatch`. On a further note, SLURM at CURC preferred `.sbatch` to `.sh` when submitting batch jobs. Scripts in this repository may include either extension.

Prior to running `execute_v9.sh`, the reference genome required indexing for RNA sequence alignment to the genome. `Hisat2-build.sh` and `STAR_indexing.sbatch` were used to index the reference genome for use in either alignment. Alignment was attempted with both `hisat2` and `STAR` to decide on the best alignment software. The final analysis only included `STAR`. 

`execute_v9.sh` was then submitted resulting in several output directories for each of the analysis programs. Quality control analys was preformed using `fastqc.sh`. `STAR_counts_v6.sh` was then used in DESeq2 as follows.

## Differential expression analysis in DESeq2

Due to the complexity of this analysis several seperate differential expression contrasts were performed. `STAR_counts_v6.txt` was partitioned into several seperate DESeq2 scripts in order to contrast between specific plant tissue types and timepoints. For example, `SalTr_DESeq2_Root_v4_24.R` and `SalTr_DESeq2_Shoot_v4_24.R` subset samples from `STAR_counts_v6.txt` to perform differential expression on root or shoot tissues only at the 0 versus 24 hours after treatment timepoints. DESeq2 scripts that contain `_EV.R` are those that were used to generate volcano plots via [Enhanced Volcano](https://github.com/kevinblighe/EnhancedVolcano). `SalTr_DESeq2_v4_24_mishra.R` was used to identify abiotic stress related genes that were reported in:

Mishra, A., & Tanna, B. (2017). Halophytes: Potential Resources for Salt Stress Tolerance Genes and Promoters. Frontiers in Plant Science, 8, 829. https://doi.org/10.3389/fpls.2017.00829


