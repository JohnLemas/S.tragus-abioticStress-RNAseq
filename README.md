# S.tragus-abioticStress-RNAseq

## Differential expression pipeline

The pipeline scripts were adapted from scripts found Erin Nishimura's repository erinosb/2023_RNAseq_pairedend.  

Major adaptations include substituting STAR for Hisat2 during alignment. 

## DESeq2 Scripts 

The DESeq2 scripts for this analysis were adapted from the RNA sequencing course offered at Colorado State University. The original author is Dr. Erin Nishimura (@erinosb).

This script was highly edited to accomidate the _Salsola tragus_ abiotic stress analysis. 

# The Experimental Design

This experiment was performed to understand abiotic stress response in the tetraploid tumbleweed _Salsola tragus_. The analysis included a transcriptome analysis of both root and shoot tissue under salinity stress and temperature stress over a 24 hour period. Each sampling time included three biological replicates. Three root and three shoot tissues were sampled at 0 hours after treatment to serve as a baseline control. The salinity stress treatment included three biological reps at 3, 8, and 24 hours after treatment. The cold stress treatment included three biological reps at 6 and 24 hours after treatment. RNA was extracted from these samples and sequenced on Illumina resulting in 72 paired-end reads. These reads wre run through the differential expression pipeline to obtain a counts file from featureCounts for use in RStudio DESeq2. 
