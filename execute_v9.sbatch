#!/usr/bin/env bash

#SBATCH --job-name=SalTr_pipeline_v9
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=12:00:00
#SBATCH --partition=amilan
#SBATCH --mail-type=END
#SBATCH --mail-user=jolemas@colostate.edu
#SBATCH --output=../03_output/log_files/log-SalTr_pipeline_v9-%j.out

######### Instructions ###########

# Modify your SLURM entries above to fit your choices
# Below, modify the SECOND argument to point to YOUR metadata.file
# Below, you don't need to change $SLURM_NTASKS. It will automatically populate whatever you put in --ntasks=# above.
# Execute this script using $ sbatch execute_RNAseq_pipeline.sbatch

metadata="/scratch/alpine/jolemas@colostate.edu/RT_stress_trial/01_input/S.tragus_metadata_v2.txt"

## Initiate conda environment:

source /curc/sw/anaconda3/latest

conda activate rnaseq

## Execute the RNA-seq_pipeline to run the pipeline
bash analyze_v9.sh ${metadata} $SLURM_NTASKS


## Execute the cleanup script to zip .fastq files and delete extra files
#bash RNAseq_cleanup_231126.sh <metadatafile>
