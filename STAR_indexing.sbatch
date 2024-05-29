#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --partition=amilan
#SBATCH --time=03:00:00
#SBATCH --qos=normal
#SBATCH --ntasks=8 
#SBATCH --job-name=STAR_index
#SBATCH --output=../03_output/log_files/log-STAR_index-%j.txt
#SBATCH --mail-type=END
#SBATCH --mail-user=jolemas@colostate.edu

DATE=`date +%Y-%m-%d`

genome="/scratch/alpine/jolemas@colostate.edu/RT_stress_trial/01_input/SalTr_Ref_v01.fasta"

annotation="/scratch/alpine/jolemas@colostate.edu/RT_stress_trial/01_input/SalTr_Ref_v01.0.gff.gff3"

outputdir="/scratch/alpine/jolemas@colostate.edu/RT_stress_trial/03_output/01_STAR_genome"

# STAR genome indexing script

echo -e	"STAR indexing>>> This script will index a genome to generate sam files using STAR."

echo -e	"STAR indexing>>> This script was run on ${genome} on ${DATE} by ${USER}."

# open conda environment

echo -e "STAR indexing>>> Activating conda environment.."

cmd1="source /curc/sw/anaconda3/latest"

echo -e "$cmd1"

$cmd1

cmd2="conda activate rnaseq"

echo -e "$cmd2"

$cmd2

# Begin STAR indexing script

echo -e "STAR indexing>>> ..environment activated.\n\tBegin indexing.."

cmd3="STAR --runMode genomeGenerate --genomeFastaFiles ${genome} --runThreadN $SLURM_NTASKS --sjdbGTFfile ${annotation} --sjdbGTFtagExonParentTranscript Parent --genomeDir ${outputdir}"

echo -e "${cmd3}"

$cmd3

echo -e "STAR indexing>>> Indexing job completed. Exiting software.."







