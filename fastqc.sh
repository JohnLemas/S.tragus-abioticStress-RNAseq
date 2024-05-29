#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --partition=amilan
#SBATCH --time=05:00:00
#SBATCH --qos=normal
#SBATCH --ntasks=8
#SBATCH --job-name=fastqc
0#SBATCH --output=../03_output/log_files/log-fastqc-%j.txt
#SBATCH --mail-type=END
#SBATCH --mail-user=jolemas@colostate.edu

# Capture the date of this run

DATE=`date +%Y-%m-%d`

# Create the output directories and the README file

outputdir="../03_output/03_fastqc"

mkdir $outputdir

touch ${outputdir}/README_fastqc.txt

# Define the input file path

inputdir="../01_input/fastq_files"

# Write information to the README file

echo -e "This operation was preformed on ${DATE} by ${USER}." >> ${outputdir}/README_fastqc.txt

echo -e "Files for these reports are found in ${inputdir}" >> ${outputdir}/README_fastqc.txt

# Activate conda environment

echo -e "..activating conda environment" >> ${outputdir}/README_fastqc.txt

cmd1="source /curc/sw/anaconda3/latest"

echo -e "${cmd1}" >> ${outputdir}/README_fastqc.txt

$cmd1

cmd2="conda activate rnaseq"

echo -e "${cmd2}" >> ${outputdir}/README_fastqc.txt

$cmd2

echo -e "..conda environment activated. Running fastqc.."

# Initiate fastqc

cmd3="fastqc -o ${outputdir} -f fastq ${inputdir}/*"

echo -e "${cmd3}" >> ${outputdir}/README_fastqc.txt

$cmd3

# Capture the version

echo -e "..fastqc completed. Output files can be found in ${outputdir}."

echo -e "Version for this run is.."

cmd4="fastqc -v" >> ${outputdir}/README_fastqc.txt

echo -e "${cmd4}"

$cmd4 

# Run multiqc on output

echo -e "Deactivating environment.."

cmd5="conda deactivate"

echo -e "${cmd5}" >> $outputdir/README_fastqc.txt

$cmd5

$cmd5

echo -e "..conda deactivated. Purging any modules.."

cmd6="module purge"

echo -e "${cmd6}" >> $outputdir/README_fastqc.txt

$cmd6

echo -e "..modules purged. Initiating multiqc.."

cmd7="module load multiqc/1.14"

echo -e "${cmd7}" >> $outputdir/README_fastqc.txt

$cmd7

echo -e "..multiqc loaded. Beginning run on files in ${outputdir}.."

cmd8="multiqc ${outputdir}"

echo -e "${cmd8}" >> $outputdir/README_fastqc.txt

$cmd8

echo -e "..fastq analysis complete. All files can be found in ${outputdir}." >> $outputdir/README_fastqc.txt

echo "Exiting software.." >> $outputdir/README_fastqc.txt