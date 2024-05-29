#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --partition=amilan
#SBATCH --time=06:00:00
#SBATCH --qos=normal
#SBATCH --ntasks=8
#SBATCH --job-name=fq_gunzip
#SBATCH --output=../03_output/log_files/log-fq_gunzip_test-%j.txt
#SBATCH --mail-type=END
#SBATCH --mail-user=jolemas@colostate.edu

DATE=`date +%Y-%m-%d`

echo -e "Assigning Paths.."

raw="/scratch/alpine/jolemas@colostate.edu/RT_stress_trial/01_input/raw_files"

files="/scratch/alpine/jolemas@colostate.edu/RT_stress_trial/01_input/fastq_files"

metadata="/scratch/alpine/jolemas@colostate.edu/RT_stress_trial/01_input/S.tragus_metadata_v2.txt"

# I need to unzip all of the fastq files and place them into a different directory
# 1st I need to do variable assignment for each sample read

#These are the sample names, R1:
samples1=( $(cut -f 1 $metadata) )

#These are the sample names, R2:
samples2=( $(cut -f 2 $metadata) )

#These are the nicknames (directory names) I want to give the files:
names=( $(cut -f 3 $metadata) )

# Test out the variable assignment

echo -e "These are the names of the directories and files that will be unzipped:"

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do
	echo -e "\t${samples1[$counter]}\t${samples2[$counter]}\t${names[$counter]}"
done



# then initiate the loop to unzip all of the files

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do
	echo -e "unzipping ${samples1[$counter]} and ${samples2[$counter]}"
	gunzip -c ${raw}/${names[$counter]}/${samples1[$counter]} > ${files}/${samples1[$counter]}
	gunzip -c ${raw}/${names[$counter]}/${samples2[$counter]} > ${files}/${samples2[$counter]}
	echo "..finished. Starting next set.."
done
