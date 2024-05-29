#!/usr/bin/env bash

################################################
# PROGRAM:
# analyze.sh
#
# DESCRIPTION:
# This is a very basic RNA-seq pipeline that I use for analyzing fastq reads. Step 1 is a
# simple wrapper that performs quality control, genome alignment, basic format
# conversions, and htseq-count tabulation for paired-end RNA-seq samples using a specified
# genome. Step2 is a clean up program that removes unnecessary files and compresses files
# to save space.
#
# AUTHOR:
# Erin Nishimura
# John Lemas
#
# START DATE:
# 5 January, 2024
#
# DEPENDENCIES:
#	Requires the installation of the follwing software: 
#		fastp
#		hisat2 (STAR)
#		featureCounts (subread)
#		samtools
#		bedtools
#		deep-tools
#		NOTE: This script will be edited to also use STAR for alignment
#		STAR
#
#	Requires access to the Nishimura Lab installed software sources on SUMMIT.
#
# REQUIRES:
#	 INPUT: .fastq files.	 For each sample, paired forward and reverse sequencing files
#								are required. These should be placed in an input
#								directory.
#
#	 INPUT: _metadata.txt file: A metadata file with two columns. The first two columns
#								are fastq file names. The third column is a "nickname"
#								of each sample. Later columns can be included with other
#								metadata information. Metadata file should be placed
#								within the inputdir directory.
#
#
#	 HISAT2 INDEXES: .ht2 files for the genome. These are produced using hisat2-build. For
#								instructions see
#			   https://ccb.jhu.edu/software/hisat2/manual.shtml#the-hisat2-build-indexer
#
#	 GENOME SEQUENCE: .fa  or .tar.gz file for the genome. This is the sequence of the 
#								 genome.
#
#	 GENOME ANNOTATION: .gtf file for the genome. This is a genome annotation file of gene
#								features. Version and coordinates must match the genome
#								sequence (.fa above).
#
# USAGE:
# $ bash analyze_RNAseq_231126.sh <metadata.txt> <number of threads>
#
# OUTPUT:
#
# KNOWN BUGS:
# When using a gff file for featureCounts make sure to use <-g 'ID'> to specify ID instead 
# of gene_id. gene_id is used in gtf files. If you dont make this alteration nothing will
# map in featureCounts.
# I think I need to use <-g 'Parent'> otherwise it will not count annotations labeled as genes.
# THINGS TO IMPROVE:
#
################################################


####### MODIFY THIS SECTION #############

#The input samples live in directory:
inputdir="/scratch/alpine/jolemas@colostate.edu/RT_stress_trial/01_input/fastq_files"

#Metadata file. This pulls the metadata path and file from the command line
metadata=$1

#This is where the ht2 files live:
hisat2path="/scratch/alpine/jolemas@colostate.edu/RT_stress_trial/03_output/02_hisat2ind/SalTr_v01"

#This is where the STAR files live
STARpath="/scratch/alpine/jolemas@colostate.edu/RT_stress_trial/03_output/01_STAR_genome"

#This is where the genome sequence lives:
genomefa="/scratch/alpine/jolemas@colostate.edu/RT_stress_trial/01_input/SalTr_Ref_v01.fasta"

#This is where the gtf file lives:
gtffile="/scratch/alpine/jolemas@colostate.edu/RT_stress_trial/01_input/SalTr_Ref_v01.0.gff.gff3"

#This is the output_directory:
#DATE=`date +%Y-%m-%d`
#OR
DATE='2024-05-13'

outputdir="../03_output/"$DATE"_output/"


########## DONE MODIFYING ###############



########## BEGIN CODE ###############

echo -e ">>> INITIATING analyzer with command:\n\t$0 $@"

#Number of threads to use:
# Note - this imports the number of threads (ntasks) given in the command line
pthread=$2

# Make output directories
echo -e ">>> MAKING output directory"
mkdir -p $outputdir



####### META DATA #############

#These are the sample names, R1:
samples1=( $(cut -f 1 --output-delimiter=' ' $metadata) )

#These are the sample names, R2:
samples2=( $(cut -f 2 --output-delimiter=' ' $metadata) )

#These are the nicknames I want to give the files:
names=( $(cut -f 3 --output-delimiter=' ' $metadata) )



####### PIPELINE ##############

# Report back to the user which files will be processed and which names they'll be given:
echo -e ">>> INPUT: This script will process files from the metafile:\n\t$metadata"
echo -e ">>> PLAN: This script will process the sample files into the following names: "
echo -e "\tSAMPLE1\tSAMPLE2\tNAMES"

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do
	echo -e "\t${samples1[$counter]}\t${samples2[$counter]}\t${names[$counter]}"
done


# FASTP to remove unwanted sequences
# FASTP to determine quality
echo -e "\n>>> FASTP: Trimming excess and low-quality sequences from .fastq file; generating quality report"
mkdir -p $outputdir"01_fastp"

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do
	samplename=${names[$counter]}
	sample1=${samples1[$counter]}
	sample2=${samples2[$counter]}

	## Echo statements##	
			##### ENTER ECHO STATEMENTS HERE #####
	
	## Make output directories
	mkdir -p $outputdir"01_fastp/"$samplename
	
	## execute fastp
	cmd1="fastp -i $inputdir/$sample1 \
-I $inputdir/$sample2 \
-o ${outputdir}01_fastp/${samplename}/${samplename}_trim_1.fastq \
-O ${outputdir}01_fastp/${samplename}/${samplename}_trim_2.fastq \
-h ${outputdir}01_fastp/${samplename}/${samplename}_report.html \
-j ${outputdir}01_fastp/${samplename}/${samplename}_report.json \
--detect_adapter_for_pe \
--thread $pthread \
-x -g "
	
	echo -e "\t$ ${cmd1}"
	time eval $cmd1

done

# HISAT2 to align to the genome
#echo -e "\n>>> HISAT2: aligning each sample to the genome"
#outhisat2=$outputdir"02_hisat2/"
#mkdir -p $outhisat2

#for (( counter=0; counter < ${#samples1[@]}; counter++ ))
#do
#	samplename=${names[$counter]}
#	sample1=${samples1[$counter]}
#	sample2=${samples2[$counter]}
#

	## execute hisat2
#	cmd2="hisat2 -x $hisat2path -1 $outputdir"01_fastp/"$samplename/$samplename"_trim_1.fastq" -2 $outputdir"01_fastp/"$samplename/$samplename"_trim_2.fastq" -S ${outhisat2}${samplename}.sam --summary-file ${outhisat2}${samplename}_summary.txt --no-unal -p $pthread"
#	echo -e "\t$ $cmd2"
#	time eval $cmd2

#done

# STAR alignment to the genenome:

echo -e "\n>>> STAR: aligning each sample to the genome"
outSTAR=$outputdir"03_STAR/"
mkdir -p $outSTAR

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do
	samplename=${names[$counter]}
	sample1=${samples1[$counter]}
	sample2=${samples2[$counter]}

	# execute STAR:
	cmd3="STAR --runThreadN ${pthread} --genomeDir ${STARpath} --readFilesIn $outputdir"01_fastp/"$samplename/$samplename"_trim_1.fastq" $outputdir"01_fastp/"$samplename/$samplename"_trim_2.fastq" --outFileNamePrefix ${outSTAR}${samplename}"
	echo -e "\t$ $cmd3"
	time eval $cmd3

done

# Rename STAR output for featureCounts:

echo -e "Renaming the STAR aligned output files for featureCounts..."

for (( counter=0; counter < ${#names[@]}; counter++ ))
do
  	echo -e "renaming ${names[$counter]}Aligned.out.sam..."
        mv ${outSTAR}${names[$counter]}Aligned.out.sam ${outSTAR}${names[$counter]}.sam
        echo "..finished. Starting next file.."
done

echo "..Renaming complete. Moving on to featureCounts.."

# FEATURECOUNTS to tabulate reads per gene:
echo -e "\n>>> FEATURECOUNTS: Run featureCounts on all files to tabulate read counts per gene"
outfeature=$outputdir"04_feature/"
mkdir -p $outfeature

# Acquire a list of .sam names
samfilePath=()
for (( counter=0; counter < ${#names[@]}; counter++ ))
do
	samfile=${names[$counter]}.sam
	samfilePath+=(${outSTAR}${samfile})

done


# Execute featureCounts
cmd4="featureCounts -p -T $pthread -a $gtffile -g 'Parent' -o ${outfeature}STAR_counts_v2.txt ${samfilePath[*]}"
echo -e "\t$ $cmd4"
time eval $cmd4




# SAMTOOLS and BAMCOVERAGE: to convert .sam output to uploadable .bam and .wg files
echo -e "\n>>> SAMTOOLS/BAMCOVERAGE: to convert files to uploadable _sort.bam and _sort.bam.bai files:"
samout=$outputdir"05_samtools/"
mkdir -p $samout

for seqname in ${names[@]}
do
	# echo
	echo -e "\tSamtools and BamCoverage convert: ${seqname}"
	
	# Samtools: compress .sam -> .bam
	cmd5="samtools view --threads $pthread -bS ${outSTAR}${seqname}.sam > ${samout}${seqname}.bam"
	echo -e "\t$ ${cmd5}"
	time eval $cmd5

	
	# Samtools: sort .bam -> _sort.bam
	cmd6="samtools sort --threads $pthread -o ${samout}${seqname}_sort.bam --reference $genomefa ${samout}${seqname}.bam"
	echo -e "\t$ ${cmd6}"
	time eval $cmd6
	
	
	# Samtools: index _sort.bam -> _sort.bam.bai
	cmd7="samtools index ${samout}${seqname}_sort.bam"
	echo -e "\t$ ${cmd7}"
	time eval $cmd7
	
	
	# bamCoverage: Create a .bw file that is normalized. This can be uploaded to IGV or UCSC
	cmd8="bamCoverage -b ${samout}${seqname}_sort.bam -o ${samout}${seqname}_sort.bw --outFileFormat bigwig -p $pthread --normalizeUsing CPM --binSize 1"
	echo -e "\t$ ${cmd8}"
	time eval $cmd8
	
done




######## VERSIONS #############
echo -e "\n>>> VERSIONS:"
echo -e "\n>>> FASTP VERSION:"
fastp --version
echo -e "\n>>> HISAT2 VERSION:"
hisat2 --version
echo -e "\n>>> SAMTOOLS VERSION:"
samtools --version
echo -e "\n>>> FEATURECOUNTS VERSION:"
featureCounts -v
echo -e "\n>>> BAMCOVERAGE VERSION:"
bamCoverage --version
echo -e ">>> END: Analayzer complete."
