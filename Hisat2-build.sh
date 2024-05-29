#!/usr/env/bin bash
#SBATCH --nodes=1
#SBATCH --partition=amilan
#SBATCH --time=05:00:00
#SBATCH --qos=normal
#SBATCH --ntasks=12
#SBATCH --job-name=Hisat2-build_v1
#SBATCH --output=../03_output/log_files/log-Hisat2-build_v1-%j.txt
#SBATCH --mail-type=END
#SBATCH --mail-user=jolemas@colostate.edu

DATE=`date +%Y-%m-%d`

outputdir="../03_output/01_hisat2"

cmd1="source /curc/sw/anaconda3/latest"

echo -e "${cmd1}"

cmd1

cmd2="conda activate rnaseq"

echo -e "${cmd2}"

cmd2



