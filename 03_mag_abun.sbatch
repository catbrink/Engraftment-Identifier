#!/bin/bash
#SBATCH -J mag_abun
#SBATCH -N1 --ntasks-per-node=18
#SBATCH --mem=100G
#SBATCH -t 24:00:00
#SBATCH -q 
#SBATCH -A 

source deactivate
source activate CoverM
module load samtools/1.16.1

coverm genome -1 ${metagenome1} -2 ${metagenome2} --genome-fasta-directory ${genome_dir} -m trimmed_mean --trim-min 10 --trim-max 90 -o ${out} --min-read-percent-identity 95 --min-read-aligned-percent 70 --exclude-supplementary --output-format sparse

##########USAGE############
#sbatch -o ./log/mag_abun_out.${name} --export metagenome1=/path/to/forward,metagenome2=/path/to/reverse,genome_dir=/path/to/fasta/files,out=${output_name} 03_mag_abun.sbatch
