#!/bin/bash
#SBATCH -J microbecensus
#SBATCH -N1 --ntasks-per-node=22
#SBATCH --mem=150G
#SBATCH -t 08:00:00
#SBATCH -q 
#SBATCH -A 

source deactivate
source activate MicrobeCensus

python run_microbe_census.py -n 100000000 -t 20 -l 120 ${dir}/${name}.1.fastq.gz,${dir}/${name}.2.fastq.gz ./microbecensus_results/${name}.txt

########USAGE#########
#sbatch -o ./log/microbecensus_out.name --export dir=/path/to/reads,name=sample_name 04_microbecensus.sbatch
