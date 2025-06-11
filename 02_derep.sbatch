#!/bin/bash
#SBATCH -J derep
#SBATCH -N1 --ntasks-per-node=22
#SBATCH --mem=250G
#SBATCH -t 48:00:00
#SBATCH -q 
#SBATCH -A 

source deactivate
source activate dRep

dRep dereplicate --completeness 50 --contamination 10 --S_algorithm skani -pa 0.95 -sa 0.995 --processors 20 derep_out -g ${dir}/*.fna

########USAGE#########
#sbatch -o ./log/derep_out.name --export dir=/path/to/MAGs 02_derep.sbatch
