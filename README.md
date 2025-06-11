# Engraftment-Identifier
Scripts and programmes used to identify engraftment in genomic features as used in the paper: Identifying a prospective consortium of bacterial species and functions for the reduction of enteric multidrug resistant organisms by microbiota transplantation


The basic function of this repository is to descrive how to use the python script "identify_engrafting_features.py".

#STEP 00: gather data
Theoretically, the script should work on any genomic data, but it was specifically tested on shotgun metagenomic sequences.
Metagenomic sequences should be trimmed to remove adapter sequences prior to analysis. There are several programmes that can be used, including bbduk and trimmomatic.
#STEP 01: calculate genome equivalents
Genome equivalents are an estimate of the expected species diversity in a sample. This gives us a way to normalise coverage values across different metagenomic samples. We calculate genome equivalents using MicrobeCensus. You can use the following script and command to calculate GEQs:
>sbatch -o ./log/microbecensus_out.${name} --export dir=/path/to/reads,name=sample_name 01_microbecensus.sbatch 
#patient_name is the filename assigned to the paired metagenomic reads. This script assumes reads are saved as: ${dir}/${sample_name}.[12].fastq.gz
#this script will result in a folder called "microbecensus_results" with a file for each paired sample. You can create a combined file using the following bash command if you provide a file containing a list of sample names:
>while read sample_name; do geq=`grep "genome_equivalents" ${sample_name}.txt`; echo -e "${sample_name}\t${geq}" >> GEQ_all.tsv; done < sample_name.list
#STEP 02: Assemble and bin metagenomes
The paper worked with metagenome assembled genomes (MAGs) as the primary genomic feature that was begin tracked. MAGs are created from metagenomic reads by assembling (using e.g. SPAdes) and binning (using e.g. Maxbin2) metagenomes. 
The resulting assemblies were further processed using dRep to capture species representatives using the following script:
>sbatch -o ./log/derep_out.${name}  --export dir=/path/to/MAGs 02_derep.sbatch
#the script saves medium to high-quality MAGs defined as completeness >50%, and contamination <10%, and clusters MAGs as species representatives at the ANI threshold of 95%.
#STEP 03: Calculate MAG/feature abundance
>sbatch -o ./log/mag_abun_out.${name} --export metagenome1=/path/to/forward_reads,metagenome2=/path/to/reverse_reads,genome_dir=/path/to/drep/fasta/files,out=${output_name} 03_mag_abun.sbatch
#this script uses CoverM to calculate truncated average depth (TAD) excluding the highest and lowest 10% of matches (TAD80)
#STEP 03.5: Concatenate patient files
#concatenate all the files belonging to a single patient by using the following command:
>while read uniqueID; do x=`basename $uniqueID | cut -d- -f1`; y=`basename $uniqueID | cut -d- -f2`; tail -q -n +2 ./coverm_out/${x}.tsv >> ./concat_mag_abun/${y}_abundance.tsv; done < patient_concat.list
#where "patient_concat.list" is a matching sample names to patient_IDs with the format ${sample_name}-${patient_id} as follows:
    SRR14868472_1_kneaddata_paired-PM01
    SRR14868471_1_kneaddata_paired-PM01
    SRR14868467_1_kneaddata_paired-PM01
    SRR14868466_1_kneaddata_paired-PM01
    SRR14868465_1_kneaddata_paired-PM01
#STEP 04: Identify engrafting MAGs
This script requires 3 input parameters:
- Abundance file: the concatenated CoverM output produced in Step 03.5
- Metadata file: this is a user-supplied tab-separated file that summarises the sample names, patient IDs and days post FMT for each sample file. With the following format:
  filename  patient_id  days_post_fmt
  Ge_R11_V2  RGMGeR11  0
  Ge_R11_V3  RGMGeR11  28
  #Days-post-FMT are expected to be labelled as follows:
    "0" = baseline samples
    "15" = 15 days after baseline in a patient that recieved FMT
    "-15" = 15 days after baseline in a patient that did not recieve FMT (control/placebo sample)
- GEQ file: the combined GEQ_all.tsv file produced in STEP 01
- Sample type: specify whether the patient is a "Recipient" (recieved FMT) or "Control" (no FMT/placebo). Default: Recipient

The python script is run as follows:
>python ~/scratch/2025/python_scripts/01-1a_identify_engrafting_mags_control.py -a ${abundance_file} -m ${metadata_file} -g ${geq_file} -s ${Recipient/Control} -o ${out}

It will produce two output files:
${out}_sorted_filtered.tsv #a horizontal matrix of the abundance of all MAGs over time in one patient
${out}_engrafting.tsv #a horizontal matrix of the abundance of MAGs identified as engrafting in a patient
  
