#!/bin/bash
#-------------------------------------------------------------------------------------------------------------
# MGG_GH5 - INBREEDING
# script_get_data.sh
#-------------------------------------------------------------------------------------------------------------
#SBATCH -J get_data
#SBATCH -o get_data_%j.o 
#SBATCH -e get_data_%j.e
#SBATCH -t 00:10:00 # execution time
#SBATCH --mem=16GB # memory needed
#SBATCH -c 4
#-------------------------------------------------------------------------------------------------------------
# Load modules
module load cesga/2020 bcftools/1.19
#-------------------------------------------------------------------------------------------------------------
# Check number of arguments
if [ $# -ne 1 ]  
then
	echo "Usage: $0 <POP>" 
	exit 1
fi
# Set arguments
POP=$1 # population name

#-------------------------------------------------------------------------------------------------------------

### WORKING DIRECTORY ###########################

WDIR=$PWD

### OUTPUT DIRECTORY ############################

if [ ! -d data ] ; then mkdir data ; fi
cd data

### 1000 GENOMES PROJECT ########################

# DOWNLOAD COMPRESSED VCF AND PANEL WITH IDS AND POPULATION CODES

wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

# TABLE WITH POPULATION INFORMATION

awk 'NR>1{print $1"\t"$2"\t"$3}' integrated_call_samples_v3.20130502.ALL.panel > sample_pop.tsv

# EXTRACT IDs FROM POP

awk -v p="$POP" '$2==p {print $1}' sample_pop.tsv > ${POP}.ids

# VCF SUBSET BY SAMPLES

bcftools view -S ${POP}.ids -Oz -o chr22_${POP}.vcf.gz ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#bcftools index chr22_${POP}.vcf.gz

# Clean
rm ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
rm integrated_call_samples_v3.20130502.ALL.panel






