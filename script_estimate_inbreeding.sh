#!/bin/bash
#-------------------------------------------------------------------------------------------------------------
# MGG_GH5 - INBREEDING
# script_estimate_inbreeding.sh
#-------------------------------------------------------------------------------------------------------------
#SBATCH -J inbreeding
#SBATCH -o inbreeding_%j.o 
#SBATCH -e inbreeding_%j.e
#SBATCH -t 00:10:00 # execution time
#SBATCH --mem=1GB # memory needed
#SBATCH -c 4
#-------------------------------------------------------------------------------------------------------------
# Load modules
module load cesga/2020 plink/1.9b5
#-------------------------------------------------------------------------------------------------------------
# Check number of arguments
if [ $# -ne 1 ]  
then
	echo "Usage: $0 <DATASET>" 
	exit 1
fi
# Set arguments
DATASET=$1 # Data name without file extension
#-------------------------------------------------------------------------------------------------------------

### DIRECTORIES #####################################

# Working directory
WDIR=$PWD

# Output directories
if [ ! -d results_inbreeding/$DATASET ] ; then mkdir -p results_inbreeding/$DATASET ; fi
cd results_inbreeding/$DATASET

### GENOMIC INBREEDING ##############################
# SNP-by-SNP basis (FLH, Fhom, FhatI, FhatII, FhatII)
#  * --het
#  * --ibc
# Runs Of Homozygosity (FROH)
#  * --homozyg

# ******************************
# --het (FLH, Fhom)
# ******************************

plink --bfile $WDIR/data_pruned/$DATASET --het --out $DATASET.Fhet

sed '1d' $DATASET.Fhet.het > temp1       # remove headers
awk '{print ($3/$5)}' temp1 > temp2      # calculate average number of homozygous SNPs
sed -i '1i Fhom' temp2                   # add header
paste $DATASET.Fhet.het temp2 > temp3    # add new column 'Fhom'
rm temp1 temp2 $DATASET.Fhet.het         # delete unnecessary files
mv temp3 $DATASET.Fhet.het               # rename the final output file

# ******************************
# --ibc (FhatI, FhatII, FhatII)
# ******************************

plink --bfile $WDIR/data_pruned/$DATASET --ibc --out $DATASET.Fibc

# ******************************
# FROH (ROH > 0.1 Mb)
# ******************************

plink --bfile $WDIR/data_pruned/$DATASET --homozyg-kb 100 --homozyg --out $DATASET.FROH01

sed '1d' $DATASET.FROH01.hom.indiv > temp1     # remove headers
awk '{print ($5/100000)}' temp1 > temp2       # calculate FROH (genome length 100Mb = 100000Kb)
sed -i '1i FROH01' temp2                       # add header
paste $DATASET.FROH01.hom.indiv temp2 > temp3  # add new column 'FROH'
rm temp1 temp2 $DATASET.FROH01.hom.indiv       # delete unnecessary files
mv temp3 $DATASET.FROH01.hom.indiv             # rename the final output file

# ******************************
# FROH (ROH > 1 Mb)
# ******************************

plink --bfile $WDIR/data_pruned/$DATASET --homozyg-kb 1000 --homozyg --out $DATASET.FROH1

sed '1d' $DATASET.FROH1.hom.indiv > temp1     # remove headers
awk '{print ($5/100000)}' temp1 > temp2       # calculate FROH (genome length 100Mb = 100000Kb)
sed -i '1i FROH1' temp2                       # add header
paste $DATASET.FROH1.hom.indiv temp2 > temp3  # add new column 'FROH'
rm temp1 temp2 $DATASET.FROH1.hom.indiv       # delete unnecessary files
mv temp3 $DATASET.FROH1.hom.indiv             # rename the final output file

# ******************************
# FROH (ROH > 5 Mb)
# ******************************

plink --bfile $WDIR/data_pruned/$DATASET --homozyg-kb 5000 --homozyg --out $DATASET.FROH5

sed '1d' $DATASET.FROH5.hom.indiv > temp1     # remove headers
awk '{print ($5/100000)}' temp1 > temp2       # calculate FROH (genome length 100Mb = 100000Kb)
sed -i '1i FROH5' temp2                       # add header
paste $DATASET.FROH5.hom.indiv temp2 > temp3  # add new column 'FROH'
rm temp1 temp2 $DATASET.FROH5.hom.indiv       # delete unnecessary files
mv temp3 $DATASET.FROH5.hom.indiv             # rename the final output file


### MERGE F RESULTS #################################

awk '{print $2 "\t" $7}' $DATASET.Fhet.het > temp1     # FHOM
awk '{print $6}' $DATASET.Fibc.ibc > temp2             # Fhat3 (or FYan)
awk '{print $7}' $DATASET.FROH01.hom.indiv > temp3     # FROH01
awk '{print $7}' $DATASET.FROH1.hom.indiv > temp4      # FROH1
awk '{print $7}' $DATASET.FROH5.hom.indiv > temp5      # FROH5
paste temp1 temp2 temp3 temp4 temp5 > $DATASET.inbreeding
rm temp1 temp2 temp3 temp4 temp5

