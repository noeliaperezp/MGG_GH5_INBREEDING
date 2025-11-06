#!/bin/bash
#-------------------------------------------------------------------------------------------------------------
# MGG_GH5 - INBREEDING
# script_run_SLiM.sh
#-------------------------------------------------------------------------------------------------------------
#SBATCH -J slim
#SBATCH -o slim_%j.o 
#SBATCH -e slim_%j.e
#SBATCH -t 00:20:00 # execution time
#SBATCH --mem=1GB # memory needed
#SBATCH -c 4
#-------------------------------------------------------------------------------------------------------------
# Load modules
module load cesga/2020 gcccore/system slim/3.7.1
#-------------------------------------------------------------------------------------------------------------
# Check number of arguments
if [ $# -ne 1 ]  
then
	echo "Usage: $0 <INPUT>" 
	exit 1
fi
# Set arguments
INPUT=$1 # SLIM input
#-------------------------------------------------------------------------------------------------------------

### DIRECTORIES ###################

# Working directory
WDIR=$PWD

# Output directories
if [ ! -d data ] ; then mkdir data ; fi
cd data

### SLIM ##########################

slim $WDIR/$INPUT








