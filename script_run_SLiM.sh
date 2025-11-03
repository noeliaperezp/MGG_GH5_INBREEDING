#!/bin/bash
#-------------------------------------------------------------------------------------------------------------
# MGG_GH5 - INBREEDING
# script_run_SLiM.sh
#-------------------------------------------------------------------------------------------------------------
#SBATCH -J slim
#SBATCH -o slim_%j.o 
#SBATCH -e slim_%j.e
#SBATCH -t 00:10:00 # execution time
#SBATCH --mem=10GB # memory needed
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
if [ ! -d results_slim/$INPUT ] ; then mkdir -p results_slim/$INPUT ; fi
cd results_slim/$INPUT

### SLIM ##########################

slim $WDIR/$INPUT








