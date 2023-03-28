#! /usr/bin/env bash

# Sergio Alías, 20230323
# Last modified 20230328

# daemon.sh

# Script for controlling the single-cell Cell Ranger workflow


module=$1


## STAGE EXECUTION
#######################################################################

if [ "$module" == "0" ] ; then
	#STAGE 0 CONVERTING BCL FILES INTO FASTQ
	echo "Launching stage 0: Converting BCL files into FASTQ"
	if [ $launch_login == TRUE ]; then	
		cellranger_mkfastq.sh
	else
		sbatch cellranger_mkfastq.sh
	fi

elif [ "$module" == "1" ] ; then
	#STAGE 1 OBTAINING COUNTS FROM FASTQ FILES
	echo "Launching stage 1: Obtaining counts from FASTQ files"
	if [ $launch_login == TRUE ]; then	
		cellranger_count.sh
	else
		sbatch cellranger_count.sh
	fi
