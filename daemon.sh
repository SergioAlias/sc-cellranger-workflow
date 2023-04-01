#! /usr/bin/env bash

# Sergio Alías, 20230323
# Last modified 20230401

# daemon.sh

# Script for controlling the single-cell Cell Ranger workflow


CONFIG_DAEMON='config_daemon'
module=$1
source $CONFIG_DAEMON


## STAGE EXECUTION
#######################################################################

if [ "$module" == "0" ] ; then
	#STAGE 0 CONVERTING BCL FILES INTO FASTQ
	echo "Launching stage 0: Converting BCL files into FASTQ"
	if [ $launch_login == TRUE ]; then	
		$aux_folder/cellranger_mkfastq.sh $experiment_name
	else
		sbatch $aux_folder/cellranger_mkfastq.sh $experiment_name
	fi

elif [ "$module" == "1" ] ; then
	#STAGE 1 OBTAINING COUNTS FROM FASTQ FILES
	echo "Launching stage 1: Obtaining counts from FASTQ files"
	if [ $launch_login == TRUE ]; then	
		$aux_folder/cellranger_count.sh $experiment_name
	else
		sbatch $aux_folder/cellranger_count.sh $experiment_name
	fi
