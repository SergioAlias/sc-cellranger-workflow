#! /usr/bin/env bash

# Sergio Alías, 20230323
# Last modified 20230421

# daemon.sh

# Script for controlling the single-cell Cell Ranger workflow


framework_dir=`dirname $0`
export CODE_PATH=$(readlink -f $framework_dir )
CONFIG_DAEMON=$CODE_PATH'/config_daemon'
export module=$1 # For setting global vars from config_daemon according to the stage
source $CONFIG_DAEMON
export PATH=$CODE_PATH'/reports:'$PATH
export PATH=$CODE_PATH'/aux_R:'$PATH
export PATH=$CODE_PATH'/aux_sh:'$PATH


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
elif [ "$module" == "2" ] ; then
    #STAGE 2 PREPROCESSING
    echo "Launching stage 2: Preprocessing"
    if [ $launch_login == TRUE ]; then  
        preprocess.sh
    else
        sbatch preprocess.sh
    fi
fi