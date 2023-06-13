#! /usr/bin/env bash

# Sergio Alías, 20230323
# Last modified 20230612

# daemon.sh

# Script for controlling the single-cell Cell Ranger workflow


framework_dir=`dirname $0`
export CODE_PATH=$(readlink -f $framework_dir )
CONFIG_DAEMON=$CODE_PATH'/config_daemon'
export module=$1 # For setting global vars from config_daemon according to the stage
source $CONFIG_DAEMON
export PATH=$CODE_PATH'/reports:'$PATH
export PATH=$CODE_PATH'/scripts:'$PATH
export PATH=$CODE_PATH'/aux_sh:'$PATH


## STAGE EXECUTION
#######################################################################

if [ "$module" == "0" ] ; then
    # STAGE 0 CONVERTING BCL FILES INTO FASTQ
    echo "Launching stage 0: Converting BCL files into FASTQ"
    if [ $launch_login == TRUE ]; then  
        cellranger_mkfastq.sh
    else
        sbatch aux_sh/cellranger_mkfastq.sh
    fi

elif [ "$module" == "1" ] ; then
    # STAGE 1 OBTAINING COUNTS FROM FASTQ FILES
    echo "Launching stage 1: Obtaining counts from FASTQ files"
    autoflow_launcher.sh count

elif [ "$module" == "2a" ] ; then
    # STAGE 2 QUALITY CONTROL AND TRIMMING
    echo "Launching stage 2: Quality control and trimming"
    autoflow_launcher.sh qc

elif [ "$module" == "2b" ] ; then
    # STAGE 2 SAMPLES COMPARISON
    echo "Launching stage 2: Samples comparison"
    if [ $launch_login == TRUE ]; then  
        compare_samples.sh
    else
        sbatch aux_sh/compare_samples.sh
    fi

elif [ "$module" == "3" ] ; then
    # STAGE 3 PREPROCESSING
    echo "Launching stage 3: Preprocessing"
    autoflow_launcher.sh preproc

fi