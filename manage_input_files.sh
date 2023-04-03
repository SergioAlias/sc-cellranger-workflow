#! /usr/bin/env bash

# Sergio Alías, 20230403
# Last modified 20230403

# manage_input_files.sh

# Script for creating soft links for the input files
# Adapted from seoanezonjic/DEG_workflow/manage_input_files_example.sh 

source_folder=/path/to/raw/data
output_folder=`pwd`'/raw_data'


# For BCL files #########
# Remove / comment if you are using FASTQ files

ln -s $source_folder $output_folder


#########################



# For FASTQ files #######

# Remove / comment if you are using BCL files
# Make sure they follow this naming convention -> https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm

mkdir -p $output_folder


ln -s $source_folder'/raw_1_index.fastq.gz' $output_folder'/SampleName_S1_L001_I1_001.fastq.gz'
ln -s $source_folder'/raw_1_direct.fastq.gz' $output_folder'/SampleName_S1_L001_R1_001.fastq.gz'
ln -s $source_folder'/raw_1_reverse.fastq.gz' $output_folder'/SampleName_S1_L001_R2_001.fastq.gz'
ln -s $source_folder'/raw_2_index.fastq.gz' $output_folder'/SampleName_S1_L002_I1_001.fastq.gz'
ln -s $source_folder'/raw_2_direct.fastq.gz' $output_folder'/SampleName_S1_L002_R1_001.fastq.gz'
ln -s $source_folder'/raw_2_reverse.fastq.gz' $output_folder'/SampleName_S1_L002_R2_001.fastq.gz'


#########################