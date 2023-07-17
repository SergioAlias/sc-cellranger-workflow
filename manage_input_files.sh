#! /usr/bin/env bash

# Sergio Alías, 20230403
# Last modified 20230716

# manage_input_files.sh

# Script for creating soft links for the input files
# Adapted from seoanezonjic/DEG_workflow/manage_input_files_example.sh 

source_folder=/mnt2/fscratch/users/bio_267_uma/sergioalias/all_LUADs
output_folder=`pwd`'/raw_data'


# For BCL files #########
# Remove / comment if you are using FASTQ files

# ln -s $source_folder $output_folder


#########################



# For FASTQ files #######

# Remove / comment if you are using BCL files
# Make sure they follow this naming convention -> https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm

mkdir -p $output_folder


ln -s $source_folder'/CRR073022_f1.fastq.gz' $output_folder'/LUAD1_S1_L001_R1_001.fastq.gz'
ln -s $source_folder'/CRR073022_r2.fastq.gz' $output_folder'/LUAD1_S1_L001_R2_001.fastq.gz'
ln -s $source_folder'/CRR073023_f1.fastq.gz' $output_folder'/LUAD2_S1_L001_R1_001.fastq.gz'
ln -s $source_folder'/CRR073023_r2.fastq.gz' $output_folder'/LUAD2_S1_L001_R2_001.fastq.gz'
ln -s $source_folder'/CRR073024_f1.fastq.gz' $output_folder'/LUAD3_S1_L001_R1_001.fastq.gz'
ln -s $source_folder'/CRR073024_r2.fastq.gz' $output_folder'/LUAD3_S1_L001_R2_001.fastq.gz'
ln -s $source_folder'/CRR073025_f1.fastq.gz' $output_folder'/LUAD4_S1_L001_R1_001.fastq.gz'
ln -s $source_folder'/CRR073025_r2.fastq.gz' $output_folder'/LUAD4_S1_L001_R2_001.fastq.gz'
ln -s $source_folder'/CRR073026_f1.fastq.gz' $output_folder'/LUAD5_S1_L001_R1_001.fastq.gz'
ln -s $source_folder'/CRR073026_r2.fastq.gz' $output_folder'/LUAD5_S1_L001_R2_001.fastq.gz'
ln -s $source_folder'/CRR073027_f1.fastq.gz' $output_folder'/N1_S1_L001_R1_001.fastq.gz'
ln -s $source_folder'/CRR073027_r2.fastq.gz' $output_folder'/N1_S1_L001_R2_001.fastq.gz'
ln -s $source_folder'/CRR073028_f1.fastq.gz' $output_folder'/N2_S1_L001_R1_001.fastq.gz'
ln -s $source_folder'/CRR073028_r2.fastq.gz' $output_folder'/N2_S1_L001_R2_001.fastq.gz'
ln -s $source_folder'/CRR073029_f1.fastq.gz' $output_folder'/N3_S1_L001_R1_001.fastq.gz'
ln -s $source_folder'/CRR073029_r2.fastq.gz' $output_folder'/N3_S1_L001_R2_001.fastq.gz'
ln -s $source_folder'/CRR073030_f1.fastq.gz' $output_folder'/N4_S1_L001_R1_001.fastq.gz'
ln -s $source_folder'/CRR073030_r2.fastq.gz' $output_folder'/N4_S1_L001_R2_001.fastq.gz'
ln -s $source_folder'/CRR073031_f1.fastq.gz' $output_folder'/N5_S1_L001_R1_001.fastq.gz'
ln -s $source_folder'/CRR073031_r2.fastq.gz' $output_folder'/N5_S1_L001_R2_001.fastq.gz'

#########################