#! /usr/bin/env bash

# Sergio Alias, 20230711
# Last modified 20230717

# STAGE 2 QUALIMAP

#SBATCH -J qualimap.sh
#SBATCH --cpus-per-task=10
#SBATCH --mem='60gb'
#SBATCH --constraint=cal
#SBATCH --time=3-00:00:00
#SBATCH --error=job.qmap.%J.err
#SBATCH --output=job.qmap.%J.out

# Setup

hostname

module load qualimap

mkdir -p $QMAP_RESULTS_FOLDER

cd $QMAP_RESULTS_FOLDER

while IFS= read -r name; do
  if [ -d "$COUNT_RESULTS_FOLDER/$name" ]; then
    echo -e "$name\t$COUNT_RESULTS_FOLDER/$name/cellranger_0000/$name/outs/possorted_genome_bam.bam" >> $QMAP_RESULTS_FOLDER"/qualimap_input_data"
  fi
done < $SAMPLES_FILE

unset DISPLAY

# Main

/usr/bin/time qualimap multi-bamqc --run-bamqc \
                                   -d $QMAP_RESULTS_FOLDER"/qualimap_input_data" \
                                   -outdir $QMAP_RESULTS_FOLDER
 
