# Sergio Alías, 20230516
# Last modified 20230517

# STAGE 1 OBTAINING COUNTS FROM FASTQ FILES

. ~soft_bio_267/initializes/init_autoflow


while IFS= read sample; do
    
    AF_VARS=`echo "
    \\$sample=$sample
    " | tr -d [:space:]`

    AutoFlow -w $TEMPLATE -V "$AF_VARS" -o "$COUNT_RESULTS_FOLDER"/"$sample" "$RESOURCES"

done < $SAMPLES_FILE