# Sergio Alías, 20230516
# Last modified 20230516

# STAGE 1 OBTAINING COUNTS FROM FASTQ FILES

. ~soft_bio_267/initializes/init_autoflow


while IFS= read sample; do
AutoFlow -w $TEMPLATE -o "$COUNT_RESULTS_FOLDER"/"$sample" "$RESOURCES"
done < $SAMPLES_FILE
