# Sergio Alías, 20230516
# Last modified 20230920

# Generic Autoflow launcher

. ~soft_bio_267/initializes/init_autoflow


if [ "$1" == "count" ] ; then # STAGE 1 OBTAINING COUNTS FROM FASTQ FILES
    export RESULTS_FOLDER=$COUNT_RESULTS_FOLDER
    export TEMPLATE=$COUNT_TEMPLATE
elif [ "$1" == "qc" ] ; then # STAGE 2 QUALITY CONTROL AND TRIMMING
    export RESULTS_FOLDER=$QC_RESULTS_FOLDER
    export TEMPLATE=$QC_TEMPLATE
elif [ "$1" == "preproc" ] ; then # STAGE 3 PREPROCESSING
    export RESULTS_FOLDER=$PREPROC_RESULTS_FOLDER
    export TEMPLATE=$PREPROC_TEMPLATE
fi

mkdir -p $RESULTS_FOLDER

PATH=$LAB_SCRIPTS:$PATH

export S_NUMBER=6

while IFS= read sample; do
    if [ "$multi_lane" == "FALSE" ] || [ "$1" != "qc" ] ; then
        AF_VARS=`echo "
        \\$sample=$sample,
        \\$s_num=$S_NUMBER,
        \\$lane_num=1
        " | tr -d [:space:]`
        AutoFlow -w $TEMPLATE -V "$AF_VARS" -o $RESULTS_FOLDER/$sample #$RESOURCES
    elif [ "$multi_lane" == "TRUE" ] ; then
        total_fastq_files=$(ls $read_path -1 | grep $sample | wc -l)
        number_of_lanes=$((total_fastq_files / 2))
        for (( i = 1; i <= $number_of_lanes; i++ )) ; do
            if [ "$sample" = "Undetermined" ] ; then
                S_NUMBER=$(( S_NUMBER - S_NUMBER ))
            fi
            AF_VARS=`echo "
            \\$sample=$sample,
            \\$s_num=$S_NUMBER,
            \\$lane_num=$i
            " | tr -d [:space:]`
            AutoFlow -w $TEMPLATE -V "$AF_VARS" -o $RESULTS_FOLDER/$sample #$RESOURCES
        done
        S_NUMBER=$(( S_NUMBER + 1 ))
    fi
done < $SAMPLES_FILE