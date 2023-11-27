# Sergio Alías, 20230516
# Last modified 20231019

# Generic Autoflow launcher

. ~soft_bio_267/initializes/init_autoflow


if [ "$1" == "count" ] ; then # STAGE 1 OBTAINING COUNTS FROM FASTQ FILES
    export TEMPLATE=$COUNT_TEMPLATE
    export RESULTS_FOLDER=$COUNT_RESULTS_FOLDER
elif [ "$1" == "qc" ] ; then # STAGE 2 QUALITY CONTROL AND TRIMMING
    export TEMPLATE=$QC_TEMPLATE
    export RESULTS_FOLDER=$QC_RESULTS_FOLDER
elif [ "$1" == "preproc" ] ; then # STAGE 3a PREPROCESSING  
    export TEMPLATE=$PREPROC_TEMPLATE
    if [ "$integrative_analysis" == "TRUE" ] ; then
        export RESULTS_FOLDER=$INTEGRATE_RESULTS_FOLDER"/"$subset_column
    else
        export RESULTS_FOLDER=$PREPROC_RESULTS_FOLDER
    fi
fi

mkdir -p $RESULTS_FOLDER

if [ "$integrative_analysis" == "TRUE" ] ; then
    Rscript subset.R # TODO add args
    SAMPLES_FILE=
fi

PATH=$LAB_SCRIPTS:$PATH

export S_NUMBER=1

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