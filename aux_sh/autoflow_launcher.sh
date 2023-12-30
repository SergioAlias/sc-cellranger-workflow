# Sergio Alías, 20230516
# Last modified 20231226

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
    export RESULTS_FOLDER=$PREPROC_RESULTS_FOLDER

fi

mkdir -p $RESULTS_FOLDER

if [ "$1" == "preproc" ] && [ "$integrative_analysis" == "TRUE" ] ; then # TODO uncomment below
    . ~soft_bio_267/initializes/init_R
    Rscript scripts/prior_integration.R --exp_design $exp_design \
                                        --output $RESULTS_FOLDER \
                                        --condition $subset_column \
                                        --integration_file $integration_file \
                                        --experiment_name $experiment_name \
                                        --count_folder $COUNT_RESULTS_FOLDER
    export SAMPLES_FILE=$integration_file
fi

PATH=$LAB_SCRIPTS:$PATH

while IFS= read sample; do
    if [ "$1" == "qc" ] ; then # we need to iterate just in case we use FASTQC with multiple lanes (that means more FASTQ per sample)
        number_of_lanes=$(ls $read_path -1 | grep $sample | grep '_R1_001.fastq.gz' | wc -l)
    else # when not using the FASTQ files we set the iteration variable to 1 so we run Autoflow once per sample
        number_of_lanes=1
    fi
    s_num=$(ls $read_path -1 | grep $sample | grep -Eo "_S._" | cut -f 2 -d "_" | sort -u) # Extract sample number (Sn being n the number)
    for (( i = 1; i <= $number_of_lanes; i++ )) ; do
        AF_VARS=`echo "
        \\$sample=$sample,
        \\$SNUM=$s_num,
        \\$lane_num=$i
        " | tr -d [:space:]`
        AutoFlow -w $TEMPLATE -V "$AF_VARS" -o $RESULTS_FOLDER/$sample #$RESOURCES
    done
done < $SAMPLES_FILE