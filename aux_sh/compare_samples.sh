#! /usr/bin/env bash

# Sergio Alias, 20230530
# Last modified 20230721

# STAGE 2 SAMPLES COMPARISON

#SBATCH -J compare_samples.sh
#SBATCH --cpus-per-task=3
#SBATCH --mem='5gb'
#SBATCH --constraint=cal
#SBATCH --time=0-01:00:00
#SBATCH --error=job.comp.%J.err
#SBATCH --output=job.comp.%J.out

# Setup

hostname

mkdir -p $report_folder


## TODO fix code below

while IFS= read -r name; do
  if [ -d "$COUNT_RESULTS_FOLDER/$name" ]; then
    input_file="$COUNT_RESULTS_FOLDER/$name/cellranger_0000/$name/outs/metrics_summary.csv"
    cat $input_file | perl -pe 's/(\d),(\d)/$1$2/g'| sed '1 s/ /_/g' | sed 's/%//g' | sed 's/"//g' | sed 's/ /\n/g' | sed 's/,/\t/g' | awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' | awk -v var="$name" 'BEGIN {FS=OFS="\t"} {print var, $0}' | sed 's/ /\t/g' >> $experiment_folder'/cellranger_metrics'
  fi
done < $SAMPLES_FILE


. ~soft_bio_267/initializes/init_ruby
. ~soft_bio_267/initializes/init_R
create_metric_table.rb $experiment_folder'/metrics' sample $experiment_folder'/metric_table'
create_metric_table.rb $experiment_folder'/cellranger_metrics' sample $experiment_folder'/cellranger_metric_table'

# Main

/usr/bin/time $CODE_PATH/scripts/compare_samples.R -o $report_folder \
                                                   -m $experiment_folder'/metric_table' \
                                                   -l $experiment_folder'/metrics' \
                                                   -e $experiment_name \
                                                   --cellranger_metrics $experiment_folder'/cellranger_metric_table' \
                                                   --cellranger_long_metrics $experiment_folder'/cellranger_metrics'