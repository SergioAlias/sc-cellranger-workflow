#!/bin/bash
# Sergio Alías, 20230308
# Last modified 20230315

# Test code for Cell Ranger mkfast

# export PATH=~sergioalias/cellranger-7.1.0:$PATH
echo `which cellranger`


cellranger mkfastq --id=tiny-bcl \
                   --run=cellranger-tiny-bcl-1.2.0 \
                   --csv=cellranger-tiny-bcl-simple-1.2.0.csv #\
                   # --rc-i2-override=true
                   # --samplesheet=cellranger-tiny-bcl-samplesheet-1.2.0.csv