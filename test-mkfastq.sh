#!/bin/bash
# Sergio Alías, 20230308
# Last modified 20230308

# Test code for Cell Ranger mkfast
# rc-i2-override for solving the I2 issue (version < 7.0.1)

# export PATH=~sergioalias/cellranger-7.1.0:$PATH
echo `which cellranger`


cellranger mkfastq --id=tiny-bcl-7-0-0 \
                   --run=cellranger-tiny-bcl-1.2.0 \
                   --csv=cellranger-tiny-bcl-simple-1.2.0.csv #\
                   # --rc-i2-override=true
                   # --samplesheet=cellranger-tiny-bcl-samplesheet-1.2.0.csv