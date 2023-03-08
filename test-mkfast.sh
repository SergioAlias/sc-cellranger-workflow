#!/bin/bash
# Sergio Alías, 20230308
# Last modified 20230308

# Test code for Cell Ranger mkfast

cellranger mkfastq --id=tiny-bcl \
                   --run=cellranger-tiny-bcl-1.2.0 \
                   --csv=cellranger-tiny-bcl-simple-1.2.0.csv