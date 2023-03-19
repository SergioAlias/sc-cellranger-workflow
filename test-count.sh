#!/bin/bash
# Sergio Alías, 20230315
# Last modified 20230319

# Test code for Cell Ranger count

echo `which cellranger`

cellranger count --id=tiny-bcl \
                 --transcriptome=refdata-gex-GRCh38-2020-A \
                 --fastqs=../testdata/100k_reads_hiseq/TESTX \
                 --sample=TESTX_H7YRLADXX
