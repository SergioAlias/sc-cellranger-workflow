#!/bin/bash
# Sergio Alías, 20230315
# Last modified 20230320

# Test code for Cell Ranger count

echo `which cellranger`

cellranger count --id=pbmc1k-results \
                 --transcriptome=refdata-gex-GRCh38-2020-A \
                 --fastqs=pbmc_1k_v3_fastqs/ \
                 --sample=pbmc_1k_v3 \
      		 --localcores=10 \
		 --localmem=60
