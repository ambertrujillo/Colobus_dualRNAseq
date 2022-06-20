#!/bin/bash

module load gcc/10.2.0 star/intel/2.7.6a

DATA_DIR=$SCRATCH/colobus_paper/SnakeMake/data
IN_BAM=$1

READ_1="$DATA_DIR/${IN_BAM}"

STAR --genomeDir $SCRATCH/colobus_paper/SnakeMake/genomes/index \
--runThreadN 15 \
--readFilesIn $READ_1 \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 40712737493 \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterMismatchNoverReadLmax 0.04 \
--twopassMode Basic \
--outFileNamePrefix $SCRATCH/colobus_paper/SnakeMake/EMMREML_results/${IN_BAM}
