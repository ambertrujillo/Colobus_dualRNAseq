#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=48GB
#SBATCH --job-name=staralign_reads
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aet359@nyu.edu
#SBATCH --output=align_%A_%a.out


# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

working_dir=$SLURM_SUBMIT_DIR

cd $working_dir

module load gcc/9.1.0 star/intel/2.5.2b

gunzip $SCRATCH/colobus_project/STAR/data/Rcolobus/*TRIM_*.fastq.gz

DATA_DIR=$SCRATCH/colobus_project/STAR/data/Rcolobus
THIS_IND=`ls $DATA_DIR/*_1.fastq | sed -e "s:.*/::" | grep -vi TRIM | head -n $SLURM_ARRAY_TASK_ID | tail -n1 | sed -e "s/_1.*//"`


READ_1="$DATA_DIR/${THIS_IND}_TRIM_R1.fastq"
READ_2="$DATA_DIR/${THIS_IND}_TRIM_R2.fastq"


STAR --genomeDir $SCRATCH/colobus_project/STAR/genomes/index/ \
--runThreadN 15 \
--readFilesIn $READ_1 $READ_2 \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 40712737493 \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterMismatchNoverReadLmax 0.04 \
--outFileNamePrefix $SCRATCH/colobus_project/STAR/results/${THIS_IND}

