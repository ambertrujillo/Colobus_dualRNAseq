#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=48GB
#SBATCH --job-name=trim_reads
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aet359@nyu.edu
#SBATCH --output=trim_%A_%a.out

# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

# to submit an array job put sbatch --array=1-[number of individuals]

working_dir=$SLURM_SUBMIT_DIR

# ------------------------------------------------------------------------------
# Trim sequences
# ------------------------------------------------------------------------------

cd $working_dir

DATA_DIR=/scratch/aet359/colobus_paper/SnakeMake/data
THIS_IND=`ls $DATA_DIR/*_001.fastq | sed -e "s:.*/::" | grep -vi TRIM | head -n $SLURM_ARRAY_TASK_ID | tail -n1 | sed -e "s/_001.*//"`

module load trimmomatic/0.36

ADAPTERS=$TRIMMOMATIC_ROOT/adapters/TruSeq2-SE.fa

echo "Trimming for individual $THIS_IND";

java -jar $TRIMMOMATIC_ROOT/trimmomatic-0.36.jar SE -phred33 \
    -threads 8 \
    $DATA_DIR/${THIS_IND}_001.fastq \
    $DATA_DIR/${THIS_IND}_TRIM_R1.fastq.gz \
    ILLUMINACLIP:$ADAPTERS:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
exit;

