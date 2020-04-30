#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --mem=48GB
#SBATCH --job-name=index_genome
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aet359@nyu.edu
#SBATCH --output=starindex_%j.out


# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

working_dir=$SLURM_SUBMIT_DIR

cd $working_dir

module load gcc/9.1.0 star/intel/2.5.2b

mkdir -p genomes/index

STAR --runThreadN 2 \
--runMode genomeGenerate \
--genomeDir genomes/index \
--genomeFastaFiles genomes/combined.fa \
--sjdbGTFfile genomes/combined.gtf \
--sjdbOverhang 99 \
--outFileNamePrefix genomes/index/combinedindex \
--limitGenomeGenerateRAM=40712737493

cd ..
