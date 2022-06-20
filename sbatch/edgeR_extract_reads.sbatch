#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=124GB
#SBATCH --job-name=extract_reads
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aet359@nyu.edu
#SBATCH --output=extract_%A_%a.out

# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

working_dir=$SLURM_SUBMIT_DIR

cd $working_dir

module load bamtools/intel/2.5.1
module load samtools/intel/1.12

DATA_DIR=/scratch/aet359/colobus_paper/SnakeMake/data
THIS_IND=`ls $DATA_DIR/*_001.fastq | sed -e "s:.*/::" | grep -vi TRIM | head -n $SLURM_ARRAY_TASK_ID | tail -n1 | sed -e "s/_001.*//"`
READ_DIR=/scratch/aet359/colobus_paper/SnakeMake/edgeR_results

echo "Extracting primary reads and reads from desired chromosomes..."

#Only output primary reads
samtools view -F 256 -b $READ_DIR/${THIS_IND}Aligned.sortedByCoord.out.bam > $READ_DIR/${THIS_IND}.primary.bam

samtools view -b -L /scratch/aet359/colobus_paper/SnakeMake/genomes/to_include.bed $READ_DIR/${THIS_IND}.primary.bam > $READ_DIR/${THIS_IND}.colobus_hepato.bam

echo "Splitting bam files..."

bamtools split -in $READ_DIR/${THIS_IND}.colobus_hepato.bam -reference

echo "Removing unnecessary mitochondrial and apicomplexin chromosomes..."

#Remove unmapped and unplaced scaffolds -- No need for this
rm $READ_DIR/${THIS_IND}.colobus_hepato.REF_LR699* #LR chromosomes are mitochondrial and apicomplexin. Removed because they contain essentially nothing

echo "Merging bams..."

#-->Colobus
ls $READ_DIR/${THIS_IND}.colobus_hepato.REF_Colchr*.bam > $READ_DIR/${THIS_IND}.Colbamlist.txt
bamtools merge -list $READ_DIR/${THIS_IND}.Colbamlist.txt -out $READ_DIR/${THIS_IND}.colobus.bam

#-->Hepatocystis
ls $READ_DIR/${THIS_IND}.colobus_hepato.REF_CABPSV020*.bam > $READ_DIR/${THIS_IND}.hepatobamlist.txt
bamtools merge -list $READ_DIR/${THIS_IND}.hepatobamlist.txt -out $READ_DIR/${THIS_IND}.hepato.bam

echo "Moving unique BAMS..."

#Move unique colobus reads to mapped directory
mv $READ_DIR/${THIS_IND}.colobus.bam $READ_DIR/mapped_reads/colobus
mv $READ_DIR/${THIS_IND}.hepato.bam $READ_DIR/mapped_reads/hepatocystis

echo "Removing mess you made ..."

#Remove mess you made
rm $READ_DIR/${THIS_IND}.colobus_hepato.REF_Colchr*.bam
rm $READ_DIR/${THIS_IND}.colobus_hepato.REF_CABPSV*.bam


echo "DONE"
