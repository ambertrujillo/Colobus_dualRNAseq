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

module load bamtools/intel/2.4.1
module load samtools/intel/1.9

DATA_DIR=$SCRATCH/colobus_hep/data
THIS_IND=`ls $DATA_DIR/*_1.fastq | sed -e "s:.*/::" | grep -vi TRIM | head -n $SLURM_ARRAY_TASK_ID | tail -n1 | sed -e "s/_1.*//"`
READ_DIR=$SCRATCH/colobus_hep/results

echo "Extracting primary reads..."

#Only output primary reads
samtools view -F 256 -b $READ_DIR/${THIS_IND}Aligned.sortedByCoord.out.bam > $READ_DIR/${THIS_IND}.primary.bam


echo "Splitting bam files..."

#Splits bam files into chromosomes
bamtools split -in $READ_DIR/${THIS_IND}.primary.bam -reference 

echo "Removing unmapped reads and unplaced scaffolds..."

#Remove unmapped and unplaced scaffolds
rm $READ_DIR/${THIS_IND}.primary.REF_NW_022*
rm $READ_DIR/${THIS_IND}.primary.REF_LR699* #LR chromosomes are mitochondrial and apicomplexin. Removed because they contain essentially nothing
rm $READ_DIR/${THIS_IND}.primary.REF_unmapped.bam

echo "Merging bams..."

#-->Colobus
ls $READ_DIR/${THIS_IND}.primary.REF_Colchr*.bam > $READ_DIR/${THIS_IND}.Colbamlist.txt
bamtools merge -list $READ_DIR/${THIS_IND}.Colbamlist.txt -out $READ_DIR/${THIS_IND}.colobus.bam

#-->Hepatocystis
ls $READ_DIR/${THIS_IND}.primary.REF_CABPSV*.bam > $READ_DIR/${THIS_IND}.hepatobamlist.txt
bamtools merge -list $READ_DIR/${THIS_IND}.hepatobamlist.txt -out $READ_DIR/${THIS_IND}.hepato.bam

echo "Removing mess you made ..."

#Remove mess you made
rm $READ_DIR/${THIS_IND}.primary.REF_Colchr*.bam
rm $READ_DIR/${THIS_IND}.primary.REF_CABPSV*.bam

echo "Moving unique BAMS..."

#Move unique colobus reads to mapped directory
mv $READ_DIR/${THIS_IND}.colobus.bam $READ_DIR/mapped_reads/colobus
mv $READ_DIR/${THIS_IND}.hepato.bam $READ_DIR/mapped_reads/hepatocystis


echo "DONE"
