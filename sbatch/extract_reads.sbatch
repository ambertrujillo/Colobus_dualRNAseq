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

DATA_DIR=$SCRATCH/colobus_project/colobus_vivax/data/Rcolobus
THIS_IND=`ls $DATA_DIR/*_1.fastq | sed -e "s:.*/::" | grep -vi TRIM | head -n $SLURM_ARRAY_TASK_ID | tail -n1 | sed -e "s/_1.*//"`
READ_DIR=$SCRATCH/colobus_project/colobus_vivax/results

echo "Extracting primary reads..."

#Only output primary reads
samtools view -F 256 -b $READ_DIR/${THIS_IND}Aligned.sortedByCoord.out.bam > $READ_DIR/${THIS_IND}.primary.bam

echo "Splitting bam files..."

#Splits bam files into chromosomes
bamtools split -in $READ_DIR/${THIS_IND}.primary.bam -reference 

echo "Removing unmapped reads and unplaced scaffolds..."

#Remove unmapped and unplaced scaffolds
rm $READ_DIR/${THIS_IND}.primary.REF_NW_022*
rm $READ_DIR/${THIS_IND}.primary.REF_NW_0018*
rm $READ_DIR/${THIS_IND}.primary.REF_unmapped.bam

echo "Merging bams..."

#-->Colobus
bamtools merge -in $READ_DIR/${THIS_IND}.primary.REF_Colchr1.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Colchr2.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Colchr3.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Colchr4.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Colchr5.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Colchr6.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Colchr7.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Colchr8.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Colchr9.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Colchr10.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Colchr11.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Colchr12.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Colchr13.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Colchr14.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Colchr15.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Colchr16.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Colchr17.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Colchr18.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Colchr19.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Colchr20.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Colchr21.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_ColchrX.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_ColchrY.bam \
-out $READ_DIR/${THIS_IND}.colobus.bam

#-->Plasmodium
bamtools merge -in $READ_DIR/${THIS_IND}.primary.REF_Plaschr1.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Plaschr2.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Plaschr3.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Plaschr4.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Plaschr5.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Plaschr6.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Plaschr7.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Plaschr8.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Plaschr9.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Plaschr10.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Plaschr11.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Plaschr12.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Plaschr13.bam \
-in $READ_DIR/${THIS_IND}.primary.REF_Plaschr14.bam \
-out $READ_DIR/${THIS_IND}.plasmodium.bam

echo "Removing mess you made ..."

#Remove mess you made
rm $READ_DIR/${THIS_IND}.primary.REF_Colchr*.bam
rm $READ_DIR/${THIS_IND}.primary.REF_Plaschr*.bam

echo "Moving unique BAMS..."

#Move unique colobus reads to mapped directory
mv $READ_DIR/${THIS_IND}.colobus.bam $READ_DIR/mapped_reads/colobus
mv $READ_DIR/${THIS_IND}.plasmodium.bam $READ_DIR/mapped_reads/plasmodium


echo "DONE"
