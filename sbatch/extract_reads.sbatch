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

echo "Merge SRRs into individual colobus samples"
RC100.plasmodium.bam SRR6182879.plasmodium.bam SRR6182881.plasmodium.bam SRR6182897.plasmodium.bam SRR6182931.plasmodium.bam
RC101.plasmodium.bam SRR6182880.plasmodium.bam SRR6182882.plasmodium.bam SRR6182919.plasmodium.bam SRR6182928.plasmodium.bam
RC102.plasmodium.bam SRR6182877.plasmodium.bam SRR6182920.plasmodium.bam SRR6182927.plasmodium.bam SRR6182960.plasmodium.bam
RC103.plasmodium.bam SRR6182878.plasmodium.bam SRR6182921.plasmodium.bam SRR6182945.plasmodium.bam SRR6182959.plasmodium.bam
RC104.plasmodium.bam SRR6182875.plasmodium.bam SRR6182922.plasmodium.bam SRR6182946.plasmodium.bam SRR6182958.plasmodium.bam
RC105.plasmodium.bam SRR6182876.plasmodium.bam SRR6182923.plasmodium.bam SRR6182943.plasmodium.bam SRR6182957.plasmodium.bam
RC106R.plasmodium.bam SRR6182873.plasmodium.bam SRR6182924.plasmodium.bam SRR6182944.plasmodium.bam SRR6182964.plasmodium.bam
RC107.plasmodium.bam SRR6182874.plasmodium.bam SRR6182925.plasmodium.bam SRR6182941.plasmodium.bam SRR6182963.plasmodium.bam
RC108.plasmodium.bam SRR6182871.plasmodium.bam SRR6182926.plasmodium.bam SRR6182942.plasmodium.bam SRR6182962.plasmodium.bam
RC109.plasmodium.bam SRR6182872.plasmodium.bam SRR6182917.plasmodium.bam SRR6182939.plasmodium.bam SRR6182961.plasmodium.bam
RC110.plasmodium.bam SRR6182918.plasmodium.bam SRR6182940.plasmodium.bam SRR6182953.plasmodium.bam SRR6182966.plasmodium.bam
RC111.plasmodium.bam SRR6182937.plasmodium.bam SRR6182954.plasmodium.bam SRR6182965.plasmodium.bam SRR6182980.plasmodium.bam
RC112.plasmodium.bam SRR6182938.plasmodium.bam SRR6182951.plasmodium.bam SRR6182973.plasmodium.bam SRR6182979.plasmodium.bam
RC115.plasmodium.bam SRR6182914.plasmodium.bam SRR6182952.plasmodium.bam SRR6182974.plasmodium.bam SRR6182982.plasmodium.bam
RC116.plasmodium.bam SRR6182913.plasmodium.bam SRR6182949.plasmodium.bam SRR6182975.plasmodium.bam SRR6182981.plasmodium.bam
RC117.plasmodium.bam SRR6182916.plasmodium.bam SRR6182950.plasmodium.bam SRR6182976.plasmodium.bam SRR6182984.plasmodium.bam
RC118.plasmodium.bam SRR6182915.plasmodium.bam SRR6182947.plasmodium.bam SRR6182969.plasmodium.bam SRR6182983.plasmodium.bam
RC119.plasmodium.bam SRR6182910.plasmodium.bam SRR6182948.plasmodium.bam SRR6182970.plasmodium.bam SRR6182986.plasmodium.bam
RC120.plasmodium.bam SRR6182909.plasmodium.bam SRR6182955.plasmodium.bam SRR6182971.plasmodium.bam SRR6182985.plasmodium.bam
RC121.plasmodium.bam SRR6182912.plasmodium.bam SRR6182956.plasmodium.bam SRR6182972.plasmodium.bam SRR6182978.plasmodium.bam
RC122.plasmodium.bam SRR6182906.plasmodium.bam SRR6182911.plasmodium.bam SRR6182967.plasmodium.bam SRR6182977.plasmodium.bam
RC123.plasmodium.bam SRR6182885.plasmodium.bam SRR6182905.plasmodium.bam SRR6182908.plasmodium.bam SRR6182968.plasmodium.bam
RC124.plasmodium.bam SRR6182886.plasmodium.bam SRR6182904.plasmodium.bam SRR6182907.plasmodium.bam SRR6182934.plasmodium.bam
RC125.plasmodium.bam SRR6182883.plasmodium.bam SRR6182895.plasmodium.bam SRR6182903.plasmodium.bam SRR6182933.plasmodium.bam
RC126.plasmodium.bam SRR6182884.plasmodium.bam SRR6182896.plasmodium.bam SRR6182902.plasmodium.bam SRR6182936.plasmodium.bam
RC127.plasmodium.bam SRR6182889.plasmodium.bam SRR6182893.plasmodium.bam SRR6182901.plasmodium.bam SRR6182935.plasmodium.bam
RC128.plasmodium.bam SRR6182890.plasmodium.bam SRR6182894.plasmodium.bam SRR6182900.plasmodium.bam SRR6182930.plasmodium.bam
RC129.plasmodium.bam SRR6182887.plasmodium.bam SRR6182891.plasmodium.bam SRR6182899.plasmodium.bam SRR6182929.plasmodium.bam
RC130.plasmodium.bam SRR6182888.plasmodium.bam SRR6182892.plasmodium.bam SRR6182898.plasmodium.bam SRR6182932.plasmodium.bam




echo "DONE extracting"
