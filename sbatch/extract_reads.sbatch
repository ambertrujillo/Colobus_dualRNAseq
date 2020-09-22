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

echo "Merge SRRs into individual colobus samples...Hepatocystis"

samtools merge RC100.hepato.bam SRR6182879.hepato.bam SRR6182881.hepato.bam SRR6182897.hepato.bam SRR6182931.hepato.bam
samtools merge RC101.hepato.bam SRR6182880.hepato.bam SRR6182882.hepato.bam SRR6182919.hepato.bam SRR6182928.hepato.bam
samtools merge RC102.hepato.bam SRR6182877.hepato.bam SRR6182920.hepato.bam SRR6182927.hepato.bam SRR6182960.hepato.bam
samtools merge RC103.hepato.bam SRR6182878.hepato.bam SRR6182921.hepato.bam SRR6182945.hepato.bam SRR6182959.hepato.bam
samtools merge RC104.hepato.bam SRR6182875.hepato.bam SRR6182922.hepato.bam SRR6182946.hepato.bam SRR6182958.hepato.bam
samtools merge RC105.hepato.bam SRR6182876.hepato.bam SRR6182923.hepato.bam SRR6182943.hepato.bam SRR6182957.hepato.bam
samtools merge RC106R.hepato.bam SRR6182873.hepato.bam SRR6182924.hepato.bam SRR6182944.hepato.bam SRR6182964.hepato.bam
samtools merge RC107.hepato.bam SRR6182874.hepato.bam SRR6182925.hepato.bam SRR6182941.hepato.bam SRR6182963.hepato.bam
samtools merge RC108.hepato.bam SRR6182871.hepato.bam SRR6182926.hepato.bam SRR6182942.hepato.bam SRR6182962.hepato.bam
samtools merge RC109.hepato.bam SRR6182872.hepato.bam SRR6182917.hepato.bam SRR6182939.hepato.bam SRR6182961.hepato.bam
samtools merge RC110.hepato.bam SRR6182918.hepato.bam SRR6182940.hepato.bam SRR6182953.hepato.bam SRR6182966.hepato.bam
samtools merge RC111.hepato.bam SRR6182937.hepato.bam SRR6182954.hepato.bam SRR6182965.hepato.bam SRR6182980.hepato.bam
samtools merge RC112.hepato.bam SRR6182938.hepato.bam SRR6182951.hepato.bam SRR6182973.hepato.bam SRR6182979.hepato.bam
samtools merge RC115.hepato.bam SRR6182914.hepato.bam SRR6182952.hepato.bam SRR6182974.hepato.bam SRR6182982.hepato.bam
samtools merge RC116.hepato.bam SRR6182913.hepato.bam SRR6182949.hepato.bam SRR6182975.hepato.bam SRR6182981.hepato.bam
samtools merge RC117.hepato.bam SRR6182916.hepato.bam SRR6182950.hepato.bam SRR6182976.hepato.bam SRR6182984.hepato.bam
samtools merge RC118.hepato.bam SRR6182915.hepato.bam SRR6182947.hepato.bam SRR6182969.hepato.bam SRR6182983.hepato.bam
samtools merge RC119.hepato.bam SRR6182910.hepato.bam SRR6182948.hepato.bam SRR6182970.hepato.bam SRR6182986.hepato.bam
samtools merge RC120.hepato.bam SRR6182909.hepato.bam SRR6182955.hepato.bam SRR6182971.hepato.bam SRR6182985.hepato.bam
samtools merge RC121.hepato.bam SRR6182912.hepato.bam SRR6182956.hepato.bam SRR6182972.hepato.bam SRR6182978.hepato.bam
samtools merge RC122.hepato.bam SRR6182906.hepato.bam SRR6182911.hepato.bam SRR6182967.hepato.bam SRR6182977.hepato.bam
samtools merge RC123.hepato.bam SRR6182885.hepato.bam SRR6182905.hepato.bam SRR6182908.hepato.bam SRR6182968.hepato.bam
samtools merge RC124.hepato.bam SRR6182886.hepato.bam SRR6182904.hepato.bam SRR6182907.hepato.bam SRR6182934.hepato.bam
samtools merge RC125.hepato.bam SRR6182883.hepato.bam SRR6182895.hepato.bam SRR6182903.hepato.bam SRR6182933.hepato.bam
samtools merge RC126.hepato.bam SRR6182884.hepato.bam SRR6182896.hepato.bam SRR6182902.hepato.bam SRR6182936.hepato.bam
samtools merge RC127.hepato.bam SRR6182889.hepato.bam SRR6182893.hepato.bam SRR6182901.hepato.bam SRR6182935.hepato.bam
samtools merge RC128.hepato.bam SRR6182890.hepato.bam SRR6182894.hepato.bam SRR6182900.hepato.bam SRR6182930.hepato.bam
samtools merge RC129.hepato.bam SRR6182887.hepato.bam SRR6182891.hepato.bam SRR6182899.hepato.bam SRR6182929.hepato.bam
samtools merge RC130.hepato.bam SRR6182888.hepato.bam SRR6182892.hepato.bam SRR6182898.hepato.bam SRR6182932.hepato.bam

mv RC* $READ_DIR/mapped_reads/hepatocystis/merged_reads

echo "Merge SRRs into individual colobus samples...Colobus"

samtools merge RC100.colobus.bam SRR6182879.colobus.bam SRR6182881.colobus.bam SRR6182897.colobus.bam SRR6182931.colobus.bam
samtools merge RC101.colobus.bam SRR6182880.colobus.bam SRR6182882.colobus.bam SRR6182919.colobus.bam SRR6182928.colobus.bam
samtools merge RC102.colobus.bam SRR6182877.colobus.bam SRR6182920.colobus.bam SRR6182927.colobus.bam SRR6182960.colobus.bam
samtools merge RC103.colobus.bam SRR6182878.colobus.bam SRR6182921.colobus.bam SRR6182945.colobus.bam SRR6182959.colobus.bam
samtools merge RC104.colobus.bam SRR6182875.colobus.bam SRR6182922.colobus.bam SRR6182946.colobus.bam SRR6182958.colobus.bam
samtools merge RC105.colobus.bam SRR6182876.colobus.bam SRR6182923.colobus.bam SRR6182943.colobus.bam SRR6182957.colobus.bam
samtools merge RC106R.colobus.bam SRR6182873.colobus.bam SRR6182924.colobus.bam SRR6182944.colobus.bam SRR6182964.colobus.bam
samtools merge RC107.colobus.bam SRR6182874.colobus.bam SRR6182925.colobus.bam SRR6182941.colobus.bam SRR6182963.colobus.bam
samtools merge RC108.colobus.bam SRR6182871.colobus.bam SRR6182926.colobus.bam SRR6182942.colobus.bam SRR6182962.colobus.bam
samtools merge RC109.colobus.bam SRR6182872.colobus.bam SRR6182917.colobus.bam SRR6182939.colobus.bam SRR6182961.colobus.bam
samtools merge RC110.colobus.bam SRR6182918.colobus.bam SRR6182940.colobus.bam SRR6182953.colobus.bam SRR6182966.colobus.bam
samtools merge RC111.colobus.bam SRR6182937.colobus.bam SRR6182954.colobus.bam SRR6182965.colobus.bam SRR6182980.colobus.bam
samtools merge RC112.colobus.bam SRR6182938.colobus.bam SRR6182951.colobus.bam SRR6182973.colobus.bam SRR6182979.colobus.bam
samtools merge RC115.colobus.bam SRR6182914.colobus.bam SRR6182952.colobus.bam SRR6182974.colobus.bam SRR6182982.colobus.bam
samtools merge RC116.colobus.bam SRR6182913.colobus.bam SRR6182949.colobus.bam SRR6182975.colobus.bam SRR6182981.colobus.bam
samtools merge RC117.colobus.bam SRR6182916.colobus.bam SRR6182950.colobus.bam SRR6182976.colobus.bam SRR6182984.colobus.bam
samtools merge RC118.colobus.bam SRR6182915.colobus.bam SRR6182947.colobus.bam SRR6182969.colobus.bam SRR6182983.colobus.bam
samtools merge RC119.colobus.bam SRR6182910.colobus.bam SRR6182948.colobus.bam SRR6182970.colobus.bam SRR6182986.colobus.bam
samtools merge RC120.colobus.bam SRR6182909.colobus.bam SRR6182955.colobus.bam SRR6182971.colobus.bam SRR6182985.colobus.bam
samtools merge RC121.colobus.bam SRR6182912.colobus.bam SRR6182956.colobus.bam SRR6182972.colobus.bam SRR6182978.colobus.bam
samtools merge RC122.colobus.bam SRR6182906.colobus.bam SRR6182911.colobus.bam SRR6182967.colobus.bam SRR6182977.colobus.bam
samtools merge RC123.colobus.bam SRR6182885.colobus.bam SRR6182905.colobus.bam SRR6182908.colobus.bam SRR6182968.colobus.bam
samtools merge RC124.colobus.bam SRR6182886.colobus.bam SRR6182904.colobus.bam SRR6182907.colobus.bam SRR6182934.colobus.bam
samtools merge RC125.colobus.bam SRR6182883.colobus.bam SRR6182895.colobus.bam SRR6182903.colobus.bam SRR6182933.colobus.bam
samtools merge RC126.colobus.bam SRR6182884.colobus.bam SRR6182896.colobus.bam SRR6182902.colobus.bam SRR6182936.colobus.bam
samtools merge RC127.colobus.bam SRR6182889.colobus.bam SRR6182893.colobus.bam SRR6182901.colobus.bam SRR6182935.colobus.bam
samtools merge RC128.colobus.bam SRR6182890.colobus.bam SRR6182894.colobus.bam SRR6182900.colobus.bam SRR6182930.colobus.bam
samtools merge RC129.colobus.bam SRR6182887.colobus.bam SRR6182891.colobus.bam SRR6182899.colobus.bam SRR6182929.colobus.bam
samtools merge RC130.colobus.bam SRR6182888.colobus.bam SRR6182892.colobus.bam SRR6182898.colobus.bam SRR6182932.colobus.bam

mv RC* $READ_DIR/mapped_reads/colobus/merged_reads

echo "DONE extracting"
