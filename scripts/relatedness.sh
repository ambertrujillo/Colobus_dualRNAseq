#!/bin/sh

############## Call colobus SNPs for relatedness assessment ###############

############## load packages ###############

module load jvarkit/base
module load picard/2.23.8
module load gatk/4.2.0.0
module load bamtools/intel/2.5.1
module load samtools/intel/1.12
module load gcc/10.2.0
module load star/intel/2.7.6a

########################################### BELOW needs to be done prior to calling SNPs for each individual ##########################################

############## Prepare Reference Genome ###############

# Index reference genome
samtools faidx genomes/Rcolobus.fa
# Creat dictionary file for reference genome
gatk CreateSequenceDictionary -R genomes/Rcolobus.fa

########################################### ABOVE needs to be done prior to calling SNPs for each individual ##########################################

# Call with: 
for IN_FQ in `find . -name "*_R1_001*_TRIM.fastq" -type f -printf '%f\n' | sed -e "s/_TRIM.fastq//" | sort | uniq`; do 
  echo "$IN_FQ"
  sh snp_call/scripts/call_snps.sh $IN_FQ
done

######################################
######### Mark Duplicates ############
######################################

IN_BAM=$1  # call with sh snp_call.sh RC100_S1_R1_001.colobus

java -jar $PICARD_JAR MarkDuplicates \
      -I snp_call/EMMREML_results/mapped_reads/colobus/$IN_BAM.colobus.bam \
      -O snp_call/${IN_BAM}_duplicates.bam \
      -M snp_call/${IN_BAM}_marked_dup_metrics.txt  
               
##############################################################
######### Split reads for pre-RNA unspliced product ##########
##############################################################

############## Split reads ###############
gatk SplitNCigarReads \
      -R genomes/Rcolobus.fa \
      -I snp_call/${IN_BAM}_duplicates.bam \
      -O snp_call/${IN_BAM}_split.bam
      
##############################################################
################## Base Call Recalibration ###################
##############################################################

# Generate recalibration table for Base Quality Score Recalibration (BQSR)
# Not doing this step because Colobus does not have a database of known polymorphisms for the machine learning algorithm to use. It would only work on humans and model organisms. 

######################################################
################## Variant Calling ###################
######################################################

# Inspect the new bam file for any errors
#java -jar $PICARD_JAR ValidateSamFile \
#      I=snp_call/RC100_S1_R1_001.colobus_split.fixed.bam \
#      MODE=SUMMARY
      
# Missing read group 1 error?
#java -jar $PICARD_JAR AddOrReplaceReadGroups \
#	INPUT=snp_call/${IN_BAM}_split.bam \
#	OUTPUT=snp_call/${IN_BAM}_split.fixed.bam \
#	RGLB=${IN_BAM} \
#	RGPL=Illumina \
#	RGPU=Group1 \
#	RGSM=${IN_BAM}
	
# Call germline SNPs and indels via local re-assembly of haplotypes
# Single-sample GVCF calling (outputs intermediate GVCF)
#gatk --java-options "-Xmx4g" HaplotypeCaller  \
#   -R Homo_sapiens_assembly38.fasta \
#   -I input.bam \
#   -O output.g.vcf.gz \
#   -ERC GVCF
   
# Single-sample GVCF calling with allele-specific annotations
#gatk --java-options "-Xmx4g" HaplotypeCaller  \
#   -R Homo_sapiens_assembly38.fasta \
#   -I input.bam \
#   -O output.g.vcf.gz \
#   -ERC GVCF \
#   -G Standard \
#   -G AS_Standard

# Make sure input bam is okay for snp calling
# java -jar $PICARD_JAR ValidateSamFile \
#      -I snp_call/${IN_BAM}_split.bam \
#      -MODE SUMMARY
      
# Odds are, you are missing read group 1. Fix with the following
# Missing read group 1 error?
java -jar $PICARD_JAR AddOrReplaceReadGroups \
	-I snp_call/${IN_BAM}_split.bam \
	-O snp_call/${IN_BAM}_split.fixed.bam \
	-RGLB ${IN_BAM} \
	-RGPL Illumina \
	-RGPU Group1 \
	-RGSM ${IN_BAM}

# Index fixed input bam
samtools index /scratch/aet359/colobus_hep/snp_call/${IN_BAM}_split.fixed.bam

# Variant calling with bamout to show realigned reads one chromosome at a time (for computing efficiency)

for CHR in `seq 1 21`; do
    gatk --java-options "-Xmx4G" HaplotypeCaller \
         -I snp_call/${IN_BAM}_split.fixed.bam \
         -O snp_call/${IN_BAM}_${CHR}_output.vcf.gz \
         -R genomes/Rcolobus.fa \
         -bamout snp_call/${IN_BAM}_${CHR}.bam \
         -L Colchr$CHR
done

cat snp_call/${IN_BAM}_*_output.vcf.gz > snp_call/${IN_BAM}_all.vcf.gz

########################################################
################## Variant Filtering ###################
########################################################
# Index variant file in order to enable quaries by interval (-window)
gatk IndexFeatureFile \
     -I snp_call/${IN_BAM}_all.vcf.gz
     
# Filter variants
gatk VariantFiltration \
    -R genomes/Rcolobus.fa \
    -V snp_call/${IN_BAM}_all.vcf.gz \
    -window 35 \
    -cluster 3 \
    --filter-name FS \
    -filter "FS > 30.0" \
    --filter-name QD \
    -filter "QD < 2.0" \
    -O snp_call/${IN_BAM}_filtered_output.vcf 

########################################### After all individuals have their SNPs called, run the following ##########################################

########################################################################################################################################
################## Merge VCFs because you have multiple individuals. Not simply trying to concatenate them together ####################
########################################################################################################################################

module load bamtools/intel/2.5.1
module load samtools/intel/1.12
module load gcc/10.2.0

# Compress and index vcf files
for IN_FQ in `ls *_R1_001_filtered_output.vcf`; do 
  echo "compressing $IN_FQ"
  bgzip $IN_FQ
done

for IN_FQ in `ls *_R1_001_filtered_output.vcf.gz`; do 
  echo "indexing $IN_FQ"
  tabix -p vcf $IN_FQ
done

# use vcf-merge and call with sbatch merge_vcf.sbatch
module load vcftools/intel/0.1.16

vcf-merge `ls snp_call/*_R1_001_filtered_output.vcf.gz` | bgzip -c > snp_call/relatedness/all_individuals_out.vcf.gz

# TEST
vcf-merge snp_call/RC130_S29_R1_001_filtered_output.vcf.gz snp_call/RC123_S22_R1_001_filtered_output.vcf.gz | bgzip -c > snp_call/relatedness/123_130_out.vcf.gz

####################################################
################## Call kinship ####################
####################################################

#mkdir -p lcmlkin
#lcmlkin -i snp_call/macaque_all.vcf \
#	-o lcmlkin/macaque_all.txt \
#	-l phred \
#	-g all 

module load plink2/20210103

# Make bed file
plink2 --vcf snp_call/relatedness/all_individuals_out.vcf.gz --make-bed --out snp_call/relatedness/colobus.relatedness --allow-extra-chr --max-alleles 2 --chr Colchr1-Colchr21

# Estimate relatedness using Plink
plink2 --vcf snp_call/relatedness/all_individuals_out.vcf.gz --distance square 1-ibs --out snp_call/relatedness/colobus.king.relatedness --const-fid 0 --allow-extra-chr 
# Estimate relatedness using King
plink2 --vcf snp_call/relatedness/all_individuals_out.vcf.gz --make-king square --out snp_call/relatedness/colobus.king.relatedness --const-fid 0 --allow-extra-chr 

# download KING (Doesn't seem to work because we used arbitrary chromosome names
wget https://www.kingrelatedness.com/Linux-king.tar.gz
tar -xzvf Linux-king.tar.gz


# Estimate kindship
./king -b snp_call/relatedness/colobus.relatedness.bed --related

