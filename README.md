# Colobus_dualRNAseq


# Pipeline
## Prepare Reference Genome and Colobus Reads for Analysis

1. Download Reference Genomes and Annotation files from RefSeq assembly in ftp directory
> Necessary module(s): kent/385
  * _Piliocolobus tephrosceles_
```bash
scripts/download_Rcolobus.sh
```
```bash
scripts/download_Rcolobus_annotations.sh
```
  * _Hepatocystis sp._
```bash
scripts/download_hep.sh
```
```bash
scripts/download_hep_annotations.sh
```
2. Concatenate Host-Pathogen Reference Genomes and Annotation files
  * Reference Genomes
```bash
cd genomes
cat *.fa > combined.fa
```
  * Annotation files
```bash
cat Rcolobus.fix.gtf hep.gtf > combined.gtf
```
3. Index concatenated Reference Genome (SBATCH JOB)
> Necessary module(s): gcc/10.2.0, star/intel/2.7.6a

> To submit sbatch job: `sbatch sbatch/name_of_job.sbatch`

```sbatch
sbatch/index_genomes.sbatch
```
4. Trim Colobus reads (ARRAY JOB)
> Necessary module(s): trimmomatic/0.36

> To submit an sbatch array job: `sbatch --array=1-[number of individuals] sbatch/name_of_job.sbatch`

```bash
sbatch/trim_reads.sbatch
```
## 5a. EdgeR pipeline
### Align Colobus reads to concatenated Host-Pathogen Referance sequence (ARRAY JOB)
> Necessary module(s): gcc/10.2.0 star/intel/2.7.6a
```bash
mkdir edgeR_results
gunzip $SCRATCH/colobus_hep/data/*TRIM_*.fastq.gz
```
```bash
sbatch/edgeR_align_genome.sbatch
```
### Obtain "Unique" Host and Pathogen Data (ARRAY JOB)
> Necessary module(s): bamtools/intel/2.5.1, samtools/intel/1.12
```bash
# Create index of reference genome
cd /scratch/aet359/colobus_paper/SnakeMake/genomes/
samtools faidx combined.fa
# Create bed file of chromosomes of interest for unique host and pathogen data
awk 'BEGIN {OFS="\t"}; { print $1,1,$2 }' combined.fa.fai | grep -v "NW_022" > to_include.bed
cd ..
# Create necessary results files
mkdir edgeR_results/mapped_reads
cd edgeR_results/mapped_reads
mkdir colobus
cd colobus
mkdir merged_reads
cd ..
mkdir hepatocystis
cd hepatocystis
mkdir merged_reads
cd ../../..
```
```bash
sbatch/edgeR_extract_reads.sbatch
```

### Obtain Read Count Matrix and Calculate Percent Parasitemia
> Necessary module(s): r/intel/4.0.4

> Necessary R package(s): BiocManager, Rsubread
```bash
module load r/intel/4.0.4
R
```
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread")

library(Rsubread)

bams = list.files(path = "edgeR_results/mapped_reads/colobus", pattern = "*.bam$", full.names=TRUE)
gtf.file = "genomes/combined.gtf"
colobus.fc = featureCounts(bams, annot.ext=gtf.file,
    GTF.featureType="gene",
    isGTFAnnotationFile=TRUE,
    isPairedEnd=FALSE,
    nthreads=8,
    allowMultiOverlap=FALSE)

save.image("colobus.fc.Rdata")
```
  * As file is running, create percent_parasitemia table in excel (enter Colobus_Reads_Mapped):
 > Colobus_Reads_Mapped = "Successfully assigned alignments"
 
  * Do same for "unique" pathogen data 
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread")

library(Rsubread)

bams = list.files(path = "edgeR_results/mapped_reads/hepatocystis", pattern = "*.bam$", full.names=TRUE)
gtf.file = "genomes/combined.gtf"
hepato.fc = featureCounts(bams, annot.ext=gtf.file,
    GTF.featureType="gene",
    isGTFAnnotationFile=TRUE,
    isPairedEnd=FALSE,
    nthreads=8,
    allowMultiOverlap=FALSE)

save.image("hepato.fc.Rdata")
```
  * As file is running, create percent_parasitemia table in excel (enter Hepatocystis_Reads_Mapped):
 > Plasmodium_Reads_Mapped = "Successfully assigned alignments"
 
  * Calculate Total_Reads:
 > Total_Reads = sum(Colobus_Reads_Mapped, Hepatocystis_Reads_Mapped)
  * Calculate Percent Parasitemia:
 > Percent Parasitemia = Hepatocystis_Reads_Mapped / Total_Reads
 ### Run Immune cell composition PCA for M<sub>immune</sub> model
 ```bash
 module load r/intel/4.0.4
 Rscript scripts/immune_cell_PCR.R
 ```
 ### Find DE genes associated with parasitemia
 > Necessary module(s):  r/intel/4.0.4
 
 > Necessary R package(s): BiocManager, edgeR, GO.db, gprofiler2
```bash
module load r/intel/4.0.4
# Colobus
Rscript scripts/edgeR_find_DE_host_genes.R
# Hepatocystis
Rscript scripts/edgeR_find_DE_pathogen_genes.R
```
## 5b. EMMREML pipeline
### Align Colobus reads to concatenated Host-Pathogen Referance sequence (ARRAY JOB)
 > Necessary module(s): gcc/10.2.0, star/intel/2.7.6a
```bash
mkdir EMMREML_results
```
```bash
sbatch/EMMREML_align_genome.sbatch
```
### Obtain "Unique" Host and Pathogen Data (ARRAY JOB)
> Necessary module(s): bamtools/intel/2.5.1, samtools/intel/1.12
```bash
# Create necessary results files
mkdir EMMREML_results/mapped_reads
cd EMMREML_results/mapped_reads
mkdir colobus
cd colobus
mkdir merged_reads
cd ..
mkdir hepatocystis
cd hepatocystis
mkdir merged_reads
cd ../../..
```
```bash
sbatch/EMMREML_extract_reads.sbatch
```

### Obtain Read Count Matrix and Calculate Percent Parasitemia
> Necessary module(s): r/intel/4.0.4

> Necessary R package(s): BiocManager, Rsubread
```bash
module load r/intel/4.0.4
R
```
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread")

library(Rsubread)

bams = list.files(path = "EMMREML_results/mapped_reads/colobus", pattern = "*.bam$", full.names=TRUE)
gtf.file = "genomes/combined.gtf"
colobus.EMMREML.fc = featureCounts(bams, annot.ext=gtf.file,
    GTF.featureType="gene",
    isGTFAnnotationFile=TRUE,
    isPairedEnd=FALSE,
    nthreads=8,
    allowMultiOverlap=FALSE)

save.image("colobus.EMMREML.fc.Rdata")
```
  * As file is running, create percent_parasitemia table in excel (enter Colobus_Reads_Mapped):
 > Colobus_Reads_Mapped = "Successfully assigned alignments"
 
  * Do same for "unique" pathogen data 
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread")

library(Rsubread)

bams = list.files(path = "EMMREML_results/mapped_reads/hepatocystis", pattern = "*.bam$", full.names=TRUE)
gtf.file = "genomes/combined.gtf"
hepato.EMMREML.fc = featureCounts(bams, annot.ext=gtf.file,
    GTF.featureType="gene",
    isGTFAnnotationFile=TRUE,
    isPairedEnd=FALSE,
    nthreads=8,
    allowMultiOverlap=FALSE)

save.image("hepato.EMMREML.fc.Rdata")
```
  * As file is running, create percent_parasitemia table in excel (enter Hepatocystis_Reads_Mapped):
 > Plasmodium_Reads_Mapped = "Successfully assigned alignments"
 
  * Calculate Total_Reads:
 > Total_Reads = sum(Colobus_Reads_Mapped, Hepatocystis_Reads_Mapped)
  * Calculate Percent Parasitemia:
 > Percent Parasitemia = Hepatocystis_Reads_Mapped / Total_Reads
 ### Calculate Relatedness
 > Necessary module(s): jvarkit/base, picard/2.23.8, gatk/4.2.0.0, bamtools/intel/2.5.1, samtools/intel/1.12, gcc/10.2.0, star/intel/2.7.6a
 ```bash
 module load r/intel/4.0.4
 Rscript scripts/relatedness.sh
 ```
 ### Find DE genes associated with parasitemia and relatedness
 > Necessary module(s):  r/intel/4.0.4
 
 > Necessary R package(s): BiocManager, edgeR, GO.db, gprofiler2
```bash
module load r/intel/4.0.4
Rscript scripts/EMMREML_find_DE_host_genes.R
```
### Correlate _Hepatocystis sp._ Life Stage and Colobus Immune Cell Type
> Necessary R package(s): Hmisc
```bash
module load r/intel/3.6.0
Rscript script/explore_hepato_life_stage_data.R
```
