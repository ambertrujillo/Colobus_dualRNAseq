# Colobus_dualRNAseq


# Pipeline
## Align RNAseq Data to Reference Genome

1. Download Reference Genomes and Annotation files from RefSeq assembly in ftp directory
> Necessary module(s): kent/328
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
cat *.gtf > combined.gtf
```
3. Index concatenated Reference Genome (SBATCH JOB)
> Necessary module(s): gcc/9.1.0 and star/intel/2.5.2b

> To submit sbatch job: `sbatch sbatch/name_of_job.sbatch`

```sbatch
sbatch/index_genomes.sbatch
```
4. Download Ugandan Red Colobus Reads
> Necessary module(s): edirect/20181128, sra-tools/intel/2.9.6, and parallel/20171022

> After creating SRR.numbers, prefetch and dump are SBATCH JOBS

```bash
mkdir data/Rcolobus
cd data/Rcolobus
esearch -db sra -query PRJNA413051 | efetch --format runinfo | cut -d "," -f 1 > SRR.numbers
```
```
sbatch/prefetch.sbatch
```
```
sbatch/dump.sbatch
```
5. Trim Colobus reads (ARRAY JOB)
> Necessary module(s): trimmomatic/0.36

> To submit an sbatch array job: `sbatch --array=1-[number of individuals] sbatch/name_of_job.sbatch`

```bash
sbatch/trim_reads.sbatch
```
6. Align Colobus reads to concatenated Host-Pathogen Referance sequence (ARRAY JOB)
> Necessary module(s): gcc/9.1.0 and star/intel/2.5.2b
```bash
mkdir results
```
```bash
sbatch/align_genome.sbatch
```
## Obtain "Unique" Host and Pathogen Data (ARRAY JOB)
```bash
mkdir results/mapped_reads
cd results/mapped_reads
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
sbatch/extract_reads.sbatch
```

## Obtain Read Count Matrix and Calculate Percent Parasitemia
> Necessary module(s): r/intel/3.6.0

> Necessary R package(s): BiocManager, Rsubread
```bash
module load r/intel/3.6.0
R
```
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread")

library(Rsubread)

bams = list.files(path = "results/mapped_reads/colobus/merged_reads", pattern = "*.bam$", full.names=TRUE)
gtf.file = "genomes/combined.gtf"
colobusfc = featureCounts(bams, annot.ext=gtf.file,
    isGTFAnnotationFile=TRUE,
    isPairedEnd=TRUE,
    nthreads=8,
    allowMultiOverlap=TRUE)

save.image("colobusfc.Rdata")
```
  * As file is running, create percent_parasitemia table in excel (enter Colobus_Reads_Mapped):
 > Colobus_Reads_Mapped = "Successfully assigned alignments"
 
  * Do same for "unique" pathogen data 
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread")

library(Rsubread)

bams = list.files(path = "results/mapped_reads/plasmodium/merged_reads", pattern = "*.bam$", full.names=TRUE)
gtf.file = "genomes/combined.gtf"
Hepatofc = featureCounts(bams, annot.ext=gtf.file,
    GTF.featureType="gene",
    isGTFAnnotationFile=TRUE,
    isPairedEnd=TRUE,
    nthreads=8,
    allowMultiOverlap=FALSE)

save.image("Hepatofc.Rdata")
```
  * As file is running, create percent_parasitemia table in excel (enter Hepatocystis_Reads_Mapped):
 > Plasmodium_Reads_Mapped = "Successfully assigned alignments"
 
  * Calculate Total_Reads:
 > Total_Reads = sum(Colobus_Reads_Mapped, Hepatocystis_Reads_Mapped)
  * Calculate Percent Parasitemia:
 > Percent Parasitemia = Hepatocystis_Reads_Mapped / Total_Reads
 
 ## Find DE genes associated with parasitemia
 > Necessary module(s): r/intel/3.6.0
 
 > Necessary R package(s): BiocManager, edgeR, GO.db, gprofiler2
```bash
module load r/intel/3.6.0
Rscript scripts/find_DE_host_genes.R
```
 ## Single-Cell RNAseq Deconvolution
 > Necessary module(s): r/intel/3.6.0
 
 > Necessary R package(s): Seurat, limma, Biobase, reshape2, bseqsc, xbioc, ggplot2, tidyr
 
 > Necessary Files: phenotype_data.csv, single_cell_gene_matrix.csv, genehepato.ind.fc.Rdata (all available in data_files)
 ```bash
 module load r/intel/3.6.0
 Rscript script/deconvolution.R
 ```
## Correlate _Hepatocystis sp._ Life Stage and Colobus Immune Cell Type
> Necessary R package(s): Hmisc
```bash
module load r/intel/3.6.0
Rscript script/explore_hepato_life_stage_data.R
```
