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
  * _Plasmodium vivax_
```bash
scripts/download_vivax.sh
```
```bash
scripts/download_vivax_annotations.sh
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
mkdir plasmodium
cd ../..
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

bams = list.files(path = "results/mapped_reads/colobus", pattern = "*.bam$", full.names=TRUE)
gtf.file = "genomes/combined.gtf"
fc = featureCounts(bams, annot.ext=gtf.file,
    isGTFAnnotationFile=TRUE,
    isPairedEnd=TRUE,
    nthreads=8,
    allowMultiOverlap=TRUE)

save.image("fc.Rdata")
```
  * As file is running, create percent_parasitemia table in excel (enter Colobus_Reads_Mapped):
 > Colobus_Reads_Mapped = "Successfully assigned alignments"
 
  * Do same for "unique" pathogen data 
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread")

library(Rsubread)

bams = list.files(path = "results/mapped_reads/plasmodium", pattern = "*.bam$", full.names=TRUE)
gtf.file = "genomes/combined.gtf"
fc = featureCounts(bams, annot.ext=gtf.file,
    isGTFAnnotationFile=TRUE,
    isPairedEnd=TRUE,
    nthreads=8,
    allowMultiOverlap=TRUE)

save.image("Plasfc.Rdata")
```
  * As file is running, create percent_parasitemia table in excel (enter Plasmodium_Reads_Mapped):
 > Plasmodium_Reads_Mapped = "Successfully assigned alignments"
 
  * Calculate Total_Reads:
 > Total_Reads = sum(Colobus_Reads_Mapped, Plasmodium_Reads_Mapped)
  * Calculate Percent Parasitemia:
 > Percent Parasitemia = Plasmodium_Reads_Mapped / Total_Reads
 
 ## Find DE genes associated with parasitemia
 > Necessary module(s): r/intel/3.6.0
 
 > Necessary R package(s): BiocManager, edgeR, GO.db, gprofiler2
```bash
module load r/intel/3.6.0
R
```
```R
options(stringsAsFactors=FALSE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("GO.db")

library(limma)
library(edgeR)

library(ggplot2)
library(ggrepel)

library(plyr)
library(xtable)

#library(gprofiler2)

library("GO.db")

# --- Load feature count data

load("fc.Rdata")

dim(fc$counts) # number of genes and individuals
# [1] 39942   116

# --- Load parasitemia info

parasitemia = read.table("percent_parasitemia.txt", header=TRUE)
parasitemia$prop.Plasmodium = parasitemia$Plasmodium_Reads_Mapped / parasitemia$Colobus_Reads_Mapped

# --- Normalize and remove low count genes

# Can add lib.sizes=c() to this to not recalculate
fc.dge  = DGEList(counts=fc$counts, genes=fc$annotation)

# Filter for genes with low counts across conditions

keep = rowSums(cpm(fc.dge)>5) >= 2

# Make reference table list
ref.list = data.frame(Colobus.id = names(keep[keep == TRUE]))

write.table(ref.list, file="uniqresults/reference_list.txt",
        sep="\t", row.names=FALSE)

fc.dge = fc.dge[keep, , keep.lib.sizes=FALSE]

# Normalize for library compositional bias
fc.dge.norm  = calcNormFactors(fc.dge)
fc.dge.norm$samples

# Get CPM just normalized by library composition
cpm  = cpm(fc.dge.norm)

# --- Estimate dispersion the complicated way (using CR method)

# Ensure samples are in correct order (i.e. make sure fcdge.norm$samples and parasitemia$Individual have same name)
table(gsub(".bam", "", row.names(fc.dge.norm$samples)) == parasitemia$Individual)
# TRUE
#   116

design = model.matrix(~ parasitemia$prop.Plasmodium)

disp = estimateDisp(fc.dge.norm, design, robust = TRUE)  # fc.dge.disp3  

# CPM
cpm.disp = cpm(disp)

fit = glmFit(disp, design, robust = TRUE)

lrt = glmLRT(fit, coef = 2)

tt = topTags(lrt, n=Inf, adjust.method = "BH", p.value = 0.05)

# Write list of significant up-regulated genes
sig.genes = tt$table$GeneID[tt$table$FDR < 1e-5 & tt$table$logFC > 0]
length(sig.genes)
# [1] 378

# Write background list of genes
bg.genes = tt$table$GeneID
length(bg.genes)
# [1] 3534

# Write results to file
write.table(sig.genes, file="reports/edgeR_colobus_upreg_sig_genes.txt",
    sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

write.table(bg.genes, file="reports/edgeR_colobus_bg_genes.txt",
    sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
```
