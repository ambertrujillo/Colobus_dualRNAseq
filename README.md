# Colobus_dualRNAseq


# Pipeline

1. Download Reference Genomes and Annotation files from RefSeq assembly in ftp directory
> Necessary module(s): kent/328
  * _Piliocolobus tephrosceles_
```bash
scripts/download_Rcolobus.sh
```
```bash
scripts/download_Rcolobus_annotations.sh
```
  * _Plasmodium gonderi_
```bash
scripts/download_gonderi.sh
```
```bash
scripts/download_gonderi_annotations.sh
```
2. Concatenate Host-Pathogen Reference Genomes and Annotation files
  * Reference Genomes
```bash
cd genomes
cat *.fa > combined.fa
```
  * Annotation files
```bash
cat gonderi.gtf Rcolobus.fix.gtf > combined.gtf
```
3. Index concatenated Reference Genome
> Necessary module(s): gcc/9.1.0 and star/intel/2.5.2b
```sbatch
sbatch/index_genomes.sbatch
```
4. Download Ugandan Red Colobus Reads
> Necessary module(s): edirect/20181128, sra-tools/intel/2.9.6, and parallel/20171022
```bash
mkdir data
cd data
esearch -db sra -query PRJNA413051 | efetch --format runinfo | cut -d "," -f 1 > SRR.numbers
```
```
sbatch/prefetch.sbatch
```
```
sbatch/dump.sbatch
```
5. Trim Colobus reads
> Necessary module(s): trimmomatic/0.36
```bash
sbatch/trim_reads.sbatch
```
