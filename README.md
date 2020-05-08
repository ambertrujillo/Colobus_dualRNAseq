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
sbatch/align_genome.sbatch
```
