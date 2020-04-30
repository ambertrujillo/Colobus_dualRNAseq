# Colobus_dualRNAseq


# Pipeline

1. Download Reference Genomes and Annotation files from RefSeq assembly in ftp directory
> Necessary module: kent/328
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
