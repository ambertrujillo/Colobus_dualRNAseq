#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Download Plasmodium gonderi annotations
# ----------------------------------------------------------------------------------------

# --- Download annotations for Plasmodium gonderi genome, gonderi

cd papio_gonderi/genomes/gonderi #Think I should be in /home/aet359/baboon_malaria not papAnu3?

PGONDERI_URL=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/157/705

ANNO1=GCF_002157705.1_Pgonderi_assembly01/GCF_002157705.1_Pgonderi_assembly01_genomic.gtf
ANNO2=GCF_002157705.1_Pgonderi_assembly01/GCF_002157705.1_Pgonderi_assembly01_genomic.gff

for ANNO in $ANNO1 $ANNO2; do
    wget $PGONDERI_URL/$ANNO.gz \
        -O `basename $ANNO`.gz
    gunzip `basename $ANNO`.gz
done

#rename gtf file
mv GCF_002157705.1_Pgonderi_assembly01_genomic.gtf gonderi.gtf 
