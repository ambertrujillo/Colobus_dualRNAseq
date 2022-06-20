#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Download Hepatocystis annotations
# ----------------------------------------------------------------------------------------

# --- Download annotations for Hepatocystis, Hep

cd genomes/ 

HEPATOCYSTIS_URL=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/902/459/845

ANNO1=GCA_902459845.2_HEP1/GCA_902459845.2_HEP1_genomic.gtf
ANNO2=GCA_902459845.2_HEP1/GCA_902459845.2_HEP1_genomic.gff

for ANNO in $ANNO1 $ANNO2; do
    wget $HEPATOCYSTIS_URL/$ANNO.gz \
        -O `basename $ANNO`.gz
    gunzip `basename $ANNO`.gz
done

mv GCA_902459845.2_HEP1_genomic.gtf hep.gtf #rename gtf file

rm GCA_*

cd ..
