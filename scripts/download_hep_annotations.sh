#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Download Plasmodium gonderi annotations
# ----------------------------------------------------------------------------------------

# --- Download annotations for Plasmodium gonderi genome, gonderi


cd colobus_hep/genomes #Think I should be in /home/aet359/baboon_malaria not papAnu3?

HEPATOCYSTIS_URL=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/902/459/845

ANNO1=GCA_902459845.2_HEP1/GCA_902459845.2_HEP1_genomic.gtf
ANNO2=GCA_902459845.2_HEP1/GCA_902459845.2_HEP1_genomic.gff

for ANNO in $ANNO1 $ANNO2; do
    wget $HEPATOCYSTIS_URL/$ANNO.gz \
        -O `basename $ANNO`.gz
    gunzip `basename $ANNO`.gz
done

mv GCA_902459845.2_HEP1_genomic.gtf hep.gtf #rename gtf file

cd ..
