#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Download Plasmodium vivax annotations
# ----------------------------------------------------------------------------------------

# --- Download annotations for Plasmodium vivax genome, Pvivax


cd genomes #Think I should be in /home/aet359/baboon_malaria not papAnu3?

PVIVAX_URL=ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/Plasmodium_vivax/latest_assembly_versions

ANNO1=GCF_000002415.2_ASM241v2/GCF_000002415.2_ASM241v2_genomic.gtf
ANNO2=GCF_000002415.2_ASM241v2/GCF_000002415.2_ASM241v2_genomic.gff

for ANNO in $ANNO1 $ANNO2; do
    wget $PVIVAX_URL/$ANNO.gz \
        -O `basename $ANNO`.gz
    gunzip `basename $ANNO`.gz
done

sed '
    s/^\NC_009906.1/Plaschr1/ 
    s/^\NC_009907.1/Plaschr2/
    s/^\NC_009908.2/Plaschr3/
    s/^\NC_009909.1/Plaschr4/
    s/^\NC_009910.1/Plaschr5/
    s/^\NC_009911.1/Plaschr6/
    s/^\NC_009912.1/Plaschr7/
    s/^\NC_009913.1/Plaschr8/
    s/^\NC_009914.1/Plaschr9/
    s/^\NC_009915.1/Plaschr10/
    s/^\NC_009916.1/Plaschr11/
    s/^\NC_009917.1/Plaschr12/
    s/^\NC_009918.1/Plaschr13/
    s/^\NC_009919.1/Plaschr14/g' GCF_000002415.2_ASM241v2_genomic.gtf > GCF_000002415.2_ASM241v2_genomic.fix.gtf

mv GCF_000002415.2_ASM241v2_genomic.fix.gtf vivax.gtf #rename gtf file

cd ..
