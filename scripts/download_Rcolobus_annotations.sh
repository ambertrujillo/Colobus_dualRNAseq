#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Download Colobus annotations
# ----------------------------------------------------------------------------------------

# --- Download annotations for Piliocolobus tephrosceles genome, Rcolobus

RCOLOBUS_URL=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/591936/102

ANNO1=GCF_002776525.3_ASM277652v3/GCF_002776525.3_ASM277652v3_genomic.gtf
ANNO2=GCF_002776525.3_ASM277652v3/GCF_002776525.3_ASM277652v3_genomic.gff

for ANNO in $ANNO1 $ANNO2; do
    wget $RCOLOBUS_URL/$ANNO.gz \
        -O `basename $ANNO`.gz
    gunzip `basename $ANNO`.gz
done

# Replacing names to match those on the reference assembly

sed '
    s/^\NC_045434.1/Colchr1/ 
    s/^\NC_045435.1/Colchr2/
    s/^\NC_045436.1/Colchr3/
    s/^\NC_045437.1/Colchr4/
    s/^\NC_045438.1/Colchr5/
    s/^\NC_045439.1/Colchr6/
    s/^\NC_045440.1/Colchr7/
    s/^\NC_045441.1/Colchr8/
    s/^\NC_045442.1/Colchr9/
    s/^\NC_045443.1/Colchr10/
    s/^\NC_045444.1/Colchr11/
    s/^\NC_045445.1/Colchr12/
    s/^\NC_045446.1/Colchr13/
    s/^\NC_045447.1/Colchr14/
    s/^\NC_045448.1/Colchr15/
    s/^\NC_045449.1/Colchr16/
    s/^\NC_045450.1/Colchr17/
    s/^\NC_045451.1/Colchr18/
    s/^\NC_045452.1/Colchr19/
    s/^\NC_045453.1/Colchr20/
    s/^\NC_045454.1/Colchr21/g' GCF_002776525.3_ASM277652v3_genomic.gtf > GCF_002776525.3_ASM277652v3_genomic.fix.gtf

#rename gtf file
mv GCF_002776525.3_ASM277652v3_genomic.fix.gtf Rcolobus.fix.gtf 
