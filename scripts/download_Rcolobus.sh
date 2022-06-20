#!/bin/sh
  
# Script to download Colobus genome, ASM241v2

module load kent/385

mkdir -p genomes
cd genomes

GENOME_FA=Rcolobus.fa

RCOLOBUS_URL=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/591936/102
RCOLOBUS_URL=$RCOLOBUS_URL/GCF_002776525.3_ASM277652v3/GCF_002776525.3_ASM277652v3_genomic.fna.gz

wget $RCOLOBUS_URL \
    -O ${GENOME_FA}.gz
gunzip -c ${GENOME_FA}.gz | \
    sed -e "s/^>.*chromosome \([^,]*\),.*/>Colchr\\1/" > \
    $GENOME_FA

echo "Getting rid of unassembled stuff..." >&2

LAST_OK_LINE=$((`grep -n "^>[^c]" $GENOME_FA | head -n 1 | cut -d":" -f 1` - 1))
if [ $LAST_OK_LINE -gt 0 ]; then
    mv $GENOME_FA ${GENOME_FA}.backup
    head -n $LAST_OK_LINE ${GENOME_FA}.backup > ${GENOME_FA}
    rm ${GENOME_FA}.backup
fi

mkdir tmp_for_sort
faSplit byname ${GENOME_FA} tmp_for_sort/
cd tmp_for_sort/;
ls -v | xargs cat > ../${GENOME_FA}
cd ..
rm -r tmp_for_sort
