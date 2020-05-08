#!/bin/sh

# Script to download Plasmodium vivax genome, Pvivax

module load kent/328

mkdir -p genomes
cd genomes

GENOME_FA=vivax.fa

PVIVAX_URL=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/415
PVIVAX_URL=$PVIVAX_URL/GCF_000002415.2_ASM241v2/GCF_000002415.2_ASM241v2_genomic.fna.gz

wget $PVIVAX_URL \
    -O ${GENOME_FA}.gz
gunzip -c ${GENOME_FA}.gz | \
    sed -e "s/^>.*chromosome \([^,]*\),.*/>Plaschr\\1/" > \
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
