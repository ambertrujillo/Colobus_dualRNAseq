#!/bin/sh

# Script to download Hepatocystis genome, hep

module load kent/385

cd genomes/

GENOME_FA=hep.fa

HEPATOCYSTIS_URL=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/902/459/845
HEPATOCYSTIS_URL=$HEPATOCYSTIS_URL/GCA_902459845.2_HEP1/GCA_902459845.2_HEP1_genomic.fna.gz 

wget $HEPATOCYSTIS_URL \
    -O ${GENOME_FA}.gz
gunzip -c ${GENOME_FA}.gz | \
    sed -e "s/^>.*chromosome \([^,]*\),.*/>Hepchr\\1/" > \
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

cd ..
