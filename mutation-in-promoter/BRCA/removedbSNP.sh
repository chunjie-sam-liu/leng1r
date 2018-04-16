#########################################################################
# File Name: removedbSNP.sh
# Author: C.J. Liu
# Mail: chunjie.sam.liu@gmail.com
# Created Time: Tue 31 Jan 2017 02:36:36 PM CST
#########################################################################
#!/bin/bash



function usage {
  if [ "${para}" -ne 1 ]; then
    echo "Error: Input 1 avinput file"
    echo "Example: bash removedbSNP.sh vcf.avinput"
    exit 1
  fi
  if [ ! -f "${avinput}" ]; then
    echo "Error: Input is not file"
    echo "Example: bash removedbSNP.sh vcf.avinput"
    exit 1
  fi
  if [ ${avinput##*.} != "avinput" ]; then
    echo "Error: Input must be ANNOVAR recoganized file"
    echo "Example: bash removedbSNP.sh vcf.avinput"
    exit 1
  fi
  
}

para=$#
avinput=$1

usage

annotate_variation='/home/cliu18/liucj/pipelines/exome_pipeline/software/annovar/annotate_variation.pl'
humandb38='/home/cliu18/liucj/pipelines/exome_pipeline/data/hg38/humandb38'
sortchrom='/home/cliu18/liucj/scripts/sortChromosomeAndPosition.py'

perl ${annotate_variation} ${avinput} ${humandb38} -filter -build hg38 -dbtype avsnp147

droped=${avinput}.hg38_avsnp147_dropped
awk '{print $3,$4,$6,$7,$2,$8,$9,$10,$11,$12}' ${droped} |sed -e "s/  / /g" -e "s/ /\t/g" > ${droped}.reorder

filtered=${avinput}.hg38_avsnp147_filtered
awk '{print $1,$2,$4,$5,".",$6,$7,$8,$9,$10}' ${filtered} |sed -e "s/  / /g" -e "s/ /\t/g" > ${filtered}.reorder

echo "NOTICE: Combined dbSNP"
cat ${droped}.reorder ${filtered}.reorder > ${avinput}.dbsnp

echo "NOTICE: Sort dbSNP"
python ${sortchrom} ${avinput}.dbsnp

rm ${droped}* ${filtered}* ${avinput}.log ${avinput}.dbsnp 

exit 1