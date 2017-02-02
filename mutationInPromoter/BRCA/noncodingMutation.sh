#########################################################################
# File Name: main.sh
# Author: C.J. Liu
# Mail: chunjie.sam.liu@gmail.com
# Created Time: Tue 31 Jan 2017 02:36:36 PM CST
#########################################################################
#!/bin/bash

ROOT=`dirname $0`
# echo $ROOT
# Rscript for somatic mutation
Rscript ${ROOT}/02.getSomatic.R

tmp_fifofile='/tmp/$$.fifo'
mkfifo $tmp_fifofile
exec 6<>$tmp_fifofile
rm -rf $tmp_fifofile

for (( i=0; i < 50; i++ ))
do
	echo ""
done>&6

inputDir='/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis/03.somatic/01.annotation'

avinputs=(`find $inputDir -name "*avinput" -type f`)
for avinput in ${avinputs[@]}
do
    read -u6
    {
        bash ${ROOT}/removedbSNP.sh ${avinput}
        python ${ROOT}/regulatoryFeature.py ${avinput}.dbsnp.sortByChrom > ${avinput}.dbsnp.sortByChrom.region
        echo "">&6
    }&
done
wait
exec 6>&-

Rscript ${ROOT}/03.loadAnnoSomatic.R








