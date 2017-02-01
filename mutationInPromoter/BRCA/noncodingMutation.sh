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

# bash script for annovar filter dbSNP.
avinput=$1
bash ${ROOT}/removedbSNP.sh ${avinput}

# python script for regulatory build filter.
echo "NOTICE: filter regulatory region"
python ${ROOT}/regulatoryFeature.py ${avinput}.dbsnp.sortByChrom > ${avinput}.dbsnp.sortByChrom.region