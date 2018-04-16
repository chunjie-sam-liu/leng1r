#########################################################################
# File Name: mutect2Call.sh
# Author: C.J. Liu
# Mail: chunjie.sam.liu@gmail.com
# Created Time: Wed 01 Feb 2017 09:55:26 AM CST
#########################################################################
#!/bin/bash

GATK="/home/cliu18/liucj/tools/GenomeAnalysisTK-3.7/GenomeAnalysisTK-3.7.jar"
REF="/home/cliu18/liucj/pipelines/exome_pipeline/data/hg38/genomeBuild/hg38.fasta"
tumor="/extraspace/TCGA/WXS_RAW/BRCA/regulatoryBam/forTest/mutect2/552d2edb157ccf877109a7fa1b5e21b7_gdc_realn.extracted.bam"
normal="/extraspace/TCGA/WXS_RAW/BRCA/regulatoryBam/forTest/mutect2/2f1234605b512497e713e21f7978ff2e_gdc_realn.extracted.bam"
output="/extraspace/TCGA/WXS_RAW/BRCA/regulatoryBam/forTest/mutect2/output.vcf"

java -jar $GATK \
 -T MuTect2 \
 -R $REF \
 -I:tumor ${tumor} \
 -I:normal ${normal} \
 -o ${output}