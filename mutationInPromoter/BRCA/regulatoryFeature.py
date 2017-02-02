#!/usr/bin/env python
#-*- coding:utf-8 -*-
################################################
#File Name: regulatoryFeature.py
#Author: C.J. Liu
#Mail: chunjie.sam.liu@gmail.com
#Created Time: Tue 31 Jan 2017 04:21:57 PM CST
################################################

import os,sys
RegulatoryFeaturesBed = "/home/cliu18/liucj/reference/ENSEMBL/release-86/regulations/RegulatoryFeatures.gff.sortByChrom.bed"
RegulatoryFeaturesBedIndex = "/home/cliu18/liucj/reference/ENSEMBL/release-86/regulations/RegulatoryFeatures.gff.sortByChrom.bed.idx"

def usage():
    description = '''ERROR: 
    Only one input, the input must be vcf.avinput.dbsnp.sortByChrom file
Usage:
    python addFeatureToSNV.py snp.vcf.avinput.dbsnp.sortByChrom > result.vcf
    '''
    if len(sys.argv) < 2:
        print(description)
        sys.exit(1)
    
    if not sys.argv[1].endswith("dbsnp.sortByChrom"):
        print(description)
        sys.exit(1)

def getInfoFromVCFandRegulatoryBed(s):
    arr = s.rstrip().split("\t")
    info = arr
    
    # Read bed index and read file
    index = open(RegulatoryFeaturesBedIndex, 'r')
    bed = open(RegulatoryFeaturesBed, 'r')
    
    indexArr = index.readline().rstrip().split("\t")
    
    feature = [".","."]
    
    while indexArr[0] != arr[0] or (int(arr[1]) - int(indexArr[1])) > 10000:
        indexArr = index.readline().rstrip().split("\t")
    
    seek1 = int(indexArr[2])
    seek2 = int(indexArr[3])
    
    bed.seek(seek1)
    
    while True:
        if bed.tell() > seek2: break
        
        flag = 0
        
        bedArr = bed.readline().rstrip().split("\t")
        
        if int(bedArr[1]) <= int(arr[1]) <= int(bedArr[2]):
            feature = [bedArr[3], bedArr[7]]
            flag = 1
        
        if flag == 1: break
    
    info.extend(feature)
    print(*info, sep="\t")
    
def run(vcf):
    title = "\t".join("chrom,pos,ref,alt,id,depth,refdepth,altdepth,mq,barcode,ensr,feature".split(","))
    print(title)
    with open(vcf ,'r') as foo:
        for line in foo:
            line = line.rstrip()
            getInfoFromVCFandRegulatoryBed(line)
    
def main():
    usage()
    run(sys.argv[1])

if __name__ == "__main__":
    main()