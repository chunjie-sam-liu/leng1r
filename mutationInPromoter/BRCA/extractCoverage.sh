#########################################################################
# File Name: extractCoverage.sh
# Author: C.J. Liu
# Mail: samliu@hust.edu.cn
# Created Time: Wed 07 Dec 2016 12:46:30 PM CST
#########################################################################
#!/bin/bash

somatic='/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis/03.somatic/03.somaticForAnalysis.saveFiles/realSomaticMutation.recur5ForAnalysis.refinePositions.tsv'

mpileups='/extraspace/TCGA/WXS_RAW/BRCA/regulatoryBam/tumor'

output="/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis/03.somatic/03.somaticForAnalysis.recheckPositions/tumor"

pipe='/tmp/$$.fifo'
mkfifo $pipe
exec 6<>$pipe
rm -rf $pipe

for (( i=0; i<50; i++ ))
do
	echo ""
done >&6

cat $somatic | while read pos
do 
	read -u6
	{
		arr=($pos)
		key="${arr[0]}\t${arr[1]}\t"
		name="${arr[0]}_${arr[1]}_${arr[2]}_${arr[3]}"
		cmd="grep --perl-regex \"$key\"  ${mpileups}/*mpileup > ${output}/${name}"
        # echo $cmd
		eval $cmd
		echo "" >&6
		
	} &
done
wait
exec 6>&-

