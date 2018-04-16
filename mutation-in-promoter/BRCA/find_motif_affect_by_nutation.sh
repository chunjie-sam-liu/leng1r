#########################################################################
# File Name: find_motif_affect_by_nutation.sh
# Author: C.J. Liu
# Mail: chunjie.sam.liu@gmail.com
# Created Time: Mon 13 Feb 2017 12:45:43 PM CST
#########################################################################
#!/bin/bash
ROOT=`dirname $0`
gene_expression="/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis/03.somatic/06.geneExpression/candidate_mutation_with_normal.tsv"

bind_region_dir='/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis/03.somatic/09.binding_region'


echo "Notice: Get mutation list"

cut -f 1 ${gene_expression} |sort|uniq|grep -v 'mutation' > ${bind_region_dir}/mutation.list

find_motif='/home/cliu18/liucj/reference/TFs/findMotif.py'

cmd="python ${find_motif} ${bind_region_dir}/mutation.list | uniq > ${bind_region_dir}/mutation_affect_binding.tsv"
echo $cmd
eval $cmd
