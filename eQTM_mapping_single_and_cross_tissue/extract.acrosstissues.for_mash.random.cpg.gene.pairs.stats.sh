num_random=100000
echo 'Generate results/random.cpg.gene.pairs.txt'
date
awk {'print $1":"$2'} results/*.cor.stats.txt | grep -v cpg | sort | uniq -c | grep '8 cg'| awk '{print $2"\t"$3}' | shuf | head -n $num_random |  sort -k1,1 | tr ':' '\t' > results/random.cpg.gene.pairs.txt
echo 'Done'
date

for tissue in $(grep -v Kidney ../coloc/data/tissues.txt ); do
	TISSUE=$(echo $tissue | sed 's/_//g');
	echo $TISSUE
	echo 'cpg:gene','CpG_to_TSS_dist','rho_lower','rho','row_upper','pvalue','padj','ttest','ttest_se' | tr ',' '\t' > results/$TISSUE.cor.stats.random.cpg.gene.pairs.txt;
	join --nocheck-order -a 1 <(awk '{print $1":"$2}' results/random.cpg.gene.pairs.txt) <(grep -m $num_random -F -f results/random.cpg.gene.pairs.txt results/$TISSUE.cor.stats.txt | awk '{print $1":"$2"\t"$0}' | sort -k1,1) | tr ' ' '\t' | cut -f 1,4-11 >> results/$TISSUE.cor.stats.random.cpg.gene.pairs.txt;
done
