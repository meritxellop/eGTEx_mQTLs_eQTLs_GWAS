for tissue in $(grep -v Kidney ../coloc/data/tissues.txt); do
	TISSUE=$(echo $tissue | sed 's/_//g');
	echo $TISSUE
	echo 'cpg:gene','CpG_to_TSS_dist','rho_lower','rho','row_upper','pvalue','padj','ttest','ttest_se' | tr ',' '\t' > results/$TISSUE.cor.stats.signif.cpg.gene.pairs.qvalue.005.txt;
	num_signif=$(tail -n+2 results/signif.cpg.gene.pairs.qvalue.005.txt | wc -l)
	join --nocheck-order -a 1 <(tail -n+2 results/signif.cpg.gene.pairs.qvalue.005.txt | awk '{print $1":"$2}' | sort -k1,1) <(grep -m $num_signif -F -f <(tail -n+2 results/signif.cpg.gene.pairs.qvalue.005.txt | awk '{print $1"\t"$2}') results/$TISSUE.cor.stats.txt | awk '{print $1":"$2"\t"$0}' | sort -k1,1) | tr ' ' '\t' | cut -f 1,4-11 >> results/$TISSUE.cor.stats.signif.cpg.gene.pairs.qvalue.005.txt;
done
