echo -e "cpg\tgene\tCpG_to_TSS_dist\tsignif_tissues" > results/signif.cpg.gene.pairs.qvalue.005.txt
cat <(for tissue in $(grep -v Kidney ../coloc/data/tissues.txt); do
	TISSUE=$(echo $tissue | sed 's/_//g');
	#thresh=$(awk '{if($5<0.05) {print $3}}' results/$TISSUE.cor.stats.min.pval.txt | sort -rgk1 | head -n 1)
	#thresh=$(sort -gk5 results/$TISSUE.cor.stats.min.pval.txt | awk '{if($5<0.05) {print $3}}'| tail -n 1);
	#awk -v thresh=$thresh '{if($7<thresh) {print $1"\t"$2"\t"$3}}' results/$TISSUE.cor.stats.txt
	#awk -v thresh=$thresh '{if($7<thresh) {print $1"\t"$2"\t"$3}}' results/$TISSUE.cor.stats.txt
	# From qvalue<0.05 CpGs, select CpG-gene correlations at Bonferroni-corrected p-value < 0.05
	grep -F -f <(awk '{if($5<0.05) {print $1}}' results/$TISSUE.cor.stats.min.pval.txt) <(awk '{if($8<0.05) {print $1"\t"$2"\t"$3}}' results/$TISSUE.cor.stats.txt)
done) | sort | uniq -c | awk '{print $2"\t"$3"\t"$4"\t"$1}' >> results/signif.cpg.gene.pairs.qvalue.005.txt

head -n 1 results/Lung.cor.stats.txt | sed "s/^/Tissue\t/g" > results/all.signif.cpg.gene.pairs.qvalue.005.txt
for tissue in $(grep -v Kidney ../coloc/data/tissues.txt); do
        TISSUE=$(echo $tissue | sed 's/_//g');
        rm results/$TISSUE.signif.cpg.gene.pairs.qvalue.005.txt
        echo $TISSUE
        # From qvalue<0.05 CpGs, select CpG-gene correlations at Bonferroni-corrected p-value < 0.05
        grep -F -f <(awk '{if($5<0.05) {print $1}}' results/$TISSUE.cor.stats.min.pval.txt) <(awk -v tissue=$TISSUE '{if($8<0.05) {print tissue"\t"$0}}' results/$TISSUE.cor.stats.txt) >> results/all.signif.cpg.gene.pairs.qvalue.005.txt
done
