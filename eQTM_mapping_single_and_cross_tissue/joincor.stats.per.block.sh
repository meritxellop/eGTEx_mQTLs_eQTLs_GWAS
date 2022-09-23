for tissue in $(grep -v Kidney ../coloc/data/tissues.txt); do
	TISSUE=$(echo $tissue | sed 's/_//g');
	echo 'cpg','gene','CpG_to_TSS_dist','rho_lower','rho','row_upper','pvalue','padj','ttest','ttest_se' | tr ',' '\t' > results/$TISSUE.cor.stats.txt;
	for num in $(cat ../data/support_files/fileids.txt); do 
		echo $tissue $num; 
		file=/scratch/mpavia/fastQTL/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles/$TISSUE.$num.cor.stats.txt
		awk 'FNR==NR{a[$1]=$2;next}{ if ($7*a[$1] < 1) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$7*a[$1]"\t"$8"\t"$9 } else { print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"1"\t"$8"\t"$9 } }' <(awk '{print $1}' $file | uniq -c | awk '{print $2"\t"$1}') $file >> results/$TISSUE.cor.stats.txt;
	done;
	awk ' NR==1 { cpg=$1;cpgpval[cpg]=$7;line=$1"\t"$2"\t"$7"\t"$8 } { if ($1 == cpg) { if ($7 < cpgpval[cpg]) { cpgpval[cpg]=$7; line=$1"\t"$2"\t"$7"\t"$8  } } else { print line; cpg=$1; cpgpval[cpg]=$7; line=$1"\t"$2"\t"$7"\t"$8  } } END { print line }' results/$TISSUE.cor.stats.txt > results/$TISSUE.cor.stats.min.pval.txt;
	echo DONE
done
