module load htslib;
> coloc/data/fileids.all_eQTLs.txt;
for tissue in $(cat data/support_files/tissues.txt); do
	TISSUE=$tissue
	echo $tissue;
	if [[ $tissue =~ .*_.* ]] ; then tissue=$(echo $tissue | sed 's/_//g'); fi
	echo $tissue;
	#mkdir -p /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles
	rm /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/all_eQTLs."$tissue"_*
	#rm /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/all_eQTLs."$tissue"_[0-9][0-9][0-9][0-9].gz*
	header=$(zcat /gpfs/data/gtex-group/v8/59349/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices/$TISSUE.v8.normalized_expression.bed.gz | head -n 1) 
	zcat /gpfs/data/gtex-group/v8/59349/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices/$TISSUE.v8.normalized_expression.bed.gz | tail -n+2 | split - /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/all_eQTLs."$tissue"_ -l 200 -d -a 4
	for fileid in $(ls /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/all_eQTLs."$tissue"_* | cut -d'_' -f 3); do
		echo -e "$tissue\t$fileid" >> coloc/data/fileids.all_eQTLs.txt;
		echo $fileid
		file=/scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/all_eQTLs."$tissue"_"$fileid"
		cmd=$(echo 'sed -i "1 i\'$header'"' $file);
		eval "$cmd"
		sed -i 's/ /\t/g' $file
                tail -n+2 $file | cut -f 4 > /gpfs/data/gtex-group/mpavia/methylation/data/support_files/lists/cpgs_by_fileid/$tissue.$fileid.genes.all_eQTLs.txt
		bgzip -f $file;
		tabix -f -p bed $file.gz;
	done
done
