mode=regular
# To make colocalizations comparable to mQTL-derived ones, we consider 500kb cis eGene region instead of 1Mb. The bed files are generated accordingly
#echo -e 'gene_or_cpg\ttissue\tfdr' > /gpfs/data/gtex-group/mpavia/methylation/mashr/eQTLs/fdr.005.txt
for TISSUE in $(cat data/tissues.txt); do 
	tissue=$TISSUE
	if [[ $TISSUE =~ .*_.* ]] ; then tissue=$(echo $TISSUE | sed 's/_//g'); fi
	echo $tissue; 
#	rm /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.all_eQTLs.bed;
#	zcat /gpfs/data/gtex-group/v8/59349/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL/$TISSUE.v8.egenes.txt.gz | awk -v tissue=$tissue '{ if( $29<0.05 ) { print $1"\t"tissue"\t"$29 }}' >> /gpfs/data/gtex-group/mpavia/methylation/mashr/eQTLs/fdr.005.txt
	zgrep -w -F -f <(cat /gpfs/data/gtex-group/mpavia/methylation/mashr/eQTLs/fdr.005.txt| grep -F $tissue | sort | uniq |  cut -f 1 ) /gpfs/data/gtex-group/v8/59349/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL/$TISSUE.v8.egenes.txt.gz | awk '{ if( $29<0.05 ) { if ($6 == "+") { if($4-500000-1 < 0) { print $3,1,$4+500000,$1,"+" } else { print $3,$4-500000-1,$4+500000,$1,"+" } } else { if($5-500000-1 < 0) { print $3,1,$5+500000,$1,"-" } else { print $3,$5-500000-1,$5+500000,$1,"-" } } } }' | tr ' ' '\t' > /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.all_eQTLs.bed
	head /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.all_eQTLs.bed
done
