echo mQTL_shared.summary
for hit in $(sed 's/\t/,/g' data/mQTL.shared.colocs.GWAS.hits.txt | tail -n+2); do
        gwas=$(echo $hit | awk -F ',' '{print $1}');
        gwas_region=$(echo $hit | awk -F ',' '{print $2}');
# a) Select top-colocalizing (largest PP4) mCpG-tissue case
        top_tissue=$(echo $hit | awk -F ',' '{print $4}');
        top_cpg=$(echo $hit | awk -F ',' '{print $5}');
        if [[ `find results/mQTL_shared -name "$gwas.$gwas_region.$top_cpg.$top_tissue.*.HyPrColoc.txt"` ]]; then
                for file in $(find results/mQTL_shared -name "$gwas.$gwas_region.$top_cpg.$top_tissue.*.HyPrColoc.txt"); do awk -v gwas=$gwas -v region=$gwas_region -v tissue=$top_tissue -F "\t" '{if($12 ~ /$thres$/ && $2 ~ /ENSG/ && $2 ~ /, cg/ && $4 > 0.8 && $3 > 0.5) {print "mQTL_shared\t",gwas,"\t",region,"\t",tissue,"\t",$0}}' $file;done | tee >([ $(wc -m) -gt 0 ] || echo -e "eQTL_no_coloc_mQTL_shared\t$gwas\t$gwas_region\t$top_tissue\t$top_cpg";)
        else
                echo -e "eQTL_no_tested_coloc_mQTL_shared\t$gwas\t$gwas_region\t$top_tissue\t$top_cpg";
        fi
done > results/mQTL_shared.summary.$thres.txt

echo mQTL_specific.summary
for hit in $(sed 's/\t/,/g' data/mQTL.specific.colocs.GWAS.hits.txt | tail -n+2); do
        gwas=$(echo $hit | awk -F ',' '{print $1}');
        gwas_region=$(echo $hit | awk -F ',' '{print $2}');
# a) Select top-colocalizing (largest PP4) mCpG-tissue case
        top_tissue=$(echo $hit | awk -F ',' '{print $4}');
        top_cpg=$(echo $hit | awk -F ',' '{print $5}');
	if [[ `find results/mQTL_specific -name "$gwas.$gwas_region.$top_cpg.$top_tissue.*.HyPrColoc.txt"` ]]; then
		for file in $(find results/mQTL_specific -name "$gwas.$gwas_region.$top_cpg.$top_tissue.*.HyPrColoc.txt"); do awk -v gwas=$gwas -v region=$gwas_region -v tissue=$top_tissue -F "\t" '{if($12 ~ /$thres$/ && $2 ~ /ENSG/ && $2 ~ /, cg/ && $4 > 0.8 && $3 > 0.5) {print "mQTL_specific\t",gwas,"\t",region,"\t",tissue,"\t",$0}}' $file;done | tee >([ $(wc -m) -gt 0 ] || echo -e "eQTL_no_coloc_mQTL_specific\t$gwas\t$gwas_region\t$top_tissue\t$top_cpg";)
	else
		echo -e "eQTL_no_tested_coloc_mQTL_specific\t$gwas\t$gwas_region\t$top_tissue\t$top_cpg";
	fi
done > results/mQTL_specific.summary.$thres.txt
