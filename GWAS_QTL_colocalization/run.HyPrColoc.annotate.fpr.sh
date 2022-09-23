echo mQTL_shared.summary
for hit in $(sed 's/\t/,/g' data/mQTL.shared.colocs.GWAS.hits.txt | tail -n+2); do
        gwas=$(echo $hit | awk -F ',' '{print $1}');
        gwas_region=$(echo $hit | awk -F ',' '{print $2}');
# a) Select top-colocalizing (largest PP4) mCpG-tissue case
        top_tissue=$(echo $hit | awk -F ',' '{print $4}');
        top_cpg=$(echo $hit | awk -F ',' '{print $5}');
        if [[ `find results/mQTL_shared -name "$gwas.$gwas_region.$top_cpg.$top_tissue.*.HyPrColoc.txt"` ]]; then
                for file in $(find results/mQTL_shared -name "$gwas.$gwas_region.$top_cpg.$top_tissue.*.HyPrColoc.txt"); do eqtltrait=$(echo $file | cut -d'.' -f 5); if [[ ! -f results_random_eQTL/mQTL_shared/$gwas.$gwas_region.$top_cpg.$top_tissue.$eqtltrait.HyPrColoc.random_eQTL.txt || $(awk -F "\t" '{print $15}' results_random_eQTL/mQTL_shared/$gwas.$gwas_region.$top_cpg.$top_tissue.$eqtltrait.HyPrColoc.random_eQTL.txt | tail -n 1) != 1000 ]];then fpr='NA';else fpr=$(awk -F "\t" '{if($12 ~ /0.98$/ && $2 ~ /ENSG/ && $2 ~ /, cg/ && $4 > 0.8 && $3 > 0.5) {print}}' results_random_eQTL/mQTL_shared/$gwas.$gwas_region.$top_cpg.$top_tissue.$eqtltrait.HyPrColoc.random_eQTL.txt | wc -l); fi;awk -v fpr=$fpr -v gwas=$gwas -v region=$gwas_region -v tissue=$top_tissue -F "\t" '{if($12 ~ /0.98$/ && $2 ~ /ENSG/ && $2 ~ /, cg/ && $4 > 0.8 && $3 > 0.5) {print fpr"\tmQTL_shared\t",gwas,"\t",region,"\t",tissue,"\t",$0}}' $file;done | tee >([ $(wc -m) -gt 0 ] || echo -e "eQTL_no_coloc_mQTL_shared\t$gwas\t$gwas_region\t$top_tissue\t$top_cpg";)
        else
                echo -e "eQTL_no_tested_coloc_mQTL_shared\t$gwas\t$gwas_region\t$top_tissue\t$top_cpg";
        fi
done > results/mQTL_shared.summary.$thres.fpr.txt

echo mQTL_specific.summary
for hit in $(sed 's/\t/,/g' data/mQTL.specific.colocs.GWAS.hits.txt | tail -n+2); do
        gwas=$(echo $hit | awk -F ',' '{print $1}');
        gwas_region=$(echo $hit | awk -F ',' '{print $2}');
# a) Select top-colocalizing (largest PP4) mCpG-tissue case
        top_tissue=$(echo $hit | awk -F ',' '{print $4}');
        top_cpg=$(echo $hit | awk -F ',' '{print $5}');
	if [[ `find results/mQTL_specific -name "$gwas.$gwas_region.$top_cpg.$top_tissue.*.HyPrColoc.txt"` ]]; then
		for file in $(find results/mQTL_specific -name "$gwas.$gwas_region.$top_cpg.$top_tissue.*.HyPrColoc.txt"); do eqtltrait=$(echo $file | cut -d'.' -f 5); if [[ ! -f results_random_eQTL/mQTL_specific/$gwas.$gwas_region.$top_cpg.$top_tissue.$eqtltrait.HyPrColoc.random_eQTL.txt || $(awk -F "\t" '{print $15}' results_random_eQTL/mQTL_specific/$gwas.$gwas_region.$top_cpg.$top_tissue.$eqtltrait.HyPrColoc.random_eQTL.txt | tail -n 1) != 1000 ]];then fpr='NA';else fpr=$(awk -F "\t" '{if($12 ~ /0.98$/ && $2 ~ /ENSG/ && $2 ~ /, cg/ && $4 > 0.8 && $3 > 0.5) {print}}' results_random_eQTL/mQTL_specific/$gwas.$gwas_region.$top_cpg.$top_tissue.$eqtltrait.HyPrColoc.random_eQTL.txt | wc -l); fi;awk -v fpr=$fpr -v gwas=$gwas -v region=$gwas_region -v tissue=$top_tissue -F "\t" '{if($12 ~ /0.98$/ && $2 ~ /ENSG/ && $2 ~ /, cg/ && $4 > 0.8 && $3 > 0.5) {print fpr"\tmQTL_specific\t",gwas,"\t",region,"\t",tissue,"\t",$0}}' $file;done | tee >([ $(wc -m) -gt 0 ] || echo -e "eQTL_no_coloc_mQTL_specific\t$gwas\t$gwas_region\t$top_tissue\t$top_cpg";)
	else
		echo -e "eQTL_no_tested_coloc_mQTL_specific\t$gwas\t$gwas_region\t$top_tissue\t$top_cpg";
	fi
done > results/mQTL_specific.summary.$thres.fpr.txt
