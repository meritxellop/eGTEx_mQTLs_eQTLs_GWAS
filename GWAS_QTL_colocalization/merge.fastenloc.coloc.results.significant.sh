pvaluethres=5e-08
mode=regular


echo GWAS,GWAS_region,QTL_Tissue,QTL_eGene_or_mCpG,QTL_eGene_or_mCpG_rank,coloc_PP4,fastenloc_RCP,fastenloc_Num_SNP,fastenloc_CPIP_qtl,fastenloc_CPIP_gwas_marginal,fastenloc_CPIP_gwas_qtl_prior,fastenloc_SNP,fastenloc_PIP_qtl,fastenloc_PIP_gwas_marginal,fastenloc_PIP_gwas_qtl_prior,fastenloc_SCP,QTL_variant_id_top,QTL_rs_id_top,QTL_Signal_id,QTL_maf_top,QTL_slope_top,QTL_slope_se_top,QTL_pval_top,QTL_pval_beta_top,QTL_qval_top,GWAS_variant_id_top,GWAS_pval_top | tr ',' '\t' > results/summary/standard_GWAS.all_eQTLs.pp4_010.rcp_010.txt;
#echo GWAS,GWAS_region,QTL_Tissue,QTL_eGene_or_mCpG,QTL_eGene_or_mCpG_rank,coloc_PP4,fastenloc_RCP,fastenloc_Num_SNP,fastenloc_CPIP_qtl,fastenloc_CPIP_gwas_marginal,fastenloc_CPIP_gwas_qtl_prior,fastenloc_SNP,fastenloc_PIP_qtl,fastenloc_PIP_gwas_marginal,fastenloc_PIP_gwas_qtl_prior,fastenloc_SCP,QTL_variant_id_top,QTL_rs_id_top,QTL_Signal_id,QTL_maf_top,QTL_slope_top,QTL_slope_se_top,QTL_pval_top,QTL_pval_beta_top,QTL_qval_top,GWAS_variant_id_top,GWAS_pval_top | tr ',' '\t' > results/summary/standard_GWAS.all_mQTLs.pp4_010.rcp_010.txt;

for gwas in $(cat dap/data/GWASes.87.txt);do
	echo $gwas
	for tissue in $(cat data/tissues.txt|sed 's/_//g'); do
		echo $gwas $tissue
		file1=results/summary/coloc.standard_GWAS.all_eQTLs.pp4_010.txt
		file2=/gpfs/data/pierce-lab/GTEx/mQTLs/colocalizations/fastenloc/$tissue/$tissue.$gwas.all.eQTLs.summary.txt.gz
		join <(awk '{print $1":"$2":"$3":"$4"\t"$0}' $file1 | sort -k1,1) <(zcat $file2 | tr "," "\t" | awk '{if($9>0.1) {print $1":"$2":"$3":"$4"\t"$0}}' | sort -k1,1) | awk '{print $1,$6,$15,$11,$12,$13,$14,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31}' | tr ':' ' ' | tr ' ' '\t' >> results/summary/standard_GWAS.all_eQTLs.pp4_010.rcp_010.txt;
		file1=results/summary/coloc.standard_GWAS.all_mQTLs.pp4_010.txt
		file2=/gpfs/data/pierce-lab/GTEx/mQTLs/colocalizations/fastenloc/$tissue/$tissue.$gwas.all.summary.txt.gz
#		join <(awk '{print $1":"$2":"$3":"$4"\t"$0}' $file1 | sort -k1,1) <(zcat $file2 | tr "," "\t" | awk '{if($9>0.1) {print $1":"$2":"$3":"$4"\t"$0}}' | sort -k1,1) | awk '{print $1,$6,$15,$11,$12,$13,$14,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31}' | tr ':' ' ' | tr ' ' '\t' >> results/summary/standard_GWAS.all_mQTLs.pp4_010.rcp_010.txt;
	done;
done;
awk '{if($6>0.5 && $7>0.3) {print}}' results/summary/standard_GWAS.all_eQTLs.pp4_010.rcp_010.txt > results/summary/standard_GWAS.all_eQTLs.pp4_050.rcp_030.txt
awk '{if($6>0.3 && $7>0.3) {print}}' results/summary/standard_GWAS.all_eQTLs.pp4_010.rcp_010.txt > results/summary/standard_GWAS.all_eQTLs.pp4_030.rcp_030.txt
#awk '{if($6>0.5 && $7>0.3) {print}}' results/summary/standard_GWAS.all_mQTLs.pp4_010.rcp_010.txt > results/summary/standard_GWAS.all_mQTLs.pp4_050.rcp_030.txt
#awk '{if($6>0.3 && $7>0.3) {print}}' results/summary/standard_GWAS.all_mQTLs.pp4_010.rcp_010.txt > results/summary/standard_GWAS.all_mQTLs.pp4_030.rcp_030.txt
