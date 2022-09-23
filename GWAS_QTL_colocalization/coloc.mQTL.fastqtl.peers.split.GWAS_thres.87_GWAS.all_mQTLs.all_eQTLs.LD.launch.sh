pvaluethres=5e-08
mode=regular
wdir='/gpfs/data/gtex-group/mpavia/methylation/coloc';
tmpdir='/scratch/mpavia/fastQTL'
ldblocksFile=$wdir/dap/data/eur_ld.hg38.bed.gz


for gwas in $(cat dap/data/GWASes.87.txt);do
	echo $gwas
	for region in $(grep $gwas /gpfs/data/gtex-group/mpavia/methylation/coloc/results/GWASsignifLoci/summary.independent.loci.$pvaluethres.LD.txt | cut -d' ' -f2); do
		i=$(( $(echo $region | sed 's/region//g')+1 ))
		gwasid=$region.$(zcat $ldblocksFile | tail -n+2 | sed "${i}q;d" | cut -f 1-3 | sed 's/\t/_/g')
		for tissue in $(cat data/tissues.txt|sed 's/_//g'); do
			find $wdir/results/$tissue/$gwas/colocsummary/$tissue.indep_eQTLs.$mode.$gwas.$pvaluethres.$gwasid.standard_GWAS.all_eQTLs.LD.txt ! -newermt "feb 16, 2021 09:00" -ls;
			find $wdir/results/$tissue/$gwas/colocsummary/$tissue.indep_eQTLs.$mode.$gwas.$pvaluethres.$gwasid.standard_GWAS.all_eQTLs.LD.txt ! -newermt "feb 16, 2021 09:00" | xargs rm;

			if [[ ! -f /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.regular.$gwas.merged.$pvaluethres.$gwasid.all_eQTLs.bed ]]; then

				echo TODO $gwas $region $tissue 
				p1=$(grep -w p1 /gpfs/data/pierce-lab/GTEx/mQTLs/colocalizations/coloc/priors/eQTLs_coloc_prior_est_by_enloc.txt | grep -F -w $tissue | grep $gwas | cut -d','  -f 4)
				p2=$(grep -w p2 /gpfs/data/pierce-lab/GTEx/mQTLs/colocalizations/coloc/priors/eQTLs_coloc_prior_est_by_enloc.txt | grep -F -w $tissue | grep $gwas | cut -d','  -f 4)
				p12=$(grep -w p12 /gpfs/data/pierce-lab/GTEx/mQTLs/colocalizations/coloc/priors/eQTLs_coloc_prior_est_by_enloc.txt | grep -F -w $tissue | grep $gwas | cut -d','  -f 4)
				bash scripts/coloc.mQTL.fastqtl.peers.split.GWAS_thres.87_GWAS.all_mQTLs.all_eQTLs.LD.cmd $tissue $mode $gwas $pvaluethres $region $p1 $p2 $p12

			elif [[ ! -s /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.regular.$gwas.merged.$pvaluethres.$gwasid.all_eQTLs.bed ]]; then

				echo "No significant (FDR<0.05 or LFSR<0.05) eQTLs for $tissue overlapping gwas $gwas $gwasid"

			elif [[ ! -s $wdir/results/$tissue/$gwas/colocsummary/$tissue.standard_eQTLs.$mode.$gwas.$pvaluethres.$gwasid.standard_GWAS.all_eQTLs.LD.txt ]]; then
			if [[ -f /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.$pvaluethres.$gwasid.all_eQTLs.txt && ! -s /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.$pvaluethres.$gwasid.all_eQTLs.txt ]]; then continue;fi

				echo TODO $gwas $region $tissue 
				p1=$(grep -w p1 /gpfs/data/pierce-lab/GTEx/mQTLs/colocalizations/coloc/priors/eQTLs_coloc_prior_est_by_enloc.txt | grep -F -w $tissue | grep $gwas | cut -d','  -f 4)
				p2=$(grep -w p2 /gpfs/data/pierce-lab/GTEx/mQTLs/colocalizations/coloc/priors/eQTLs_coloc_prior_est_by_enloc.txt | grep -F -w $tissue | grep $gwas | cut -d','  -f 4)
				p12=$(grep -w p12 /gpfs/data/pierce-lab/GTEx/mQTLs/colocalizations/coloc/priors/eQTLs_coloc_prior_est_by_enloc.txt | grep -F -w $tissue | grep $gwas | cut -d','  -f 4)
				bash scripts/coloc.mQTL.fastqtl.peers.split.GWAS_thres.87_GWAS.all_mQTLs.all_eQTLs.LD.cmd $tissue $mode $gwas $pvaluethres $region $p1 $p2 $p12

			fi
		done;
	done;
done;
