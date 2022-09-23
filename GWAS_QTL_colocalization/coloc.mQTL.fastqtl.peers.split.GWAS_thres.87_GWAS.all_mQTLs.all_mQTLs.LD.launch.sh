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
			if [[ -f results/$tissue/$gwas/colocsummary/standard_mQTLs.$mode.$gwas.$pvaluethres.$gwasid.standard_GWAS.all_mQTLs.LD.txt ]]; then counts=$(cat results/$tissue/$gwas/colocsummary/*_mQTLs.$mode.$gwas.$pvaluethres.$gwasid.standard_GWAS.all_mQTLs.LD.txt | grep 'PP.H4.abf$' | wc -l);fi
			if [[ ! -f /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.regular.$gwas.merged.$pvaluethres.$gwasid.all_mQTLs.bed ]]; then

				echo TODO $gwas $region $tissue 
				p1=$(grep -w p1 /gpfs/data/pierce-lab/GTEx/mQTLs/colocalizations/coloc/priors/mQTLs_coloc_prior_est_by_enloc.txt | grep -F -w $tissue | grep $gwas | cut -d','  -f 4)
				p2=$(grep -w p2 /gpfs/data/pierce-lab/GTEx/mQTLs/colocalizations/coloc/priors/mQTLs_coloc_prior_est_by_enloc.txt | grep -F -w $tissue | grep $gwas | cut -d','  -f 4)
				p12=$(grep -w p12 /gpfs/data/pierce-lab/GTEx/mQTLs/colocalizations/coloc/priors/mQTLs_coloc_prior_est_by_enloc.txt | grep -F -w $tissue | grep $gwas | cut -d','  -f 4)
				bash scripts/coloc.mQTL.fastqtl.peers.split.GWAS_thres.87_GWAS.all_mQTLs.all_mQTLs.LD.cmd $tissue $mode $gwas $pvaluethres $region $p1 $p2 $p12

			elif [[ ! -s /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.regular.$gwas.merged.$pvaluethres.$gwasid.all_mQTLs.bed ]]; then

				echo "No significant (FDR<0.05) mQTLs for $tissue overlapping gwas $gwas $gwasid"

			elif [[ ! -s $wdir/results/$tissue/$gwas/colocsummary/$tissue.standard_mQTLs.$mode.$gwas.$pvaluethres.$gwasid.standard_GWAS.all_mQTLs.LD.txt ]]; then

				echo TODO $gwas $region $tissue 
				p1=$(grep -w p1 /gpfs/data/pierce-lab/GTEx/mQTLs/colocalizations/coloc/priors/mQTLs_coloc_prior_est_by_enloc.txt | grep -F -w $tissue | grep $gwas | cut -d','  -f 4)
				p2=$(grep -w p2 /gpfs/data/pierce-lab/GTEx/mQTLs/colocalizations/coloc/priors/mQTLs_coloc_prior_est_by_enloc.txt | grep -F -w $tissue | grep $gwas | cut -d','  -f 4)
				p12=$(grep -w p12 /gpfs/data/pierce-lab/GTEx/mQTLs/colocalizations/coloc/priors/mQTLs_coloc_prior_est_by_enloc.txt | grep -F -w $tissue | grep $gwas | cut -d','  -f 4)
				bash scripts/coloc.mQTL.fastqtl.peers.split.GWAS_thres.87_GWAS.all_mQTLs.all_mQTLs.LD.cmd $tissue $mode $gwas $pvaluethres $region $p1 $p2 $p12
			elif [[ $counts -gt 0 ]]; then
				echo TODO $gwas $region $tissue 
                                p1=$(grep -w p1 /gpfs/data/pierce-lab/GTEx/mQTLs/colocalizations/coloc/priors/mQTLs_coloc_prior_est_by_enloc.txt | grep -F -w $tissue | grep $gwas | cut -d','  -f 4)
                                p2=$(grep -w p2 /gpfs/data/pierce-lab/GTEx/mQTLs/colocalizations/coloc/priors/mQTLs_coloc_prior_est_by_enloc.txt | grep -F -w $tissue | grep $gwas | cut -d','  -f 4)
                                p12=$(grep -w p12 /gpfs/data/pierce-lab/GTEx/mQTLs/colocalizations/coloc/priors/mQTLs_coloc_prior_est_by_enloc.txt | grep -F -w $tissue | grep $gwas | cut -d','  -f 4)
                                bash scripts/coloc.mQTL.fastqtl.peers.split.GWAS_thres.87_GWAS.all_mQTLs.all_mQTLs.LD.cmd $tissue $mode $gwas $pvaluethres $region $p1 $p2 $p12
			fi
		done;
	done;
done;
