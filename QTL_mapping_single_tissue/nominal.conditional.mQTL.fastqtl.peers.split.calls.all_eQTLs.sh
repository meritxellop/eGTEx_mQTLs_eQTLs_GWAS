module load gcc/6.2.0 fastqtl/2.184_gtex bcftools
module load htslib

tissue=$1 # e.g. Lung 
if [[ $tissue =~ .*_.* ]] ; then TISSUE=$(echo $tissue | sed 's/_//g'); else TISSUE=$tissue; fi
id=$2 # e.g. 0000
tmpdir=/scratch/mpavia/fastQTL/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles
wdir=/gpfs/data/gtex-group/mpavia/methylation
vcfFile=$wdir/data/Genotypes/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz
resultsFile=$tmpdir/$TISSUE.$id.regular.perm.all_eQTLs.txt
if [[ ! -f $wdir/results/permuted/$TISSUE.regular.perm.all_eQTLs.txt ]]; then zcat /gpfs/data/gtex-group/v8/59349/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL/$tissue.v8.egenes.txt.gz | awk '{ print $1,$7,"NA","NA","NA","NA",$12,$13,"NA","NA","NA","NA",$24,$25,$26,$27,$28,$29}' | tail -n+2 | tr " " "\t" | sort -gk17 > $wdir/results/permuted/$TISSUE.regular.perm.all_eQTLs.txt; fi
# Generate conditional calls for FDR < 0.05 genes and LFSR < 0.05 genes.
grep -F -f <( cat /gpfs/data/gtex-group/mpavia/methylation/mashr/eQTLs/lfsr.005.txt /gpfs/data/gtex-group/mpavia/methylation/mashr/eQTLs/fdr.005.txt| grep -F $TISSUE | sort | uniq |  cut -f 1 | grep -f $wdir/data/support_files/lists/cpgs_by_fileid/$TISSUE.$id.genes.all_eQTLs.txt) $wdir/results/permuted/$TISSUE.regular.perm.all_eQTLs.txt | cut -f 1-17 > $resultsFile
bedFile=$tmpdir/all_eQTLs."$TISSUE"_"$id".bed.gz ###
covFile=/gpfs/data/gtex-group/v8/59349/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/$tissue.v8.covariates.txt
thres=$(awk '{if($18 < 0.05) {print $17}}' $wdir/results/permuted/$TISSUE.regular.perm.all_eQTLs.txt | sort -rgk1 | head -n 1)
echo $thres
samplesFile=$wdir/data/support_files/$TISSUE.eQTL.samples.txt
if [[ ! -f $samplesFile ]]; then 
	head -n 1 /gpfs/data/gtex-group/v8/59349/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/$tissue.v8.covariates.txt | tr "\t" "\n" | tail -n+2 > $samplesFile
fi
perms=1000
window=1000000
# NO MAF=1% per tissue FILTER FOR eQTLs!!!
# Window=1Mb instead of 500kb!!!

cp $resultsFile $resultsFile.tmp
> /scratch/mpavia/fastQTL/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles/$TISSUE.$id.bwd.all_eQTLs.all_genes.txt
> /scratch/mpavia/fastQTL/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles/$TISSUE.$id.fwd.all_eQTLs.all_genes.txt
for gene in $(cut -f 1 $resultsFile.tmp); do
	echo $gene
	grep $gene $resultsFile.tmp > $resultsFile
	# Call 2ary signals. Top SNP per independent signal kept.
	timeout 5m python3.8 $wdir/scripts/multiple_eqtl_mapping-master/multiple_eqtl_mapping.mod.matched_eQTLs.py --vcf $vcfFile --bed $bedFile --results $resultsFile --backward $TISSUE.$id.bwd.all_eQTLs.txt --forward $TISSUE.$id.fwd.all_eQTLs.txt --cov $covFile --overwrite --fastQTL /apps/software/gcc-6.2.0/fastqtl/2.184_gtex/bin/fastQTL -W $window --perm $perms --include-samples $samplesFile --threshold $thres --tmpdir $tmpdir
	cat /scratch/mpavia/fastQTL/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles/$TISSUE.$id.bwd.all_eQTLs.txt >> /scratch/mpavia/fastQTL/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles/$TISSUE.$id.bwd.all_eQTLs.all_genes.txt	
	cat /scratch/mpavia/fastQTL/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles/$TISSUE.$id.fwd.all_eQTLs.txt >> /scratch/mpavia/fastQTL/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles/$TISSUE.$id.fwd.all_eQTLs.all_genes.txt	
	output=$?
	echo PROGRAM_OUTPUT $output
	if [[ $output != 0 ]]; then echo "Gene $gene did not work. Oops!"; rm $tmpdir/Covariates.$gene.$TISSUE.$id.fwd.all_eQTLs.txt; else echo "Gene $gene worked. Eureka!";fi
done
mv $resultsFile.tmp $resultsFile
mv /scratch/mpavia/fastQTL/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles/$TISSUE.$id.bwd.all_eQTLs.all_genes.txt /scratch/mpavia/fastQTL/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles/$TISSUE.$id.bwd.all_eQTLs.ALL_genes.txt
mv /scratch/mpavia/fastQTL/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles/$TISSUE.$id.fwd.all_eQTLs.all_genes.txt /scratch/mpavia/fastQTL/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles/$TISSUE.$id.fwd.all_eQTLs.ALL_genes.txt


# A) Generate Covariates for genes with >1 independent signal, for genes from block $id
#if ls $tmpdir/Covariates.*.$TISSUE.$id.fwd.all_eQTLs.txt 1> /dev/null 2>&1; then
if [[ $(grep -P '\t2$' /scratch/mpavia/fastQTL/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles/$TISSUE.$id.bwd.all_eQTLs.ALL_genes.txt | wc -l) > 0 ]]; then
	# For each gene of A)
	for gene in $(ls $tmpdir/Covariates.*.$TISSUE.$id.fwd.all_eQTLs.txt | sed "s/.*Covariates\.//g" | sed "s/\.$TISSUE.*//g" | grep -P "ENSG\d+\.\d+$" | sort | uniq); do
		rank=0;
		# For each lead SNP representing a eQTL independent signal  
		for snp in $(grep chr $tmpdir/Covariates.$gene.$TISSUE.$id.fwd.all_eQTLs.txt | cut -f 1); do
			rank=$(($rank+1));
			echo $gene $rank
			cat <(grep -v chr $tmpdir/Covariates.$gene.$TISSUE.$id.fwd.all_eQTLs.txt | grep -P '.'; \
                        grep chr $tmpdir/Covariates.$gene.$TISSUE.$id.fwd.all_eQTLs.txt | grep -v $snp | grep -P '.') > $tmpdir/Covariates.$gene.$rank.$TISSUE.$id.fwd.all_eQTLs.txt;

			# Generate p-value landscape for gene with >1 independent signal
			# The p-value landscape corresponds to a lead SNP representing a eQTL independent signal, adjusted for any lead SNP representing other eQTL independent signals found for that gene
			bash $wdir/scripts/nominal.conditional.mQTL.fastqtl.peers.split.all_eQTLs.sh $TISSUE $id regular $gene $rank
		done
	done
fi
rm $tmpdir/Covariates.*.$TISSUE.$id.*.all_eQTLs.txt
