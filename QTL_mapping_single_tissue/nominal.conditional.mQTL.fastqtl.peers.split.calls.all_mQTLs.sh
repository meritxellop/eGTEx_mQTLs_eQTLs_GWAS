module load gcc/6.2.0 fastqtl/2.184_gtex bcftools
module load htslib

tissue=$1 # e.g. Lung 
id=$2 # e.g. 0000
tmpdir=/scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles
wdir=/gpfs/data/gtex-group/mpavia/methylation
vcfFile=$wdir/data/Genotypes/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz
resultsFile=$tmpdir/$id.regular.perm.txt
# Generate conditional calls for FDR < 0.05 CpGs and LFSR < 0.05 CpGs. To obtain LFSR < 0.05 CpGs, one needs to run mash first.
grep -F -f <(cat /gpfs/data/gtex-group/mpavia/methylation/mashr/mQTLs/lfsr.005.txt /gpfs/data/gtex-group/mpavia/methylation/mashr/mQTLs/fdr.005.txt  | grep -F $tissue | sort | uniq |  cut -f 1 | grep -f $wdir/data/support_files/lists/cpgs_by_fileid/$id.cpgs.txt) $wdir/results/permuted/$tissue.regular.perm.txt> $resultsFile
bedFile=$tmpdir/"$tissue"_"$id".gz
covFile=$wdir/data/Covariates/$tissue.covariates.txt
thres=$(awk '{if($8 < 0.05) {print $7}}' $wdir/results/permuted/$tissue.regular.perm.fdr.txt | sort -rgk1 | head -n 1)
samplesFile=$wdir/data/support_files/$tissue.mQTL.samples.txt
perms=1000
window=500000
# MAF=1% per tissue FILTER!!!
# Window=500Kb!!!

cp $resultsFile $resultsFile.tmp
> /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$id.bwd.all_mQTLs.all_cpgs.txt
> /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$id.fwd.all_mQTLs.all_cpgs.txt
for cpg in $(cut -f 1 $resultsFile.tmp); do
        echo $cpg
        grep $cpg $resultsFile.tmp > $resultsFile
        # Call 2ary signals. Top SNP per independent signal kept.
        timeout 5m python3.8 $wdir/scripts/multiple_eqtl_mapping-master/multiple_eqtl_mapping.mod.py --vcf $vcfFile --bed $bedFile --results $resultsFile --backward $tissue.$id.bwd.all_mQTLs.txt --forward $tissue.$id.fwd.all_mQTLs.txt --cov $covFile --overwrite --fastQTL /apps/software/gcc-6.2.0/fastqtl/2.184_gtex/bin/fastQTL -W $window --perm $perms --include-samples $samplesFile --threshold $thres --tmpdir $tmpdir
        cat /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$id.bwd.all_mQTLs.txt >> /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$id.bwd.all_mQTLs.all_cpgs.txt
        cat /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$id.fwd.all_mQTLs.txt >> /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$id.fwd.all_mQTLs.all_cpgs.txt
        output=$?
        echo PROGRAM_OUTPUT $output
        if [[ $output != 0 ]]; then echo "Gene $gene did not work. Oops!"; rm $tmpdir/Covariates.$gene.$tissue.$id.fwd.all_mQTLs.txt; else echo "Gene $gene worked. Eureka!";fi
done
mv $resultsFile.tmp $resultsFile
mv /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$id.bwd.all_mQTLs.all_cpgs.txt /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$id.bwd.all_mQTLs.ALL_cpgs.txt
mv /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$id.fwd.all_mQTLs.all_cpgs.txt /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$id.fwd.all_mQTLs.ALL_cpgs.txt
if [[ ! -s $resultsFile ]]; then exit;fi

# A) Generate Covariates for CpGs with >1 independent signal, for CpGs from block $id
if [[ $(grep -P '\t2$' /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$id.bwd.all_mQTLs.ALL_cpgs.txt | wc -l) > 0 ]]; then
	for cpg in $(grep -P '\t2$' /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$id.bwd.all_mQTLs.ALL_cpgs.txt | cut -f 1); do
		 rank=0;
                # For each lead SNP representing a eQTL independent signal  
                for snp in $(grep chr $tmpdir/Covariates.$cpg.$tissue.$id.fwd.all_mQTLs.txt | cut -f 1); do
                        rank=$(($rank+1));
			echo $cpg $rank
			cat <(grep -v chr $tmpdir/Covariates.$cpg.$tissue.$id.fwd.all_mQTLs.txt | grep -P '.'; \
                        grep chr $tmpdir/Covariates.$cpg.$tissue.$id.fwd.all_mQTLs.txt | grep -v $snp | grep -P '.') > $tmpdir/Covariates.$cpg.$rank.$tissue.$id.fwd.txt;

			# Generate p-value landscape for CpGs with >1 independent signal
			# The p-value landscape corresponds to a lead SNP representing a mQTL independent signal, adjusted for any lead SNP representing other mQTL independent signals found for that CpG
			bash $wdir/scripts/nominal.conditional.mQTL.fastqtl.peers.split.sh $tissue $id regular $cpg $rank
		done
	done
fi

rm $tmpdir/Covariates.*.$tissue.$id.*.txt
