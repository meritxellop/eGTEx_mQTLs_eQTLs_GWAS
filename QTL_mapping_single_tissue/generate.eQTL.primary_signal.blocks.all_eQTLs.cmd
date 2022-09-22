tissue=$1
if [[ $tissue =~ .*_.* ]] ; then TISSUE=$(echo $tissue | sed 's/_//g'); else TISSUE=$tissue; fi
block=$2
tmpdir='/scratch/mpavia/fastQTL'
wdir=/gpfs/data/gtex-group/mpavia/methylation

if [[ -f $tmpdir/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles/$TISSUE.$block.regular.all_eQTLs.txt.gz && -f $tmpdir/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles/all_eQTLs."$TISSUE"_"$block".bed.gz ]]; then
	exit
fi
	echo List genes $tissue $block

echo "
#PBS -N generate_all_eQTLs
#PBS -S /bin/bash
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=2000mb
#PBS -d $wdir
#PBS -o $wdir/logs/generate_all_eQTLs_"$tissue"_"$block"_\$PBS_JOBID.out
#PBS -e $wdir/logs/generate_all_eQTLs_"$tissue"_"$block"_\$PBS_JOBID.err

module load htslib
date;

eQTLfile=/gpfs/data/gtex-group/v8/59349/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_all_associations/$tissue.allpairs.txt.gz
bedFile=/gpfs/data/gtex-group/v8/59349/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices/$tissue.v8.normalized_expression.bed.gz

# Extract eQTL stats for 200-genes block, subsetting for significant eGenes
# All variants in cis (+/-1Mb) tested for eGenes [FDR < 0.05 & LFSR < 0.05] are considered
# For colocalization purposes, a smaller cis region, e.g. 500kb can be considered. The variants will be subsetted.

zgrep -F -f <(cat /gpfs/data/gtex-group/mpavia/methylation/mashr/eQTLs/lfsr.005.txt /gpfs/data/gtex-group/mpavia/methylation/mashr/eQTLs/fdr.005.txt | grep -F $TISSUE | cut -f 1 | sort | uniq | grep -F -f $wdir/data/support_files/lists/cpgs_by_fileid/$TISSUE.$block.genes.all_eQTLs.txt ) \$eQTLfile > $tmpdir/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles/$TISSUE.$block.regular.all_eQTLs.txt
bgzip -f $tmpdir/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles/$TISSUE.$block.regular.all_eQTLs.txt;


bedFileSplit=$tmpdir/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles/all_eQTLs."$TISSUE"_"$block".bed
zcat \$bedFile | head -n 1 > \$bedFileSplit;

if [[ ! -s $wdir/data/support_files/lists/cpgs_by_fileid/$TISSUE.$block.genes.all_eQTLs.txt ]]; then
	echo List genes $TISSUE $block
	file=/scratch/mpavia/fastQTL/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles/all_eQTLs."$TISSUE"_"$block".gz
	zcat $file | tail -n+2 | cut -f 4 > $wdir/data/support_files/lists/cpgs_by_fileid/$TISSUE.$block.genes.all_eQTLs.txt
fi
zgrep -f $wdir/data/support_files/lists/cpgs_by_fileid/$TISSUE.$block.genes.all_eQTLs.txt \$bedFile | sort -k1,1V -k2,2n  >> \$bedFileSplit;
bgzip -f \$bedFileSplit;
tabix -f -p bed \$bedFileSplit.gz;

date;
 
" >  $wdir/logs/generate_all_eQTLs.$tissue.$block.sh 
#bash $wdir/logs/generate_all_eQTLs.$tissue.$block.sh
qsub $wdir/logs/generate_all_eQTLs.$tissue.$block.sh
