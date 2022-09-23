wdir='/gpfs/data/gtex-group/mpavia/methylation/coloc';
tmpdir='/scratch/mpavia/fastQTL'

tissue=$1 # e.g. Lung
mode=$2 # e.g. regular
gwas=$3 # e.g. BCAC_Overall_BreastCancer_EUR
gwasp=$4 # e.g. 5e-08
gwasr=$5 # e.g. region0
p1=$6
p2=$7
p12=$8
type=$(grep $gwas $wdir/data/GWAS.metadata.txt | cut -f 21 | awk '{if($1 == 0) {print "quant"} else {print "cc"}}')
cases=$(grep $gwas $wdir/data/GWAS.metadata.txt | cut -f 22)
ss=$(grep $gwas $wdir/data/GWAS.metadata.txt | cut -f 16)
i=$(( $(echo $gwasr | sed 's/region//g')+1 ))
ldblocksFile=$wdir/dap/data/eur_ld.hg38.bed
gwasid=$gwasr.$(cat $ldblocksFile | tail -n+2 | sed "${i}q;d" | cut -f 1-3 | sed 's/\t/_/g')
start=$(cat $ldblocksFile | tail -n+2 | sed "${i}q;d" | cut -f 2)
end=$(cat $ldblocksFile | tail -n+2 | sed "${i}q;d" | cut -f 3)
TISSUE=$(echo $tissue | sed 's/BreastMammaryTissue/Breast_Mammary_Tissue/g' | sed 's/ColonTransverse/Colon_Transverse/g' | sed 's/KidneyCortex/Kidney_Cortex/g' | sed 's/MuscleSkeletal/Muscle_Skeletal/g' | sed 's/WholeBlood/Whole_Blood/g' );
eQTLfile=/gpfs/data/gtex-group/v8/59349/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_all_associations/$TISSUE.allpairs.txt.gz

#if [[ -f /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.merged.$gwasp.$gwasid.all_eQTLs.bed ]]; then
#	exit
#fi
echo "$tissue $gwas $gwasr"
if [[ ! -d $wdir/results/$tissue/$gwas/colocsummary ]]; then mkdir -p $wdir/results/$tissue/$gwas/colocsummary; fi

echo "
#PBS -N perform_coloc_fastqtl
#PBS -S /bin/bash
#PBS -l walltime=4:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=3GB
#PBS -d $wdir
#PBS -o $tmpdir/logs/perform_coloc_"$tissue"_"$gwas"_"$mode"_\$PBS_JOBID.out
#PBS -e $tmpdir/logs/perform_coloc_"$tissue"_"$gwas"_"$mode"_\$PBS_JOBID.err

date;
module load gcc/6.2.0 
module load R/3.4.1
module load bedtools
module load htslib

echo 'Select signif ( FDR<0.05 ) eQTL loci overlapping GWAS (P<$gwasp) GWAS loci, for GWAS region $gwasr'
echo 'For that, we intersect GWAS with eQTL loci coordinates, for signif regions'

i=$i
echo 'Consider region coordinates of a) eQTL locus +/- 500Kb and b) LD-defined GWAS regions'
# Generate these FDR < 0.05 eQTL cis-window files with the generate.genes.bedFiles.all_eQTLs.sh script in QTL_mapping_single_tissue
cat $ldblocksFile | tail -n+2 | sed \"\${i}q;d\" | cut -f 1-3 | bedtools intersect -a stdin -b /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.all_eQTLs.bed -wb > /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.merged.$gwasp.$gwasid.all_eQTLs.bed 
 
echo $gwasid
if [[ ! -s /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.merged.$gwasp.$gwasid.all_eQTLs.bed ]]; then
        echo 'No significant (FDR<0.05) eQTLs for $tissue overlapping gwas $gwas $gwasid'
        exit
fi

if [[ ! -f /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.$gwasp.$gwasid.all_eQTLs.txt ]]; then
	echo \"Extract gwas $gwas $gwasid and eQTL data (for eQTLs with FDR<0.05 and GWAS loci at P < $gwasp ) for $tissue $mode. Merge\"

	# Extract GWAS stats from eQTL locus +/- 500Kb. Trim by GWAS region limits. Selected region is the intersection of eQTL locus +/- 500Kb and LD-defined GWAS region
	join --nocheck-order -1 2 -2 2 <(zgrep -F -f  <(cut -f 7 /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.merged.$gwasp.$gwasid.all_eQTLs.bed | sort | uniq) $eQTLfile | awk '{ if((\$3^2)^(1/2) <= 500000 ) {print}}' | awk -F '_' -v start=$start -v end=$end '{if (\$2 > start && \$2 < end) {print \$0}}' | sort -k2,2) <(tail -n+2 /scratch/mpavia/GWASes/$gwas.txt) > /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.$gwasp.$gwasid.all_eQTLs.txt

	echo \"Extract genes of eQTL $tissue $mode\"
	cut -d' ' -f 2 /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.$gwasp.$gwasid.all_eQTLs.txt | sort | uniq > /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.$gwasp.$gwasid.genes.all_eQTLs.txt 
fi
if [[ ! -s /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.$gwasp.$gwasid.all_eQTLs.txt ]]; then
	echo \"Genes +/- 500Kb region overlap LD-defined GWAS locus but no snp appears. PROBLEM!. Skip\"
	exit
fi

echo \"Perform coloc of gwas $gwas and eQTL data for $tissue $mode\"
#$wdir/scripts/coloc.test.GWAS_thres.87_GWAS.all_eQTLs.R --tissue $tissue --mode $mode --gwas $gwas --gwasp $gwasp --gwasid $gwasid  --sample_size $ss --cases $cases --type $type --p1 $p1 --p2 $p2 --p12 $p12


# check for secondary signals in any gene of the block
if [[ ! -f /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.bwd.all_eQTLs.ALL_genes.txt ]]; then cat /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.*.bwd.all_eQTLs.ALL_genes.txt > /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.bwd.all_eQTLs.ALL_genes.txt;fi
if [[ \$(grep -P '\t2\$' /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.bwd.all_eQTLs.ALL_genes.txt | wc -l) > 0 ]]; then
	>  /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$gwas.$gwasid.all.$mode.all_eQTLs.txt;
	# for each overlapping gene with 2ary signal
	for gene in \$(grep -P '\t2\$' /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.bwd.all_eQTLs.ALL_genes.txt | cut -f 1); do
		found=\$(awk -v gene=\$gene '{if((\$7 == gene)) {print 1}}' /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.merged.$gwasp.$gwasid.all_eQTLs.bed);
		if [[ \$found > 0 ]]; then
		id=\$(grep -m 1 \$gene /gpfs/data/gtex-group/mpavia/methylation/data/support_files/lists/cpgs_by_fileid/$tissue.*.genes.all_eQTLs.txt | awk -F '.' '{print \$2}')
			echo FOUND
			#rank=0;
			for rank in \$(grep \$gene  /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.bwd.all_eQTLs.ALL_genes.txt | cut -f 18); do
				echo \$snp
			#	rank=\$((\$rank+1));
				echo \$gene \$rank
				fastqtlfile=/scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.\$id.\$gene.\$rank.$mode.all_eQTLs.txt.gz
				# Pile up all fastqtl results of eQTLs with 2ary signals in a block
				zcat \$fastqtlfile | awk '{ if((\$3^2)^(1/2) <= 500000 ) {print}}' | sed \"s/\tchr/.\$rank\tchr/g\" >> /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$gwas.$gwasid.all.$mode.all_eQTLs.txt;
			done
		fi
	done
	if [[ ! -s /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$gwas.$gwasid.all.$mode.all_eQTLs.txt ]]; then
		echo "Skip, no matches"
		rm /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$gwas.$gwasid.all.$mode.all_eQTLs.txt
		exit
	fi
	echo \"Extract gwas $gwas $gwasid and eQTL data (for eQTLs with FDR<0.05 and GWAS loci at P < $gwasp ) for $tissue $mode. Merge\"
	# Extract GWAS stats from eQTL locus +/- 500Kb. Trim by GWAS region limits. Selected region is the intersection of eQTL locus +/- 500Kb and lumped GWAS region
	join --nocheck-order -1 2 -2 2 <(cat /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$gwas.$gwasid.all.$mode.all_eQTLs.txt | awk -F '_' -v start=$start -v end=$end '{if (\$2 > start && \$2 < end) {print \$0}}' | sort -k2,2 | uniq) <(tail -n+2 /scratch/mpavia/GWASes/$gwas.txt) > /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.indep_eQTLs.$mode.$gwas.$gwasp.$gwasid.all_eQTLs.txt
	bgzip -f /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$gwas.$gwasid.all.$mode.all_eQTLs.txt;

	echo \"Extract genes of eQTL $tissue $mode\"
	cut -d' ' -f 2 /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.indep_eQTLs.$mode.$gwas.$gwasp.$gwasid.all_eQTLs.txt | sort | uniq > /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.indep_eQTLs.$mode.$gwas.$gwasp.$gwasid.genes.all_eQTLs.txt

	echo \"Perform coloc of gwas $gwas and eQTL data for $tissue $mode\"
	$wdir/scripts/coloc.test.GWAS_thres.87_GWAS.indep_eQTLs.all_eQTLs.R --tissue $tissue --mode $mode --gwas $gwas --gwasp $gwasp --gwasid $gwasid --sample_size $ss --cases $cases --type $type --p1 $p1 --p2 $p2 --p12 $p12

else
	echo 'No significant secondary signals loci for genes in $gwas $gwasid'
	exit
fi
date;
" >  $tmpdir/logs/perform_coloc_"$tissue"_"$gwas"_"$gwasp"_"$gwasid"_"$mode".standard_GWAS.all_eQTLs.qsub
qsub $tmpdir/logs/perform_coloc_"$tissue"_"$gwas"_"$gwasp"_"$gwasid"_"$mode".standard_GWAS.all_eQTLs.qsub
#ls $tmpdir/logs/perform_coloc_"$tissue"_"$gwas"_"$gwasp"_"$gwasid"_"$mode".standard_GWAS.all_eQTLs.qsub
#bash $tmpdir/logs/perform_coloc_"$tissue"_"$gwas"_"$gwasp"_"$gwasid"_"$mode".standard_GWAS.all_eQTLs.qsub
