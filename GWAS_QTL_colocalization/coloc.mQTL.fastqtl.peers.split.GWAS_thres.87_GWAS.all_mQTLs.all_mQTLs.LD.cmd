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

#if [[ -f /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.merged.$gwasp.$gwasid.all_mQTLs.bed ]]; then
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

echo 'Select signif ( FDR<0.05 ) mQTL loci overlapping GWAS (P<$gwasp) GWAS loci, for GWAS region $gwasr'
echo 'For that, we intersect GWAS with mQTL loci coordinates, for signif regions'

i=$i
echo 'Consider region coordinates of a) mQTL locus +/- 500Kb and b) LD-defined GWAS regions'
# Generate these FDR < 0.05 mQTL cis-window files with the generate.cpgs.bedFiles.sh script in QTL_mapping_single_tissue
if [[ ! -f /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.all_mQTLs.bed ]]; then ln -s /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.bed /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.all_mQTLs.bed;fi
cat $ldblocksFile | tail -n+2 | sed \"\${i}q;d\" | cut -f 1-3 | bedtools intersect -a stdin -b /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.all_mQTLs.bed -wb > /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.merged.$gwasp.$gwasid.all_mQTLs.bed 
 
echo $gwasid
if [[ ! -s /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.merged.$gwasp.$gwasid.all_mQTLs.bed ]]; then
        echo 'No significant (FDR<0.05) mQTLs for $tissue overlapping gwas $gwas $gwasid'
        exit
fi
if [[ ! -f /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.$gwasp.$gwasid.all_mQTLs.txt ]]; then
	echo \"Extract gwas $gwas $gwasid and mQTL data (for mQTLs with FDR<0.05 and GWAS loci at P < $gwasp ) for $tissue $mode. Merge\"
	mQTLfile=/scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$gwas.$gwasp.$gwasid.$mode.mQTLs_pile.txt
	> \$mQTLfile
	for id in \$(cat /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.merged.$gwasp.$gwasid.all_mQTLs.bed | cut -f 8 | sort | uniq); do
		zcat /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.\$id.$mode.txt.gz >> \$mQTLfile
	done
	bgzip -f \$mQTLfile

	# Extract GWAS stats from mQTL locus +/- 500Kb. Trim by GWAS region limits. Selected region is the intersection of mQTL locus +/- 500Kb and LD-defined GWAS region
	join --nocheck-order -1 2 -2 2 <(zgrep -F -f  <(cut -f 7 /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.merged.$gwasp.$gwasid.all_mQTLs.bed | sort | uniq) \$mQTLfile.gz | awk '{ if((\$3^2)^(1/2) <= 500000 ) {print}}' | awk -F '_' -v start=$start -v end=$end '{if (\$2 > start && \$2 < end) {print \$0}}' | sort -k2,2) <(tail -n+2 /scratch/mpavia/GWASes/$gwas.txt) > /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.$gwasp.$gwasid.all_mQTLs.txt
	rm \$mQTLfile.gz

	echo \"Extract cpgs of mQTL $tissue $mode\"
	cut -d' ' -f 2 /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.$gwasp.$gwasid.all_mQTLs.txt | sort | uniq > /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.$gwasp.$gwasid.cpgs.all_mQTLs.txt 
fi
if [[ ! -s /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.$gwasp.$gwasid.all_mQTLs.txt ]]; then
	echo \"CpGs +/- 500Kb region overlap LD-defined GWAS locus but no snp appears. PROBLEM!. Skip\"
	exit
fi

echo \"Perform coloc of gwas $gwas and mQTL data for $tissue $mode\"
$wdir/scripts/coloc.test.GWAS_thres.87_GWAS.all_mQTLs.R --tissue $tissue --mode $mode --gwas $gwas --gwasp $gwasp --gwasid $gwasid  --sample_size $ss --cases $cases --type $type --p1 $p1 --p2 $p2 --p12 $p12
echo $p1 $p2 $p12


# check for secondary signals in any cpg of the block
if [[ ! -f /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.bwd.all_mQTLs.ALL_cpgs.txt ]]; then cat /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.*.bwd.all_mQTLs.ALL_cpgs.txt > /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.bwd.all_mQTLs.ALL_cpgs.txt;fi
if [[ \$(grep -P '\t2\$' /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.bwd.all_mQTLs.ALL_cpgs.txt | wc -l) > 0 ]]; then
	>  /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$gwas.$gwasid.all.$mode.all_mQTLs.txt;
	# for each overlapping cpg with 2ary signal
	for cpg in \$(grep -P '\t2\$' /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.bwd.all_mQTLs.ALL_cpgs.txt | cut -f 1); do
		found=\$(awk -v cpg=\$cpg '{if((\$7 == cpg)) {print 1}}' /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.$gwas.merged.$gwasp.$gwasid.all_mQTLs.bed);
		if [[ \$found > 0 ]]; then
		id=\$(grep -m 1 \$cpg /gpfs/data/gtex-group/mpavia/methylation/coloc/data/MethylationEPIC_v-1-0_B4.filtered.hg38.CpGblocks.bed | awk '{print \$5}')
			echo FOUND
			for rank in \$(grep \$cpg  /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.bwd.all_mQTLs.ALL_cpgs.txt | cut -f 18); do
				echo \$snp
				echo \$cpg \$rank
				fastqtlfile=/scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.\$id.\$cpg.\$rank.$mode.txt.gz
				# Pile up all fastqtl results of mQTLs with 2ary signals in a block
				zcat \$fastqtlfile | awk '{ if((\$3^2)^(1/2) <= 500000 ) {print}}' | sed \"s/\tchr/.\$rank\tchr/g\" >> /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$gwas.$gwasid.all.$mode.all_mQTLs.txt;
			done
		fi
	done
	if [[ ! -s /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$gwas.$gwasid.all.$mode.all_mQTLs.txt ]]; then
		echo "Skip, no matches"
		rm /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$gwas.$gwasid.all.$mode.all_mQTLs.txt
		exit
	fi
	echo \"Extract gwas $gwas $gwasid and mQTL data (for mQTLs with FDR<0.05 and GWAS loci at P < $gwasp ) for $tissue $mode. Merge\"
	# Extract GWAS stats from mQTL locus +/- 500Kb. Trim by GWAS region limits. Selected region is the intersection of mQTL locus +/- 500Kb and lumped GWAS region
	join --nocheck-order -1 2 -2 2 <(cat /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$gwas.$gwasid.all.$mode.all_mQTLs.txt | awk -F '_' -v start=$start -v end=$end '{if (\$2 > start && \$2 < end) {print \$0}}' | sort -k2,2 | uniq) <(tail -n+2 /scratch/mpavia/GWASes/$gwas.txt) > /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.indep_mQTLs.$mode.$gwas.$gwasp.$gwasid.all_mQTLs.txt
	bgzip -f /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$gwas.$gwasid.all.$mode.all_mQTLs.txt;

	echo \"Extract cpgs of mQTL $tissue $mode\"
	cut -d' ' -f 2 /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.indep_mQTLs.$mode.$gwas.$gwasp.$gwasid.all_mQTLs.txt | sort | uniq > /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.indep_mQTLs.$mode.$gwas.$gwasp.$gwasid.cpgs.all_mQTLs.txt

	echo \"Perform coloc of gwas $gwas and mQTL data for $tissue $mode\"
	#$wdir/scripts/coloc.test.GWAS_thres.87_GWAS.indep_mQTLs.all_mQTLs.R --tissue $tissue --mode $mode --gwas $gwas --gwasp $gwasp --gwasid $gwasid --sample_size $ss --cases $cases --type $type --p1 $p1 --p2 $p2 --p12 $p12

else
	echo 'No significant secondary signals loci for cpgs in $gwas $gwasid'
	exit
fi
date;
" >  $tmpdir/logs/perform_coloc_"$tissue"_"$gwas"_"$gwasp"_"$gwasid"_"$mode".standard_GWAS.all_mQTLs.qsub
#qsub $tmpdir/logs/perform_coloc_"$tissue"_"$gwas"_"$gwasp"_"$gwasid"_"$mode".standard_GWAS.all_mQTLs.qsub
#ls $tmpdir/logs/perform_coloc_"$tissue"_"$gwas"_"$gwasp"_"$gwasid"_"$mode".standard_GWAS.all_mQTLs.qsub
bash $tmpdir/logs/perform_coloc_"$tissue"_"$gwas"_"$gwasp"_"$gwasid"_"$mode".standard_GWAS.all_mQTLs.qsub
