wdirc='/gpfs/data/gtex-group/mpavia/methylation/coloc';
wdir='/gpfs/data/gtex-group/mpavia/methylation/HyPrColoc';
tmpdir='/scratch/mpavia/fastQTL'

mode='regular'
tissue=$1 # e.g. MuscleSkeletal
gwas=$2 # e.g. Astle_et_al_2016_Eosinophil_counts
gwasr=$3 # e.g. region1004
cpg=$4 # e.g cg12622986 
p1=$5 # 1e-04
p2=$6 # 0.98
mQTL_coloc_category=$7 # either mQTL_specific or mQTL_shared

type=$(grep $gwas $wdir/data/GWAS.metadata.txt | cut -f 21 | awk '{if($1 == 0) {print "quant"} else {print "cc"}}')
cases=$(grep $gwas $wdir/data/GWAS.metadata.txt | cut -f 22)
ss=$(grep $gwas $wdir/data/GWAS.metadata.txt | cut -f 16)
i=$(( $(echo $gwasr | sed 's/region//g')+1 ))
ldblocksFile=$wdirc/dap/data/eur_ld.hg38.bed
gwasid=$gwasr.$(cat $ldblocksFile | tail -n+2 | sed "${i}q;d" | cut -f 1-3 | sed 's/\t/_/g')
start=$(cat $ldblocksFile | tail -n+2 | sed "${i}q;d" | cut -f 2)
end=$(cat $ldblocksFile | tail -n+2 | sed "${i}q;d" | cut -f 3)

echo "$tissue $cpg $gwas $gwasr"
#if [[ ! -d $wdir/results/$tissue/$gwas/colocsummary ]]; then mkdir -p $wdir/results/$tissue/$gwas/colocsummary; fi
if [[ ! -d $wdir/results/$mQTL_coloc_category ]]; then mkdir -p $wdir/results/$mQTL_coloc_category; fi

echo "
#PBS -N perform_coloc_rand_eQTL_fastqtl
#PBS -S /bin/bash
#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=3GB
#PBS -d $wdir
#PBS -o $wdir/logs/perform_coloc_"$cpg"_"$tissue"_"$gwas"_"$gwasr"_\$PBS_JOBID.out
#PBS -e $wdir/logs/perform_coloc_"$cpg"_"$tissue"_"$gwas"_"$gwasr"_\$PBS_JOBID.err

date;
module load gcc/6.2.0 
module load bedtools
module load htslib

echo 'Consider mQTL-specific regions. Consider corresponding GWAS, GWAS region and top colocalizing (largest PP4) tissue-cpg pair'
echo 'Select mQTL locus $cpg overlapping GWAS locus $gwasid, for GWAS $gwas, region $gwasr'
echo 'For that, we intersect GWAS with mQTL locus coordinates'

i=$i
echo 'Consider region coordinates of a) mQTL locus +/- 500Kb and b) LD-defined GWAS regions'
if [[ ! -f /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.all_mQTLs.bed ]]; then ln -s /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.bed /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.all_mQTLs.bed;fi

cat $ldblocksFile | tail -n+2 | sed \"\${i}q;d\" | cut -f 1-3 | bedtools intersect -a stdin -b /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$mode.all_mQTLs.bed -wb | grep $cpg > /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$cpg.$gwas.$gwasr.bed 

echo $gwasid
if [[ ! -s /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$cpg.$gwas.$gwasr.bed ]]; then
        echo '$cpg does not overlap gwas $gwas $gwasid'
        exit
fi

if [[ ! -f /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/HyPrColoc.$tissue.$cpg.$gwas.$gwasr.mQTLs.txt ]]; then
	echo \"Extract gwas $gwas $gwasr $gwasid and mQTL data for $tissue. Merge\"
	mQTLfile=/scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/HyPrColoc.$tissue.$cpg.$gwas.$gwasr.mQTLs.txt
	> \$mQTLfile
	id=\$(cut -f 8 /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$cpg.$gwas.$gwasr.bed);
	zgrep -F $cpg /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.\$id.$mode.txt.gz >> \$mQTLfile
	bgzip -f \$mQTLfile

	# Extract GWAS stats from mQTL locus +/- 500Kb. Trim by GWAS region limits. Selected region is the intersection of mQTL locus +/- 500Kb and LD-defined GWAS region
	join --nocheck-order -1 2 -2 2 <(zcat \$mQTLfile.gz | awk '{ if((\$3^2)^(1/2) <= 500000 ) {print}}' | awk -F '_' -v start=$start -v end=$end '{if (\$2 > start && \$2 < end) {print \$0}}' | sort -k2,2) <(tail -n+2 /scratch/mpavia/GWASes/$gwas.txt) > /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/HyPrColoc.$tissue.$cpg.$gwas.$gwasr.mQTLs.txt
	rm \$mQTLfile.gz
fi

mQTL_chr=\$(head -n 1 /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/HyPrColoc.$tissue.$cpg.$gwas.$gwasr.mQTLs.txt | cut -d' ' -f 11 | sed \"s/chr//g\")
mQTL_start=\$(sort -nk12,12 /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/HyPrColoc.$tissue.$cpg.$gwas.$gwasr.mQTLs.txt | cut -d' ' -f 12 | head -n 1)
mQTL_end=\$(sort -rnk12,12 /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/HyPrColoc.$tissue.$cpg.$gwas.$gwasr.mQTLs.txt | cut -d' ' -f 12 | head -n 1)
coord=\$mQTL_chr":"\$mQTL_start"-"\$mQTL_end

for eQTLtrait in \$( cat $wdir/data/external_eQTLs.traits.txt); do
#eQTLtrait=\"TwinsUK_ge_blood\";
#eQTLtrait=\"i2QTL\";
#eQTLtrait=\"eQTLGen\";
echo \$eQTLtrait
join -1 1 -2 6 -a 1 /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/HyPrColoc.$tissue.$cpg.$gwas.$gwasr.mQTLs.txt <(tabix $wdir/data/external_eQTLs/\$eQTLtrait.all.tsv.gz \$coord | awk '{gsub(/\$/,\"_b38\",\$6); print}' | rev | cut -d' ' -f2- | rev | sort | uniq | sort -k6,6) > /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/HyPrColoc.$tissue.$cpg.$gwas.$gwasr.mQTLs.\$eQTLtrait.txt
ls /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/HyPrColoc.$tissue.$cpg.$gwas.$gwasr.mQTLs.\$eQTLtrait.txt

module load R/3.6.3
echo \"Perform HyPrColoc of gwas $gwas $gwasr and mQTL data for $tissue $cpg\"
$wdir/scripts/HyPrColoc.mQTL_specific.random_eQTL.2.R --cpg $cpg --tissue $tissue --gwas $gwas --gwasr $gwasr  --sample_size $ss --cases $cases --type $type --p1 $p1 --p2 $p2 --eQTLtrait \$eQTLtrait --coloc_category $mQTL_coloc_category
done

date;
" >  $tmpdir/logs/perform_coloc_"$cpg"_"$tissue"_"$gwas"_"$gwasr"_"$gwasid".standard_GWAS.all_mQTLs.qsub
qsub $tmpdir/logs/perform_coloc_"$cpg"_"$tissue"_"$gwas"_"$gwasr"_"$gwasid".standard_GWAS.all_mQTLs.qsub
#ls $tmpdir/logs/perform_coloc_"$cpg"_"$tissue"_"$gwas"_"$gwasr"_"$gwasid".standard_GWAS.all_mQTLs.qsub
#bash $tmpdir/logs/perform_coloc_"$cpg"_"$tissue"_"$gwas"_"$gwasr"_"$gwasid".standard_GWAS.all_mQTLs.qsub
