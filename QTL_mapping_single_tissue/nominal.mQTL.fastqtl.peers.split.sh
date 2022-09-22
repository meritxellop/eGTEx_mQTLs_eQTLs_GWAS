module load gcc/6.2.0 fastqtl/2.184_gtex bcftools
module load htslib

tissue=$1 # e.g. Lung
id=$2 # e.g. 0000
mode=$3 # e.g. regular

wdir=/gpfs/data/gtex-group/mpavia/methylation
tmpdir=/scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults
tmpdir2=/scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles
vcfFile=$wdir/data/Genotypes/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz
methFile=$tmpdir2/"$tissue"_"$id".gz
tabix -p bed -f $methFile
covFile=$wdir/data/Covariates/$tissue.covariates.txt
outFileTmp=$tmpdir/tmpFiles/$tissue.$id.$mode.txt
sampleFile=$wdir/data/support_files/$tissue.mQTL.samples.txt
thresholds=$wdir/results/permuted/$tissue.regular.perm.fdr005.thresholds.txt
thresholds_phenotypes=$wdir/results/permuted/$tissue.regular.perm.fdr005.thresholds.phenotypes.txt

if [ ! -d $tmpdir2 ]; then
	mkdir -p $tmpdir/tmpFiles
fi
nchr=$(zcat $methFile | cut -f 1 | uniq | wc -l)
if [[ $mode = "regular" ]]; then
	echo "Chunk $id of $tissue"
        if [ $nchr -lt 3 ]; then
                echo "Chunk $id of $tissue"
                fastQTL --vcf $vcfFile --bed $methFile --cov $covFile --out $outFileTmp --chunk 1 1 --include-samples $sampleFile --maf-threshold 0.01 -W 500000
		gzip -f $outFileTmp;
        else
                outFileTmp=$tmpdir/tmpFiles/$tissue.$id.$mode.1.txt
                fastQTL --vcf $vcfFile --bed $methFile --cov $covFile --out $outFileTmp --chunk 1 2 --include-samples $sampleFile --maf-threshold 0.01 -W 500000
                outFileTmp=$tmpdir/tmpFiles/$tissue.$id.$mode.2.txt
                fastQTL --vcf $vcfFile --bed $methFile --cov $covFile --out $outFileTmp --chunk 2 2 --include-samples $sampleFile --maf-threshold 0.01 -W 500000
		outFileTmp=$tmpdir/tmpFiles/$tissue.$id.$mode.txt
		cat $tmpdir/tmpFiles/$tissue.$id.$mode.1.txt $tmpdir/tmpFiles/$tissue.$id.$mode.2.txt > $outFileTmp
		rm $tmpdir/tmpFiles/$tissue.$id.$mode.1.txt $tmpdir/tmpFiles/$tissue.$id.$mode.2.txt
		gzip -f $outFileTmp;
        fi
fi
if [[ $mode = "conditional" ]]; then
        echo "Chunk $id of $tissue"
        if [ $nchr -lt 3 ]; then
                echo "Chunk $id of $tissue"
                fastQTL --vcf $vcfFile --bed $methFile --cov $covFile --out $outFileTmp --chunk 1 1 --include-samples $sampleFile --maf-threshold 0.01 -W 500000 --map $thresholds --include-phenotypes $thresholds_phenotypes 
		gzip -f $outFileTmp;
        else
                outFileTmp=$tmpdir/tmpFiles/$tissue.$id.$mode.1.txt
                fastQTL --vcf $vcfFile --bed $methFile --cov $covFile --out $outFileTmp --chunk 1 2 --include-samples $sampleFile --maf-threshold 0.01 -W 500000 --map $thresholds --include-phenotypes $thresholds_phenotypes
                outFileTmp=$tmpdir/tmpFiles/$tissue.$id.$mode.2.txt
                fastQTL --vcf $vcfFile --bed $methFile --cov $covFile --out $outFileTmp --chunk 2 2 --include-samples $sampleFile --maf-threshold 0.01 -W 500000 --map $thresholds --include-phenotypes $thresholds_phenotypes
		outFileTmp=$tmpdir/tmpFiles/$tissue.$id.$mode.txt
                cat $tmpdir/tmpFiles/$tissue.$id.$mode.1.txt $tmpdir/tmpFiles/$tissue.$id.$mode.2.txt > $outFileTmp
		rm $tmpdir/tmpFiles/$tissue.$id.$mode.1.txt $tmpdir/tmpFiles/$tissue.$id.$mode.2.txt
		gzip -f $outFileTmp;
        fi
fi
if [[ $mode = "interaction" ]]; then
	interaction=$4 # e.g. sex [optional]
	outFileTmp=$tmpdir/tmpFiles/$tissue.$id.$mode.$interaction.txt
	covFile=$wdir/data/Covariates/interaction_Covariates/$interaction/$tissue.covariates.no$interaction.txt
	interactionFile=$wdir/data/Covariates/interaction_Covariates/$interaction/$tissue.$interaction.txt
        if [ $nchr -lt 3 ]; then
		outFileTopTmp=$tmpdir/tmpFiles/$tissue.$id.$mode.$interaction.topSNP.txt
                echo "Chunk $id of $tissue"
                fastQTL --vcf $vcfFile --bed $methFile --cov $covFile --out $outFileTmp --chunk 1 1 --include-samples $sampleFile --maf-threshold 0.01 -W 500000 --interaction $interactionFile --report-best-only 
		cat $outFileTmp | awk ' NR==1 { gene=$1;genepval[gene]=$7;line=$0 } { if ($1 == gene) { if ($7 < genepval[gene]) { genepval[gene]=$7; line=$0 } } else { print line; gene=$1; genepval[gene]=$7; line=$0 } } END { print line }' > $outFileTopTmp
		gzip -f $outFileTmp;
        else
                echo "Chunk $id of $tissue part1"
		outFileTmp=$tmpdir/tmpFiles/$tissue.$id.$mode.$interaction.1.txt
		outFileTopTmp=$tmpdir/tmpFiles/$tissue.$id.$mode.$interaction.topSNP.1.txt
                fastQTL --vcf $vcfFile --bed $methFile --cov $covFile --out $outFileTmp --chunk 1 2 --include-samples $sampleFile --maf-threshold 0.01 -W 500000 --interaction $interactionFile --report-best-only
		cat $outFileTmp | awk ' NR==1 { gene=$1;genepval[gene]=$7;line=$0 } { if ($1 == gene) { if ($7 < genepval[gene]) { genepval[gene]=$7; line=$0 } } else { print line; gene=$1; genepval[gene]=$7; line=$0 } } END { print line }' > $outFileTopTmp
		gzip -f $outFileTmp;
                echo "Chunk $id of $tissue part2"
		outFileTmp=$tmpdir/tmpFiles/$tissue.$id.$mode.$interaction.2.txt
                outFileTopTmp=$tmpdir/tmpFiles/$tissue.$id.$mode.$interaction.topSNP.2.txt
                fastQTL --vcf $vcfFile --bed $methFile --cov $covFile --out $outFileTmp --chunk 2 2 --include-samples $sampleFile --maf-threshold 0.01 -W 500000 --interaction $interactionFile --report-best-only
		cat $outFileTmp | awk ' NR==1 { gene=$1;genepval[gene]=$7;line=$0 } { if ($1 == gene) { if ($7 < genepval[gene]) { genepval[gene]=$7; line=$0 } } else { print line; gene=$1; genepval[gene]=$7; line=$0 } } END { print line }' > $outFileTopTmp
		gzip -f $outFileTmp;
        fi
fi

#cat <(echo -e "SNP\tgene\tbeta\tt-stat\tp-value"; awk -F "\t" '{ if ($9 = 0 || $8 = 0) {print $2"\t"$1"\t"$8"\tNA\t"$7} else {print $2"\t"$1"\t"$8"\t"$8/$9"\t"$7}}' $outFileTmp) > $outFile
#ln -s $outFileTmp $outFile
#gzip -f $outFile;
#rm -rf $tmpdir/tmpFiles/$outFileTmp;
