module load gcc/6.2.0 fastqtl/2.184_gtex bcftools
module load htslib

tissue=$1 # e.g. Lung
id=$2 # e.g. 0000
npeers=$3 # e.g. 5
mode=$4 # regular

wdir=/gpfs/data/gtex-group/mpavia/methylation
tmpdir=/scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults
tmpdir2=/scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles
vcfFile=$wdir/data/Genotypes/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz
methFile=$tmpdir2/"$tissue"_"$id".gz
#covFile=$wdir/data/Covariates/$tissue.npeers.$npeers.blockcovs.final.txt
covFile=$wdir/data/Covariates/$tissue.$tissue.covariates.txt
outFileTmp=$tmpdir/tmpFiles/$tissue.$id.$mode.$npeers.perm.txt
sampleFile=$wdir/data/support_files/$tissue.mQTL.samples.txt

if [ ! -d $tmpdir2 ]; then
	mkdir -p $tmpdir/tmpFiles
fi
nchr=$(zcat $methFile | cut -f 1 | uniq | wc -l)
if [[ $mode = "regular" ]]; then
	if [ $nchr -lt 3 ]; then
		echo "Chunk $id of $tissue"
		fastQTL --vcf $vcfFile --bed $methFile --cov $covFile --out $outFileTmp --chunk 1 1 --include-samples $sampleFile -P 1000 --maf-threshold 0.01 -W 500000
	else
		echo "Chunk $id of $tissue part1"
		outFileTmp=$tmpdir/tmpFiles/$tissue.$id.$mode.$npeers.perm.1.txt
		fastQTL --vcf $vcfFile --bed $methFile --cov $covFile --out $outFileTmp --chunk 1 2 --include-samples $sampleFile -P 1000 --maf-threshold 0.01 -W 500000
		echo "Chunk $id of $tissue part2"
		outFileTmp=$tmpdir/tmpFiles/$tissue.$id.$mode.$npeers.perm.2.txt
		fastQTL --vcf $vcfFile --bed $methFile --cov $covFile --out $outFileTmp --chunk 2 2 --include-samples $sampleFile -P 1000 --maf-threshold 0.01 -W 500000
		cat $tmpdir/tmpFiles/$tissue.$id.$mode.$npeers.perm.1.txt $tmpdir/tmpFiles/$tissue.$id.$mode.$npeers.perm.2.txt > $tmpdir/tmpFiles/$tissue.$id.$mode.$npeers.perm.txt
		rm $tmpdir/tmpFiles/$tissue.$id.$mode.$npeers.perm.1.txt $tmpdir/tmpFiles/$tissue.$id.$mode.$npeers.perm.*.txt
	fi
fi
if [[ $mode = "interaction" ]]; then
	covFile=$wdir/data/Covariates/$tissue.npeers.$npeers.blockcovs.final.nosex.txt
	interactionFile=$wdir/data/Covariates/$tissue.sex.txt
        if [ $nchr -lt 3 ]; then
                echo "Chunk $id of $tissue"
                fastQTL --vcf $vcfFile --bed $methFile --cov $covFile --out $outFileTmp --chunk 1 1 --include-samples $sampleFile -P 1000 --maf-threshold 0.01 -W 500000 --interaction $interactionFile 
        else
                echo "Chunk $id of $tissue part1"
                outFileTmp=$tmpdir/tmpFiles/$tissue.$id.$mode.$npeers.perm.1.txt
                fastQTL --vcf $vcfFile --bed $methFile --cov $covFile --out $outFileTmp --chunk 1 2 --include-samples $sampleFile -P 1000 --maf-threshold 0.01 -W 500000 --interaction $interactionFile
                echo "Chunk $id of $tissue part2"
                outFileTmp=$tmpdir/tmpFiles/$tissue.$id.$mode.$npeers.perm.2.txt
                fastQTL --vcf $vcfFile --bed $methFile --cov $covFile --out $outFileTmp --chunk 2 2 --include-samples $sampleFile -P 1000 --maf-threshold 0.01 -W 500000 --interaction $interactionFile
		cat $tmpdir/tmpFiles/$tissue.$id.$mode.$npeers.perm.1.txt $tmpdir/tmpFiles/$tissue.$id.$mode.$npeers.perm.2.txt > $tmpdir/tmpFiles/$tissue.$id.$mode.$npeers.perm.txt
		rm $tmpdir/tmpFiles/$tissue.$id.$mode.$npeers.perm.1.txt $tmpdir/tmpFiles/$tissue.$id.$mode.$npeers.perm.*.txt
        fi
fi

#cat <(echo -e "SNP\tgene\tbeta\tt-stat\tp-value"; awk -F "\t" '{ if ($9 = 0 || $8 = 0) {print $2"\t"$1"\t"$8"\tNA\t"$7} else {print $2"\t"$1"\t"$8"\t"$8/$9"\t"$7}}' $outFileTmp) > $outFile
#ln -s $outFileTmp $outFile
#gzip -f $outFile;
#rm -rf $tmpdir/tmpFiles/$outFileTmp;
