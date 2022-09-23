wdir=$PWD
tmpdir=/scratch/mpavia/colocalization
vcfFile=/gpfs/data/gtex-group/mpavia/methylation/data/Genotypes/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz
### Organize CpG PIPs by tissue-chr
for tissue in $(cat /gpfs/data/gtex-group/mpavia/methylation/data/support_files/tissues.txt); do
	TISSUE=$tissue
	if [[ $tissue =~ .*_.* ]] ; then tissue=$(echo $tissue | sed 's/_//g'); fi
	eqtlFile=/gpfs/data/gtex-group/mpavia/methylation/coloc/dap/data/gtex_v8.eqtl_annot.vcf.gz
	# pre-computed eQTL annotations from GTEx (v8) data were obtained from the fastenloc github repository: https://github.com/xqwen/fastenloc
	if [[ ! -d /scratch/mpavia/colocalization/results/$tissue/fastenloc ]]; then
		mkdir -p /scratch/mpavia/colocalization/results/$tissue/fastenloc
	fi
	for gwas in $(cat $wdir/data/GWASes.87.txt); do
		if [[ -f /scratch/mpavia/colocalization/results/$tissue/fastenloc/$tissue.$gwas.all.eQTLs.enloc.enrich.out ]]; then continue; fi
		echo $tissue $gwas TODO

echo "
#PBS -N perform_coloc_fastqtl
#PBS -S /bin/bash
#PBS -l walltime=15:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=10GB
#PBS -d $wdir
#PBS -o $tmpdir/logs/run_fastenloc_eQTLs_"$tissue"_"$gwas"_\$PBS_JOBID.out
#PBS -e $tmpdir/logs/run_fastenloc_eQTLs_"$tissue"_"$gwas"_\$PBS_JOBID.err
date
module load gcc/6.2.0
module load fastenloc

gwasFile=/scratch/mpavia/GWASes/$gwas.eur_ld.hg38.pip.gz
if [[ ! -f \$gwasFile ]]; then
	module load htslib
	echo \"Slice GWAS\"
	module load python/3.8.1
        python3 /gpfs/data/im-lab/nas40t2/abarbeira/software/genomic_tools/genomic_tools/src/slice_gwas_by_region.py -region_file $wdir/data/eur_ld.hg38.bed.gz -gwas_file /gpfs/data/im-lab/nas40t2/Data/SummaryResults/imputed_gwas_hg38_1.1/imputed_"$gwas".txt.gz  -output /scratch/mpavia/GWASes/$gwas.eur_ld.hg38.torus.zval
	bgzip /scratch/mpavia/GWASes/$gwas.eur_ld.hg38.torus.zval
	echo \"Compute GWAS PIPs\"
	module load torus
	torus -d /scratch/mpavia/GWASes/$gwas.eur_ld.hg38.torus.zval.gz --load_zval -dump_pip /scratch/mpavia/GWASes/$gwas.eur_ld.hg38.pip
	cut -f1,2,4 /scratch/mpavia/GWASes/$gwas.eur_ld.hg38.pip > /scratch/mpavia/GWASes/$gwas.eur_ld.hg38.pip.tmp
	mv /scratch/mpavia/GWASes/$gwas.eur_ld.hg38.pip.tmp /scratch/mpavia/GWASes/$gwas.eur_ld.hg38.pip
	bgzip /scratch/mpavia/GWASes/$gwas.eur_ld.hg38.pip
fi

#/gpfs/data/gtex-group/mpavia/methylation/coloc/dap/src/fastenloc.static -eqtl $eqtlFile -gwas \$gwasFile -t $tissue -thread 1 --all -prefix /scratch/mpavia/colocalization/results/$tissue/fastenloc/$tissue.$gwas.all
# the static version does not have the '--all' param implemented
fastenloc -eqtl $eqtlFile -gwas \$gwasFile -t $TISSUE -thread 1 --all -prefix /scratch/mpavia/colocalization/results/$tissue/fastenloc/$tissue.$gwas.all.eQTLs

date
" > $tmpdir/logs/run_fastenloc_eQTLs_"$tissue"_"$gwas".eQTLs.qsub
qsub $tmpdir/logs/run_fastenloc_eQTLs_"$tissue"_"$gwas".eQTLs.qsub;
	done;
done
