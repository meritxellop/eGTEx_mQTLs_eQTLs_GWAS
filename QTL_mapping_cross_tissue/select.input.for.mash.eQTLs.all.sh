#PBS -N perform_coloc_fastqtl
#PBS -S /bin/bash
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=25GB
#PBS -d /gpfs/data/gtex-group/mpavia/methylation/mashr
#PBS -o /scratch/mpavia/fastQTL/logs/mash_\$PBS_JOBID.out
#PBS -e /scratch/mpavia/fastQTL/logs/mash_\$PBS_JOBID.err

wdir=/gpfs/data/gtex-group/mpavia/methylation
TISSUE=$(echo $tissue | sed 's/_//g')
join --nocheck-order -a 1 <(sort /gpfs/data/gtex-group/mpavia/methylation/results/permuted/signif.top_eQTL.fdr005.anytissue.all.top.txt) <(awk -F' ' 'NR==FNR{a[$1]=1;} NR!=FNR{ if(a[$1] == 1 ) {print $0}}' /gpfs/data/gtex-group/mpavia/methylation/results/permuted/signif.top_eQTL.fdr005.anytissue.all.top.txt <(zcat /gpfs/data/gtex-group/v8/59349/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_all_associations/$tissue.allpairs.txt.gz | awk '{print $1":"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}') | sort -k1,1) > $wdir/mashr/$TISSUE.eQTL_top.fdr005.all.txt
