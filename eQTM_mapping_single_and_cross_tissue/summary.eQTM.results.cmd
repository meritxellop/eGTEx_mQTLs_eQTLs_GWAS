#PBS -N perform_mash_fastqtl
#PBS -S /bin/bash
#PBS -l walltime=01:30:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=15GB
#PBS -d /gpfs/data/gtex-group/mpavia/methylation/CpG_gene_correlation
#PBS -o /gpfs/data/gtex-group/mpavia/methylation/CpG_gene_correlation/logs/summary_cor_results\$PBS_JOBID.out
#PBS -e /gpfs/data/gtex-group/mpavia/methylation/CpG_gene_correlation/logs/summary_cor_results\$PBS_JOBID.err

date;
module load R/3.6.3
Rscript --vanilla /gpfs/data/gtex-group/mpavia/methylation/CpG_gene_correlation/scripts/summary.eQTM.results.R TISSUE > logs/summary.eQTM.results.TISSUE.Rout
date;
