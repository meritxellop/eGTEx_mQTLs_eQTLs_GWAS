#PBS -N extract_correlation
#PBS -S /bin/bash
#PBS -l walltime=01:30:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=5GB
#PBS -d /gpfs/data/gtex-group/mpavia/methylation/CpG_gene_correlation
#PBS -o /gpfs/data/gtex-group/mpavia/methylation/CpG_gene_correlation/logs/extract_correlation_\$PBS_JOBID.out
#PBS -e /gpfs/data/gtex-group/mpavia/methylation/CpG_gene_correlation/logs/extract_correlation_\$PBS_JOBID.err

date;
Rscript --vanilla /gpfs/data/gtex-group/mpavia/methylation/CpG_gene_correlation/scripts/extract.cor.stats.per.block.R num tissue > logs/extract.cor.stats.per.block.num.tissue.Rout
