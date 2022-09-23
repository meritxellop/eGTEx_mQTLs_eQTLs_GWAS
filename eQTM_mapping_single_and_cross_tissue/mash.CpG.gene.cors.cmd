#PBS -N perform_mash_CpG_gene_cors
#PBS -S /bin/bash
#PBS -l walltime=04:15:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=2GB
#PBS -d /gpfs/data/gtex-group/mpavia/methylation/CpG_gene_correlation
#PBS -o /gpfs/data/gtex-group/mpavia/methylation/CpG_gene_correlation/logs/perform_mash_CpG_gene_cors_\$PBS_JOBID.out
#PBS -e /gpfs/data/gtex-group/mpavia/methylation/CpG_gene_correlation/logs/perform_mash_CpG_gene_cors_\$PBS_JOBID.err

date;
module load gcc/6.2.0
module load mosek/8.1
module load R/3.6.3

# Allocate 1.5GB, 15m for this step
R --vanilla <  /gpfs/data/gtex-group/mpavia/methylation/CpG_gene_correlation/scripts/mash.from.betas.to.FLASH.CpG.gene.cors.R > logs/mash.from.betas.to.FLASH.CpG.gene.cors.Rout

# Allocate 10GB, 45min-1h for this step
R --vanilla <  /gpfs/data/gtex-group/mpavia/methylation/CpG_gene_correlation/scripts/mash.from.FLASH.to.mash.model.CpG.gene.cors.R > logs/mash.from.FLASH.to.mash.model.CpG.gene.cors.Rout
# Do not execute the following line
#R --vanilla < /gpfs/data/gtex-group/mpavia/methylation/CpG_gene_correlation/scripts/mash.from.Ulist.to.mash.model.CpG.gene.cors.R > logs/mash.from.Ulist.to.mash.model.CpG.gene.cors.Rout

# Allocate 1.5GB, 1h for this step
R --vanilla <  /gpfs/data/gtex-group/mpavia/methylation/CpG_gene_correlation/scripts/mash.from.mash.model.to.posteriors.CpG.gene.cors.R > logs/mash.from.mash.model.to.posteriors.CpG.gene.cors.Rout
