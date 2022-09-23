#PBS -N perform_mash_fastqtl
#PBS -S /bin/bash
#PBS -l walltime=5:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=10GB
#PBS -d /gpfs/data/gtex-group/mpavia/methylation/mashr
#PBS -o /gpfs/data/gtex-group/mpavia/methylation/mashr/logs/perform_mash_eQTL_\$PBS_JOBID.out
#PBS -e /gpfs/data/gtex-group/mpavia/methylation/mashr/logs/perform_mash_eQTL_\$PBS_JOBID.err

date;
module load gcc/6.2.0
module load mosek/8.1
module load R/3.6.3

# Allocate 1.5GB, 15m for this step
#R --vanilla <  /gpfs/data/gtex-group/mpavia/methylation/mashr/scripts/mash.from.betas.to.FLASH.eQTLs.R > logs/mash.from.betas.to.FLASH.eQTLs.Rout

# Allocate 10GB, 45min-1h for this step
#R --vanilla <  /gpfs/data/gtex-group/mpavia/methylation/mashr/scripts/mash.from.FLASH.to.mash.model.eQTLs.R > logs/mash.from.FLASH.to.mash.model.eQTLs.Rout

# Allocate 1.5GB, 15m for this step
R --vanilla <  /gpfs/data/gtex-group/mpavia/methylation/mashr/scripts/mash.from.mash.model.to.posteriors.eQTLs.R > logs/mash.from.mash.model.to.posteriors.eQTLs.Rout
