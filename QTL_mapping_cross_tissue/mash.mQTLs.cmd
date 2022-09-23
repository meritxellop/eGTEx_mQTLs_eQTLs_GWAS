#PBS -N perform_mash_fastqtl
#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=50GB
#PBS -d /gpfs/data/gtex-group/mpavia/methylation/mashr
#PBS -o /gpfs/data/gtex-group/mpavia/methylation/mashr/logs/perform_mash_mQTL_\$PBS_JOBID.out
#PBS -e /gpfs/data/gtex-group/mpavia/methylation/mashr/logs/perform_mash_mQTL_\$PBS_JOBID.err

date;
module load gcc/6.2.0
module load mosek/8.1
module load R/3.6.3

# Allocate 3GB, 1h45m-2h for this step
R --vanilla <  /gpfs/data/gtex-group/mpavia/methylation/mashr/scripts/mash.from.betas.to.FLASH.mQTLs.R > logs/mash.from.betas.to.FLASH.mQTLs.Rout

# Allocate 10GB, 40m for this step
R --vanilla <  /gpfs/data/gtex-group/mpavia/methylation/mashr/scripts/mash.from.FLASH.to.mash.model.mQTLs.R > logs/mash.from.FLASH.to.mash.model.mQTLs.Rout
R --vanilla <  /gpfs/data/gtex-group/mpavia/methylation/mashr/scripts/mash.from.Ulist.to.mash.model.mQTLs.R > logs/mash.from.Ulist.to.mash.model.mQTLs.Rout

# Allocate XGB, Xm for this step
R --vanilla <  /gpfs/data/gtex-group/mpavia/methylation/mashr/scripts/mash.from.mash.model.to.posteriors.mQTLs.R > mash_output/mash.from.mash.model.to.posteriors.mQTLs.Rout
