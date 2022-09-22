wdir='/gpfs/data/gtex-group/mpavia/methylation';
tmpdir='/scratch/mpavia/fastQTL'

tissue=$1 # e.g. Lung
id=$2 # e.g. 0000
mode=$3 # e.g. regular
interaction=$4 # e.g. sex [optional]
echo "
#PBS -N generate_nominal_gtex_fastqtl
#PBS -S /bin/bash
#PBS -l walltime=15:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=100MB
#PBS -d $wdir
#PBS -o $tmpdir/logs/generate_nominal_gtex_"$tissue"_"$interaction"_"$mode"_\$PBS_JOBID.out
#PBS -e $tmpdir/logs/generate_nominal_gtex_"$tissue"_"$interaction"_"$mode"_\$PBS_JOBID.err

date;
bash $wdir/scripts/nominal.mQTL.fastqtl.peers.split.sh $tissue $id $mode $interaction
date;
" >  $tmpdir/logs/generate_nominal_gtex_"$tissue"_"$id"_"$interaction"_"$mode".qsub
qsub $tmpdir/logs/generate_nominal_gtex_"$tissue"_"$id"_"$interaction"_"$mode".qsub
#bash $tmpdir/logs/generate_nominal_gtex_"$tissue"_"$id"_"$interaction"_"$mode".qsub
