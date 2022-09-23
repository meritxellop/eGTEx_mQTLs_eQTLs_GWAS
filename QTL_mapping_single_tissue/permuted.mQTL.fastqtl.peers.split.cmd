wdir='/gpfs/data/gtex-group/mpavia/methylation';
tmpdir='/scratch/mpavia/fastQTL'

tissue=$1 # e.g. Lung
id=$2 # e.g. 0000
npeers=$3 # 5
mode=$4 # regular
echo "
#PBS -N generate_permuted_gtex_fastqtl
#PBS -S /bin/bash
#PBS -l walltime=15:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=100MB
#PBS -d $wdir
#PBS -o $tmpdir/logs/generate_permuted_gtex_"$tissue"_"$npeers"_"$mode"_\$PBS_JOBID.out
#PBS -e $tmpdir/logs/generate_permuted_gtex_"$tissue"_"$npeers"_"$mode"_\$PBS_JOBID.err

date;
bash $wdir/scripts/permuted.mQTL.fastqtl.peers.split.sh $tissue $id $npeers $mode
date;
" >  $tmpdir/logs/generate_permuted_gtex_"$tissue"_"$id"_"$npeers"_"$mode".qsub
qsub $tmpdir/logs/generate_permuted_gtex_"$tissue"_"$id"_"$npeers"_"$mode".qsub
#bash $tmpdir/logs/generate_permuted_gtex_"$tissue"_"$id"_"$npeers"_"$mode".qsub
