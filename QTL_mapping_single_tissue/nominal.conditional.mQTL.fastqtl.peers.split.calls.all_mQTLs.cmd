wdir='/gpfs/data/gtex-group/mpavia/methylation';
tmpdir='/scratch/mpavia/fastQTL'

tissue=$1 # e.g. Lung
id=$2 # e.g. 0000
if [[ ! -s /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$id.bwd.all_mQTLs.ALL_cpgs.txt && ! -s /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$id.bwd.all_mQTLs.all_cpgs.txt ]]; then
echo "
#PBS -N generate_nominal_gtex_fastqtl
#PBS -S /bin/bash
#PBS -l walltime=06:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=200MB
#PBS -d $wdir
#PBS -o $tmpdir/logs/generate_nominal_gtex_"$tissue"_\$PBS_JOBID.out
#PBS -e $tmpdir/logs/generate_nominal_gtex_"$tissue"_\$PBS_JOBID.err

date;
bash $wdir/scripts/nominal.conditional.mQTL.fastqtl.peers.split.calls.all_mQTLs.sh $tissue $id
date;
" >  $tmpdir/logs/generate_nominal_gtex_"$tissue"_"$id".qsub
#qsub $tmpdir/logs/generate_nominal_gtex_"$tissue"_"$id".qsub
bash $tmpdir/logs/generate_nominal_gtex_"$tissue"_"$id".qsub
#ls $tmpdir/logs/generate_nominal_gtex_"$tissue"_"$id".qsub
fi
