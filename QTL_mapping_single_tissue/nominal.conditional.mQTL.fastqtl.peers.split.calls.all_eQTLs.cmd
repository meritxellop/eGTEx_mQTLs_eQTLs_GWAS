wdir='/gpfs/data/gtex-group/mpavia/methylation';
tmpdir='/scratch/mpavia/fastQTL'

tissue=$1 # e.g. Lung
if [[ $tissue =~ .*_.* ]] ; then TISSUE=$(echo $tissue | sed 's/_//g'); else TISSUE=$tissue; fi
id=$2 # e.g. 0000
if [[ ! -s /scratch/mpavia/fastQTL/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles/$TISSUE.$id.bwd.all_eQTLs.ALL_genes.txt && ! -s /scratch/mpavia/fastQTL/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles/$TISSUE.$id.bwd.all_eQTLs.all_genes.txt ]]; then
echo $tissue $id
echo "
#PBS -N generate_nominal_gtex_fastqtl
#PBS -S /bin/bash
#PBS -l walltime=05:59:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8GB
#PBS -d $wdir
#PBS -o $tmpdir/logs/generate_nominal_gtex_"$tissue"_\$PBS_JOBID.out
#PBS -e $tmpdir/logs/generate_nominal_gtex_"$tissue"_\$PBS_JOBID.err

date;
bash $wdir/scripts/nominal.conditional.mQTL.fastqtl.peers.split.calls.all_eQTLs.sh $tissue $id
date;
" >  $tmpdir/logs/generate_nominal_gtex_"$tissue"_"$id".qsub
#qsub $tmpdir/logs/generate_nominal_gtex_"$tissue"_"$id".qsub
ls $tmpdir/logs/generate_nominal_gtex_"$tissue"_"$id".qsub
#bash $tmpdir/logs/generate_nominal_gtex_"$tissue"_"$id".qsub
fi
