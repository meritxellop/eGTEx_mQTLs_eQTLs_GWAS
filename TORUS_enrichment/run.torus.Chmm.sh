wdir=$PWD;
qtl=$1; # e.g. eqtls/mqtls
feature=$2; # e.g. Chmm
feature_category=$3 # e.g. Chmm 
tissue=$4 # e.g. Lung
egene_allSNPGenepairs_file=$wdir/data/$qtl/$tissue.allpairs_4torus.txt.gz
if [[ $qtl == 'mqtls' ]]; then egene_allSNPGenepairs_file=$wdir/data/$qtl/$tissue.allpairs_4torus.lite.txt.gz; fi

module load torus

echo "Calculate annotations enrichment, conditional on tss distance"
echo $tissue

echo "
#PBS -N generate_torus_results_torus
#PBS -S /bin/bash
#PBS -l walltime=25:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=50GB
#PBS -d $wdir
#PBS -o $wdir/logs/generate_torus_results_\$PBS_JOBID.out
#PBS -e $wdir/logs/generate_torus_results_\$PBS_JOBID.err

date;
module load torus

echo \"Run torus for $feature_category $feature in $tissue\"
torus -d $egene_allSNPGenepairs_file --fastqtl -est -annot $wdir/torus_annots/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.Chmm.active_states.$tissue.txt.gz > $wdir/output/Chmm/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.$tissue.$qtl.est.txt;

date;
" >  $wdir/logs/generate_torus_results_"$feature_category"."$feature".$tissue.$qtl.qsub
qsub $wdir/logs/generate_torus_results_"$feature_category"."$feature".$tissue.$qtl.qsub
