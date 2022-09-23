wdir=$PWD;
qtl=$1; # e.g. eqtls/mqtls
feature=$2; # e.g. Promoter 
feature_category=$3 # e.g. ENSEMBL_REGULATORY_BUILD 
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
#PBS -l mem=20GB
#PBS -d $wdir
#PBS -o $wdir/logs/generate_torus_results_\$PBS_JOBID.out
#PBS -e $wdir/logs/generate_torus_results_\$PBS_JOBID.err

date;
module load torus

echo \"Run torus for $feature_category $feature in $tissue\"
if [ $feature_category == 'ENSEMBL_REGULATORY_BUILD' ]; then
#torus -d $egene_allSNPGenepairs_file --fastqtl -est -annot $wdir/torus_annots/annotations.$feature.txt.gz > $wdir/output/$feature.$tissue.est.txt;
#torus -d $egene_allSNPGenepairs_file --fastqtl -est -annot $wdir/torus_annots/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.complemented.txt.gz > $wdir/output/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.$tissue.$qtl.est.txt;
torus -d $egene_allSNPGenepairs_file --fastqtl -est -annot $wdir/torus_annots/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.complemented.full.txt.gz > $wdir/output/full/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.$tissue.$qtl.est.txt;
elif [ $feature_category == 'JASPAR2018' ]; then
#torus -d $egene_allSNPGenepairs_file --fastqtl -est -annot $wdir/torus_annots/$feature.jaspar_score400.annot.txt.gz > $wdir/output/$feature.$tissue.est.txt;
torus -d $egene_allSNPGenepairs_file --fastqtl -est -annot $wdir/torus_annots/$feature.jaspar_scorefull.annot.txt.gz > $wdir/output/$feature.$tissue.scorefull.est.txt;
#torus -d $egene_allSNPGenepairs_file --fastqtl -est -annot $wdir/torus_annots/$feature.gfr.jaspar_scorefull.annot.txt.gz > $wdir/output/$feature.gfr.$tissue.scorefull.est.txt;
fi

date;
" >  $wdir/logs/generate_torus_results_"$feature_category"."$feature".$tissue.$qtl.qsub
qsub $wdir/logs/generate_torus_results_"$feature_category"."$feature".$tissue.$qtl.qsub
