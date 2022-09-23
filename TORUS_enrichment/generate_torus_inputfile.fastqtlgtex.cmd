wdir=$PWD
qtl=$1 # e.g. mqtls/eqtls

for tissue in $(cat ../../coloc/data/tissues.txt); do 
	echo $tissue
echo "
#PBS -N generate_torus_eqtl_format
#PBS -S /bin/bash
#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=5GB
#PBS -d $wdir
#PBS -o $wdir/logs/generate_torus_eqtl_format_"$tissue"_"$qtl"_\$PBS_JOBID.out
#PBS -e $wdir/logs/generate_torus_eqtl_format_"$tissue"_"$qtl"_\$PBS_JOBID.err

date

bash scripts/generate_torus_inputfile.fastqtlgtex.sh $tissue $qtl

date
" >  $wdir/logs/generate_torus_eqtl_format_"$tissue".$qtl.qsub
qsub $wdir/logs/generate_torus_eqtl_format_"$tissue".$qtl.qsub
done
