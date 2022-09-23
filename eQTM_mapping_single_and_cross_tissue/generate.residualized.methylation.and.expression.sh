wdir=/gpfs/data/gtex-group/mpavia/methylation/CpG_gene_correlation
for tissue in $(cat ../coloc/data/tissues.txt); do 
	echo $tissue
	TISSUE=$(echo $tissue | sed 's/_//g');
echo "
#PBS -N generate_residualized_met_geneexp
#PBS -S /bin/bash
#PBS -l walltime=00:15:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=10GB
#PBS -d $wdir
#PBS -o $wdir/logs/generate_residualized_met_geneexp_"$tissue"_\$PBS_JOBID.out
#PBS -e $wdir/logs/generate_residualized_met_geneexp_"$tissue"_\$PBS_JOBID.err

date;
Rscript --vanilla /gpfs/data/gtex-group/mpavia/methylation/CpG_gene_correlation/scripts/generate.residualized.methylation.and.expression.R $tissue $TISSUE
date;
" >  $wdir/logs/generate_residualized_met_geneexp_"$tissue".qsub
qsub $wdir/logs/generate_residualized_met_geneexp_"$tissue".qsub
done

