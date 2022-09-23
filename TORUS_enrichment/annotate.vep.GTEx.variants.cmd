wdir=$PWD
for chr in $(echo $(seq 1 22) X); do
chr=chr$chr
echo "
#PBS -N annotate_vep
#PBS -S /bin/bash
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=10GB
#PBS -d $wdir
#PBS -o $wdir/logs/annotate_vep_"$chr"_\$PBS_JOBID.out
#PBS -e $wdir/logs/annotate_vep_"$chr"_\$PBS_JOBID.err
module load gcc/6.2.0
module load ensembl-vep/100.0
module load htslib

# curl -O ftp://ftp.ensembl.org/pub/release-102/variation/indexed_vep_cache/homo_sapiens_vep_102_GRCh38.tar.gz
# perl /apps/software/gcc-6.2.0/ensembl-vep/100.0/convert_cache.pl --dir /gpfs/data/gtex-group/mpavia/methylation/post_coloc_analyses/TORUS/data --species homo_sapiens --version 102

echo 'Subsetting vcf for' $chr
tabix /gpfs/data/gtex-group/sex_biased_regulation_v8/data/Genotype_Data/raw/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz $chr > /scratch/mpavia/TORUS/$chr.vcf
echo 'Annotating...'
vep -i /scratch/mpavia/TORUS/$chr.vcf --offline --dir_cache /gpfs/data/gtex-group/mpavia/methylation/post_coloc_analyses/TORUS/data -o /scratch/mpavia/TORUS/WGS_VEP_4torus.MAF01.$chr --cache_version 102 --force_overwrite
echo 'Compressing...'
bgzip /scratch/mpavia/TORUS/WGS_VEP_4torus.MAF01.$chr
echo 'Done'
rm /scratch/mpavia/TORUS/$chr.vcf" > logs/vep.$chr.qsub
qsub logs/vep.$chr.qsub
done
