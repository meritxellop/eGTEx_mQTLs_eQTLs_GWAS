wdir=$PWD
for ft in $(grep 'Y' data/selected.VEP.fields.txt | cut -d' ' -f 3); do 

# Consider splice donor,acceptor and region as single category: splice
if [[ $ft == 'splice_donor_variant' ]];then continue;fi
if [[ $ft == 'splice_acceptor_variant' ]];then continue;fi
if [[ $ft == 'splice_region_variant' ]];then ft='splice';fi 

echo $ft;

echo "
#PBS -N vep_to_TORUS
#PBS -S /bin/bash
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=4GB
#PBS -d $wdir
#PBS -o $wdir/logs/vep_to_TORUS.$ft.\$PBS_JOBID.out
#PBS -e $wdir/logs/vep_to_TORUS.$ft.\$PBS_JOBID.err

module load htslib;
cat <(echo \"SNP \"$ft\"_d\") <(awk 'FNR==NR{a[\$1]=1;next}{ if(a[\$1] > 0) {print \$1\" \"1} else {print \$1\" \"0}}' <(zcat /scratch/mpavia/TORUS/WGS_VEP_4torus.MAF01.chr*.gz | awk -v ft=$ft '{ if (\$7 ~ ft) {print}}') <(zcat data/gtex.WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.bed.gz | cut -f 4)) | bgzip > data/VEP.$ft.txt.gz;
" > logs/VEP.$ft.qsub
#qsub logs/VEP.$ft.qsub
ls logs/VEP.$ft.qsub

done
for ft in $(grep 'CY' data/selected.VEP.fields.txt | cut -d' ' -f 3); do

# Consider CDS subtypes as single category: CDS
if [[ $ft == 'start_retained_variant' || $ft == 'protein_altering_variant' || $ft == 'incomplete_terminal_codon_variant' || $ft == 'stop_retained_variant' || $ft == 'synonymous_variant' || $ft == 'stop_lost' || $ft == 'start_lost' || $ft == 'inframe_insertion' || $ft == 'stop_gained' || $ft == 'inframe_deletion' || $ft == 'frameshift_variant' || $ft == 'missense_variant' ]];then continue;fi

echo $ft;

echo "
#PBS -N vep_to_TORUS
#PBS -S /bin/bash
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=4GB
#PBS -d $wdir
#PBS -o $wdir/logs/vep_to_TORUS.$ft.\$PBS_JOBID.out
#PBS -e $wdir/logs/vep_to_TORUS.$ft.\$PBS_JOBID.err

module load htslib;
cat <(echo \"SNP \"$ft\"_d\") <(awk 'FNR==NR{a[\$1]=1;next}{ if(a[\$1] > 0) {print \$1\" \"1} else {print \$1\" \"0}}' <(zcat /scratch/mpavia/TORUS/WGS_VEP_4torus.MAF01.chr*.gz | awk -v ft=$ft '{ if (\$7 ~ ft || \$7 ~ \"start_retained_variant\" || \$7 ~ \"protein_altering_variant\" || \$7 ~ \"incomplete_terminal_codon_variant\" || \$7 ~ \"stop_retained_variant\" || \$7 ~ \"synonymous_variant\" || \$7 ~ \"stop_lost\" || \$7 ~ \"start_lost\" || \$7 ~ \"inframe_insertion\" || \$7 ~ \"stop_gained\" || \$7 ~ \"inframe_deletion\" || \$7 ~ \"frameshift_variant\" || \$7 ~ \"missense_variant\" ) {print}}') <(zcat data/gtex.WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.bed.gz | cut -f 4)) | bgzip > data/VEP.CDS.txt.gz;
" > logs/VEP.CDS.qsub
#qsub logs/VEP.$CDS.qsub
ls logs/VEP.CDS.qsub

done
