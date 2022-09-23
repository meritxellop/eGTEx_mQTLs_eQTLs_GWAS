tissue=$1
wdir=$PWD

get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

echo "
#PBS -N generate_torus_format_tissue
#PBS -S /bin/bash
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=2GB
#PBS -d $wdir
#PBS -o $wdir/logs/generate_torus_format_tissue_"$block"_\$PBS_JOBID.out
#PBS -e $wdir/logs/generate_torus_format_tissue_"$block"_\$PBS_JOBID.err

date
module load bedtools
module load htslib
for index in \$(seq 1 18); do
	Chmm=\$(grep -w \"^\$index\"  ../data/ROADMAP/README | cut -f 2)
	echo \$Chmm \$tissue
	cat <(echo SNP \"\$Chmm\"_d) <(bedtools merge -i <(zgrep -F -w \$index'_'\$Chmm ../data/ROADMAP/original/Lung.ChmmModels.bed.gz | grep -w -P '^chr(\d+|X)' |  sort -k1,1 -k2,2n)  -c 4 -o collapse |  bedtools intersect -loj -a data/gtex.WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.bed.gz -b stdin | awk '{ if(\$6<0) {print \$4\" \"0} else {print \$4\" \"1} }' ) | uniq > data/Chmm.$tissue.tmp.\$index 
done
#paste data/Chmm.$tissue.tmp.1 <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.2) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.3) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.4) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.5) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.6) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.7) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.8) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.9) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.10) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.11) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.12) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.13) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.14) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.15) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.16) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.17) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.18) | tr ' ' '\t' | bgzip > torus_annots/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.Chmm.$tissue.txt.gz 
paste data/Chmm.$tissue.tmp.1 <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.2) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.3) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.4) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.5) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.6) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.7) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.8) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.9) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.10) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.11) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.12) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.14) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.15) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.16) <(cut -d' ' -f 2 data/Chmm.$tissue.tmp.17) | tr ' ' '\t' | bgzip > torus_annots/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.Chmm.$tissue.txt.gz 
rm data/Chmm.$tissue.tmp.*" > logs/Chmm.$tissue.qsub
qsub logs/Chmm.$tissue.qsub
