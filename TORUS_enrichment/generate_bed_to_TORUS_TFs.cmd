TFblock=$1
wdir=$PWD

get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

if [[ ! -f ../data/encRegTfbsClustered/TFs.shuf.txt ]]; then
echo Generate shuffled TF blocks of 20 TFs
shuf --random-source=<(get_seeded_random 42) ../data/encRegTfbsClustered/TFs.txt > ../data/encRegTfbsClustered/TFs.shuf.txt
for i in $(seq 1 17); do I=1;J=20;if [[ i > 1 ]]; then I=$(($(($i-1))*20+1)); J=$(($I+19));fi;awk 'NR==FNR{data[$1]; next}FNR in data' <(seq $I $J | tr ' ' '\n') ../data/encRegTfbsClustered/TFs.shuf.txt > ../data/encRegTfbsClustered/TFs.shuf.$i.txt;done
fi

echo "
#PBS -N generate_torus_format_TF
#PBS -S /bin/bash
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=2GB
#PBS -d $wdir
#PBS -o $wdir/logs/generate_torus_format_TF_"$block"_\$PBS_JOBID.out
#PBS -e $wdir/logs/generate_torus_format_TF_"$block"_\$PBS_JOBID.err

date
module load bedtools
module load htslib
for TFindex in \$(seq 1 20); do
	TF=\$(head -n \$TFindex ../data/encRegTfbsClustered/TFs.shuf.$TFblock.txt | tail -n1)
	echo $TFblock \$TFindex \$TF
	cat <(echo SNP \"\$TF\"_d) <(bedtools merge -i <(zgrep -w -F \$TF ../data/encRegTfbsClustered/original/encRegTfbsClustered.gz | cut -f 2,3,4,5 | grep -w -P '^chr(\d+|X)' |  sort -k1,1 -k2,2n)  -c 4 -o collapse |  bedtools intersect -loj -a data/gtex.WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.bed.gz -b stdin | awk '{ if(\$6<0) {print \$4\" \"0} else {print \$4\" \"1} }' ) | uniq > data/DNASeq.$TFblock.tmp.\$TFindex 
done
paste data/DNASeq.$TFblock.tmp.1 <(cut -d' ' -f 2 data/DNASeq.$TFblock.tmp.2) <(cut -d' ' -f 2 data/DNASeq.$TFblock.tmp.3) <(cut -d' ' -f 2 data/DNASeq.$TFblock.tmp.4) <(cut -d' ' -f 2 data/DNASeq.$TFblock.tmp.5) <(cut -d' ' -f 2 data/DNASeq.$TFblock.tmp.6) <(cut -d' ' -f 2 data/DNASeq.$TFblock.tmp.7) <(cut -d' ' -f 2 data/DNASeq.$TFblock.tmp.8) <(cut -d' ' -f 2 data/DNASeq.$TFblock.tmp.9) <(cut -d' ' -f 2 data/DNASeq.$TFblock.tmp.10) <(cut -d' ' -f 2 data/DNASeq.$TFblock.tmp.11) <(cut -d' ' -f 2 data/DNASeq.$TFblock.tmp.12) <(cut -d' ' -f 2 data/DNASeq.$TFblock.tmp.13) <(cut -d' ' -f 2 data/DNASeq.$TFblock.tmp.14) <(cut -d' ' -f 2 data/DNASeq.$TFblock.tmp.15) <(cut -d' ' -f 2 data/DNASeq.$TFblock.tmp.16) <(cut -d' ' -f 2 data/DNASeq.$TFblock.tmp.17) <(cut -d' ' -f 2 data/DNASeq.$TFblock.tmp.18) <(cut -d' ' -f 2 data/DNASeq.$TFblock.tmp.19) <(cut -d' ' -f 2 data/DNASeq.$TFblock.tmp.20) | tr ' ' '\t' | bgzip > torus_annots/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.TF.$TFblock.txt.gz 
rm data/DNASeq.$TFblock.tmp.*" > logs/DNASeq.$TFblock.qsub
qsub logs/DNASeq.$TFblock.qsub
