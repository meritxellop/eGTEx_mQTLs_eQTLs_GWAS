wdir=$PWD;
tissue=$1 # e.g. Lung
qtl=mqtls
module load htslib

tissue=$(echo $tissue | sed 's/_//g')
egene_allSNPGenepairs_file_torus=/scratch/mpavia/TORUS/$tissue.allpairs_4torus.$qtl.lite.txt

echo "Convert mqtl fastqtlfile-gtex to fastqtlfile-regular "

if [[ ! -f $wdir/data/fileids.21k.random_subset.txt ]]; then
get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}
sort --random-source=<(get_seeded_random 42) -R ../../coloc/data/fileids.txt  | head -42 > $wdir/data/fileids.21k.random_subset.txt
fi

> $egene_allSNPGenepairs_file_torus
for block in $(cat $wdir/data/fileids.21k.random_subset.txt); do
	echo $block
	egene_allSNPGenepairs_file=/scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$block.regular.txt.gz
	zcat $egene_allSNPGenepairs_file | tail -n +2 - | awk '{print $1,$2,$3,$7,$8,$9}' >> $egene_allSNPGenepairs_file_torus
done
bgzip -f $egene_allSNPGenepairs_file_torus
ln -s $egene_allSNPGenepairs_file_torus.gz $wdir/data/mqtls/$tissue.allpairs_4torus.lite.txt.gz
