wdir=$PWD;
tissue=$1 # e.g. Lung
qtl=$2 # e.g. mqtls/eqtls
module load htslib

if [[ $qtl == 'mqtls' ]]; then
	tissue=$(echo $tissue | sed 's/_//g')
	egene_allSNPGenepairs_file=/gpfs/data/pierce-lab/GTEx/mQTLs/mQTL_mapping/regular/$tissue.mQTLs.regular.txt.gz
	egene_allSNPGenepairs_file_torus=/scratch/mpavia/TORUS/$tissue.allpairs_4torus.$qtl.txt.gz

	echo "Convert mqtl fastqtlfile-gtex to fastqtlfile-regular "

	zcat $egene_allSNPGenepairs_file | tail -n +2 - | awk '{print $1,$2,$3,$7,$8,$9}' - | bgzip > $egene_allSNPGenepairs_file_torus
	ln -s $egene_allSNPGenepairs_file_torus $wdir/data/mqtls/$tissue.allpairs_4torus.txt.gz
elif [[ $qtl == 'eqtls' ]]; then
	egene_allSNPGenepairs_file=/gpfs/data/gtex-group/v8/59349/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_all_associations/$tissue.allpairs.txt.gz
	tissue=$(echo $tissue | sed 's/_//g')
	egene_allSNPGenepairs_file_torus=/scratch/mpavia/TORUS/$tissue.allpairs_4torus.$qtl.txt.gz

	echo "Convert eqtl fastqtlfile-gtex to fastqtlfile-regular "

	zcat $egene_allSNPGenepairs_file | tail -n +2 - | awk '{print $1,$2,$3,$7,$8,$9}' - | bgzip > $egene_allSNPGenepairs_file_torus
	ln -s $egene_allSNPGenepairs_file_torus $wdir/data/eqtls/$tissue.allpairs_4torus.txt.gz
fi
