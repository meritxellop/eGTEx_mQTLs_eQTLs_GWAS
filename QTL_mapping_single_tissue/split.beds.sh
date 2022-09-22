module load htslib;
for tissue in $(cat data/support_files/tissues.txt); do
	if [[ $tissue =~ .*_.* ]] ; then tissue=$(echo $tissue | sed 's/_//g'); fi
	echo $tissue;
	#mkdir -p /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles
	#rm /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/"$tissue"_[0-9][0-9][0-9][0-9].gz*
	header=$(zcat data/Phenotypes/bedfiles/$tissue.bed.gz | head -n 1) 
	zcat data/Phenotypes/bedfiles/$tissue.bed.gz | tail -n+2 | split - /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/"$tissue"_ -l 500 -d -a 4
	#for fileid in $(ls /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/"$tissue"_* | cut -d'_' -f 2); do
	for fileid in $( cat coloc/data/fileids.txt ); do 
		if [[ ! -f /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/"$tissue"_"$fileid".gz ]]; then
		echo $fileid
		file=/scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/"$tissue"_"$fileid"
		cmd=$(echo 'sed -i "1 i\'$header'"' $file);
		eval "$cmd"
		sed -i 's/ /\t/g' $file
		bgzip -f $file;
		tabix -f -p bed $file.gz;
	fi
	done
done
