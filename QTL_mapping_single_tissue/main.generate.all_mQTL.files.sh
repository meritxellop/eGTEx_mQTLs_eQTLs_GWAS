# Generate bed files of inverse rank-normalized methylation values
for tissue in $(cat data/support_files/tissues.txt); do
	bash scripts/generate.methylation.bedFile.sh $tissue
done

# Split bedFiles of methylation, generate CpG lists per block files
bash scripts/split.beds.sh
#tail -n+2 $file | cut -f 4 > /gpfs/data/gtex-group/mpavia/methylation/data/support_files/lists/cpgs_by_fileid/$TISSUE.$fileid.cpgs.all_mQTLs.txt

# Generate primary signal mQTL block files
for tissue in $(cat data/support_files/tissues.txt); do if [[ $tissue =~ .*_.* ]] ; then tissue=$(echo $tissue | sed 's/_//g'); fi; for id in $(cat coloc/data/fileids.txt); do echo $id; bash scripts/nominal.mQTL.fastqtl.peers.split.cmd $tissue $id regular;done;done

# Generate conditional signal mQTL block files
for tissue in $(cat data/support_files/tissues.txt); do if [[ $tissue =~ .*_.* ]] ; then tissue=$(echo $tissue | sed 's/_//g'); fi; for id in $(cat coloc/data/fileids.txt); do echo $id; bash scripts/nominal.conditional.mQTL.fastqtl.peers.split.calls.all_mQTLs.cmd $tissue $id;done;done

for tissue in $(cat data/support_files/tissues.txt); do if [[ $tissue =~ .*_.* ]] ; then tissue=$(echo $tissue | sed 's/_//g'); fi; > /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.bwd.all_mQTLs.ALL_cpgs.txt; for id in $(cat coloc/data/fileids.txt); do echo $id; cat /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$id.bwd.all_mQTLs.ALL_cpgs.txt >> /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.bwd.all_mQTLs.ALL_cpgs.txt; done;done
