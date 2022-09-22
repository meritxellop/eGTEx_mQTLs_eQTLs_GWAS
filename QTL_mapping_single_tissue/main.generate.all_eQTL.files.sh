# Split bedFiles of expression, generate gene lists per block files
bash scripts/split.all_eQTLs.beds.sh

# Generate primary signal eQTL block files
for tissue in $(cat data/support_files/tissues.txt); do TISSUE=$tissue ; if [[ $tissue =~ .*_.* ]] ; then tissue=$(echo $tissue | sed 's/_//g'); fi; for id in $(grep $tissue coloc/data/fileids.all_eQTLs.txt | cut -f 2); do echo $id; bash scripts/generate.eQTL.primary_signal.blocks.all_eQTLs.cmd $TISSUE $id;done;done

# Generate conditional signal eQTL block files
for tissue in $(cat data/support_files/tissues.txt); do TISSUE=$tissue ; if [[ $tissue =~ .*_.* ]] ; then tissue=$(echo $tissue | sed 's/_//g'); fi; for id in $(grep $tissue coloc/data/fileids.all_eQTLs.txt | cut -f 2); do echo $id; bash scripts/nominal.conditional.mQTL.fastqtl.peers.split.calls.all_eQTLs.cmd $TISSUE $id;done;done

for tissue in $(cat data/support_files/tissues.txt); do TISSUE=$tissue ; if [[ $tissue =~ .*_.* ]] ; then tissue=$(echo $tissue | sed 's/_//g'); fi; > /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.bwd.all_eQTLs.ALL_genes.txt; for id in $(grep $tissue coloc/data/fileids.all_eQTLs.txt | cut -f 2); do echo $id; cat /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.$id.bwd.all_eQTLs.ALL_genes.txt >> /scratch/mpavia/fastQTL/results/$tissue/mqtls/NonPermutedResults/tmpFiles/$tissue.bwd.all_eQTLs.ALL_genes.txt; done;done
