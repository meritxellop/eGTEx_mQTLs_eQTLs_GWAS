cat <(for tissue in $(cat ../../coloc/data/tissues.txt | sed 's/_//g'); do grep -F -f <(cat ../../mashr/mQTLs/lfsr.005.txt ../../mashr/mQTLs/fdr.005.txt | grep $tissue | cut -f 1 | sort | uniq) ../../results/permuted/$tissue.regular.perm.fdr.txt | cut -f 2;done) | sort | uniq > data/mqtls/signif.top.mQTL.variants.txt