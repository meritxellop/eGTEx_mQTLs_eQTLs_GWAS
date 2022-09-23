for tissue in $(grep -v Kidney ../coloc/data/tissues.txt); do TISSUE=$(echo $tissue | sed 's/_//g'); for num in $(cat ../data/support_files/fileids.txt); do echo $tissue $num; cat scripts/extract.cor.stats.per.block.cmd | sed "s/num/$num/g" | sed "s/tissue/$TISSUE/g" > logs/extract.cor.stats.per.block.$num.$tissue.cmd;if [[ ! -s /scratch/mpavia/fastQTL/results/$TISSUE/mqtls/NonPermutedResults/tmpFiles/$TISSUE.$num.cor.stats.txt ]]; then echo TODO $tissue $num; bash logs/extract.cor.stats.per.block.$num.$tissue.cmd; ls logs/extract.cor.stats.per.block.$num.$tissue.cmd;fi ;done; done
