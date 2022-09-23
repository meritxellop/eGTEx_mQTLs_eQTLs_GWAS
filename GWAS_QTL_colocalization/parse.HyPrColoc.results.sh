echo "Count GWAShit-eQTL instances"
mQTL_specific_eQTL_maps=$(grep -v '^eQTL' results/mQTL_specific.summary.$thres.txt | cut -f 13 | wc -l);
echo $mQTL_specific_eQTL_maps
mQTL_shared_eQTL_maps=$(grep -v '^eQTL' results/mQTL_shared.summary.$thres.txt | cut -f 13 | wc -l);
echo $mQTL_shared_eQTL_maps

echo "Count GWAS hit instances"
mQTL_specific_eQTL_maps=$(grep -v '^eQTL' results/mQTL_specific.summary.$thres.txt | cut -f 2,3 | sort | uniq |  wc -l );
echo $mQTL_specific_eQTL_maps "/" $(tail -n+2 data/mQTL.specific.colocs.GWAS.hits.txt | wc -l)
mQTL_shared_eQTL_maps=$(grep -v '^eQTL'  results/mQTL_shared.summary.$thres.txt | cut -f 2,3 | sort | uniq |  wc -l );
echo $mQTL_shared_eQTL_maps "/" $(tail -n+2 data/mQTL.shared.colocs.GWAS.hits.txt | wc -l)

# Count how many instances of each eQTL dataset exist for mQTL-shared and mQTL-specific colocalizations
join -1 1 -2 2 -a 1 <(join -1 2 -2 2 <(grep -v '^eQTL' results/mQTL_specific.summary.$thres.txt | cut -f 13 |  sort | uniq -c | sed 's/^ *//g') <(tail -n+2 data/external_eQTLs.trait.classes.2.txt | sort -k2)) <(grep -v '^eQTL'  results/mQTL_shared.summary.$thres.txt | cut -f 13 | sort | uniq -c | sed 's/^ *//g') | awk '{print $3,$1,$2,$4}' > results/eQTL.source.counts.$thres.txt
