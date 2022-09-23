wdir=/gpfs/data/gtex-group/mpavia/methylation

# Generate $wdir/results/permuted/signif.top_eQTL.fdr005.anytissue.txt. Contains basic QTL stats for any eQTL hit in any tissue (top variant per gene). Standard - not conditional - QTL analysis 
# Select top SNP per eGene
awk ' NR==1 { gene=$2;genepval[gene]=$7;line=$2":"$3 } { if ($2 == gene) { if ($7 < genepval[gene]) { genepval[gene]=$7; line=$2":"$3 } } else { print line; gene=$2; genepval[gene]=$7; line=$2":"$3 } } END { print line }' <(sort -k2,2 $wdir/results/permuted/signif.top_eQTL.fdr005.anytissue.all.txt) > $wdir/results/permuted/signif.top_eQTL.fdr005.anytissue.all.top.txt

# Gather QTL stats across tissues, for $wdir/results/permuted/signif.top_eQTL.fdr005.anytissue.all.top.txt
for tissue in $(cat ../data/tissues.txt); do sed "s/\$tissue/$tissue/g" scripts/select.input.for.mash.eQTLs.all.sh > commands/tmp.$tissue.eQTLs.all.cmd; qsub commands/tmp.$tissue.eQTLs.all.cmd; done

for tissue in $(cat ../data/tissues.txt | sed 's/_//g'); do
	rm $tissue.eQTL_top.fdr005.txt;
	echo $tissue;
	for id in $(cat ../data/fileids.txt); do 
		cat $tissue.eQTL_top.fdr005.$id.txt >> $tissue.eQTL_top.fdr005.txt; 
		rm $tissue.eQTL_top.fdr005.$id.txt;
	done;
done

bash scripts/generate.mash.input.eQTLs.sh

mkdir eQTLs
mv *.eQTL_top.fdr005.all.txt eQTLs


for id in $(cat ../data/fileids.txt); do sed "s/\$id/$id/g" scripts/select.input.for.mash.eQTLs.all.random.sh > commands/tmp.$id.eQTLs.all.random.cmd; qsub commands/tmp.$id.eQTLs.all.random.cmd; done

for tissue in $(cat ../data/tissues.txt | sed 's/_//g'); do
        rm $tissue.eQTL_top.fdr005.random.txt;
        echo $tissue;
        for id in $(cat ../data/fileids.txt); do
                cat $tissue.eQTL_top.fdr005.$id.random.txt >> $tissue.eQTL_top.fdr005.random.txt;
                rm $tissue.eQTL_top.fdr005.$id.random.txt;
        done;
done

bash scripts/generate.mash.input.eQTLs.random.sh

mv *.eQTL_top.fdr005.random.all.txt eQTLs

