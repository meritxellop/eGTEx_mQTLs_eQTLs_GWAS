wdir=/gpfs/data/gtex-group/mpavia/methylation
awk '{ if(sqrt($5^2) > effsize[$2]) { effsize[$2]=$5; snp[$2]=$3;} } END { for (key in effsize) { print key":"snp[key] } }' $wdir/results/permuted/signif.top_mQTL.fdr005.anytissue.txt > $wdir/results/permuted/signif.top_mQTL.fdr005.anytissue.top.txt

# Select top SNP per mCpG
awk ' NR==1 { gene=$2;genepval[gene]=$7;line=$2":"$3 } { if ($2 == gene) { if ($7 < genepval[gene]) { genepval[gene]=$7; line=$2":"$3  } } else { print line; gene=$2; genepval[gene]=$7; line=$2":"$3  } } END { print line }' <(sort -k2,2 $wdir/results/permuted/signif.top_mQTL.fdr005.anytissue.txt) > $wdir/results/permuted/signif.top_mQTL.fdr005.anytissue.top.txt

# Generate ../results/permuted/signif.top_mQTL.fdr005.anytissue.txt
# Any mQTL hit in any tissue (top variant per gene). Standard - not conditional - QTL analysis

for id in $(cat ../data/fileids.txt); do sed "s/\$id/$id/g" scripts/select.input.for.mash.sh > commands/tmp.$id.cmd; qsub commands/tmp.$id.cmd; done

for tissue in $(cat ../data/tissues.txt | sed 's/_//g'); do
	rm $tissue.mQTL_top.fdr005.txt;
	echo $tissue;
	for id in $(cat ../data/fileids.txt); do 
		cat $tissue.mQTL_top.fdr005.$id.txt >> $tissue.mQTL_top.fdr005.txt; 
		rm $tissue.mQTL_top.fdr005.$id.txt;
	done;
done

bash scripts/generate.mash.input.sh

mkdir mQTLs
mv *.mQTL_top.fdr005.txt mQTLs

for id in $(cat ../data/fileids.txt); do sed "s/\$id/$id/g" scripts/select.input.for.mash.random.sh > commands/tmp.$id.random.cmd; qsub commands/tmp.$id.random.cmd; done

for tissue in $(cat ../data/tissues.txt | sed 's/_//g'); do
        rm $tissue.mQTL_top.fdr005.random.txt;
        echo $tissue;
        for id in $(cat ../data/fileids.txt); do
                cat $tissue.mQTL_top.fdr005.$id.random.txt >> $tissue.mQTL_top.fdr005.random.txt;
                rm $tissue.mQTL_top.fdr005.$id.random.txt;
        done;
done

bash scripts/generate.mash.input.random.sh

mv *.mQTL_top.fdr005.random.txt mQTLs
