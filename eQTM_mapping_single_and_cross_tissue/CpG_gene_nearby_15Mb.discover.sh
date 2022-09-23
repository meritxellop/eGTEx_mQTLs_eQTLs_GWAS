module load bedtools
for num in $(cat ../data/support_files/fileids.txt); do
echo $num
	bedtools intersect -a <(grep -F -m 500 -f ../data/support_files/lists/cpgs_by_fileid/$num.cpgs.txt /gpfs/data/gtex-group/mpavia/methylation/coloc/data/MethylationEPIC_v-1-0_B4.filtered.hg38.CpGblocks.bed | awk '{ if( $2-1500000 < 0 ) { print $1"\t"0"\t"$3+1500000"\t"$4"\t"$5 } else { print $1"\t"$2-1500000"\t"$3+1500000"\t"$4"\t"$5 }}') -b <(awk '{ if($5 == "+") { print $1"\t"$2"\t"$2+1"\t"$4"\t"$5 } else { print $1"\t"$3"\t"$3+1"\t"$4"\t"$5  }}' /gpfs/data/gtex-group/mpavia/methylation/coloc/data/gencode.v26.GRCh38.genes.strand.bed) -wa -wb | awk '{ if($10 == "+") { print $4"\t"$5"\t"$9"\t"($3-1500000-1)-$7 } else { print $4"\t"$5"\t"$9"\t"$7-($3-1500000-1) } }' > CpG_gene_overlaps/CpG_gene_1500kb.$num.txt 
done
