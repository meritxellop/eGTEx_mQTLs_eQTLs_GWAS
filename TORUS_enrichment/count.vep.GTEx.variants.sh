#zgrep -v '^#' /gpfs/data/gtex-group/sex_biased_regulation_v8/data/Genotype_Data/raw/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz | wc -l
cat <(for chr in $(echo X $(seq 1 22)); do zcat /scratch/mpavia/TORUS/WGS_VEP_4torus.MAF01.chr$chr.gz;done) | cut -f 7 | grep -v '^#' | grep -v Consequence | sed 's/,/\n/g' | sort | uniq -c | sort -gk1
