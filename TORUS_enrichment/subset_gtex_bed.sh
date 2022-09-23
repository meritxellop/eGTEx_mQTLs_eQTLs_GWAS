module load htslib

# Subset variants bed file to those used in the GTEx main paper
# https://science.sciencemag.org/content/369/6509/1318

date
echo 'Subseting gtex.bed.gz'
awk 'FNR==NR{a[$4]=$0;next}{print a[$1]}' <(zcat /gpfs/data/gtex-group/sex_biased_regulation_v8/TORUS/data/gtex.bed.gz) <(zcat torus_annots/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.txt.gz | tail -n+2 | cut -f 1) | bgzip > /gpfs/data/gtex-group/sex_biased_regulation_v8/TORUS/data/gtex.WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.bed.gz
echo 'Done'
date
