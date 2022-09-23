module load bedtools; 
module load htslib;

# Generate TORUS-format CpG Island and CpG Island shore anotation
# cpgIslandExt.gz is extracted from UCSC (Jan 2021)
# CpG island shores are defined as 2-kb-long regions that lie on both sides of a CpG island [16]. CpG islands are, on average, 1000 base pairs (bp).
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0239196

echo 'Generate TORUS-format CpG Island annotation'
cat <(echo 'SNP cpgIsland_d') <(zcat ../data/cpgIslandExt/converted/cpgIslandExt.gz | bedtools intersect -loj -a data/gtex.WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.bed.gz -b stdin | awk '{ if($6<0) {print $4" "0} else {print $4" "1} }' ) | bgzip > data/cpgIslandExt.txt.gz

echo 'Generate TORUS-format CpG Island Shore annotation (-/+ 2Kb from CpG Island boundaries)'

cat <(echo 'SNP cpgIslandShore_d') <(zcat ../data/cpgIslandExt/converted/cpgShoreExt.gz | bedtools intersect -loj -a data/gtex.WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.bed.gz -b stdin | awk '{ if($6<0) {print $4" "0} else {print $4" "1} }' ) | bgzip > data/cpgIslandShoreExt.txt.gz

echo 'Generate TORUS-format CpG Island Shelf annotation (-/+ 2Kb from CpG Island Shore boundaries)'

cat <(echo 'SNP cpgIslandShelf_d') <(zcat ../data/cpgIslandExt/converted/cpgShelfExt.gz | bedtools intersect -loj -a data/gtex.WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.bed.gz -b stdin | awk '{ if($6<0) {print $4" "0} else {print $4" "1} }' ) | bgzip > data/cpgIslandShelfExt.txt.gz

# Generate TORUS-format ENCODE Ccre anotation
# encodeCcreCombined.gz is extracted from UCSC (Jan 2021). Matches ENCODE5 release.

echo 'Generate TORUS-format ENCODE Ccre annotation'
echo 'CTCF'

cat <(echo 'SNP CTCF_d') <(zcat ../data/encodeCcreCombined/converted/encodeCcreCombined.CTCF.bed.gz | bedtools intersect -loj -a data/gtex.WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.bed.gz -b stdin | awk '{ if($6<0) {print $4" "0} else {print $4" "1} }' ) | uniq | bgzip > data/encodeCcreCombined.CTCF.txt.gz

echo 'enhD'
cat <(echo 'SNP enhD_d') <(zcat ../data/encodeCcreCombined/converted/encodeCcreCombined.enhD.bed.gz | bedtools intersect -loj -a data/gtex.WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.bed.gz -b stdin | awk '{ if($6<0) {print $4" "0} else {print $4" "1} }' ) | uniq | bgzip > data/encodeCcreCombined.enhD.txt.gz

echo 'enhP'
cat <(echo 'SNP enhP_d') <(zcat ../data/encodeCcreCombined/converted/encodeCcreCombined.enhP.bed.gz | bedtools intersect -loj -a data/gtex.WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.bed.gz -b stdin | awk '{ if($6<0) {print $4" "0} else {print $4" "1} }' ) | uniq | bgzip > data/encodeCcreCombined.enhP.txt.gz

echo 'prom'
cat <(echo 'SNP prom_d') <(zcat ../data/encodeCcreCombined/converted/encodeCcreCombined.prom.bed.gz | bedtools intersect -loj -a data/gtex.WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.bed.gz -b stdin | awk '{ if($6<0) {print $4" "0} else {print $4" "1} }' ) | uniq | bgzip > data/encodeCcreCombined.prom.txt.gz

# Generate TORUS-format DNASe-seq anotation
# narrow Peaks at FDR=5% are extracted from ENCODE official website (Jan 2021), ENCODE5 release

cat <(echo 'SNP DNASeqBreast_d') <(zcat ../data/ENCODE5/converted/ENCODE5.Breast.DNASeq.bed.gz | bedtools intersect -loj -a data/gtex.WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.bed.gz -b stdin | awk '{ if($6<0) {print $4" "0} else {print $4" "1} }' ) | uniq | bgzip > data/DNASeq.BreastMammaryTissue.txt.gz
cat <(echo 'SNP DNASeqColon_d') <(zcat ../data/ENCODE5/converted/ENCODE5.Colon.DNASeq.bed.gz | bedtools intersect -loj -a data/gtex.WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.bed.gz -b stdin | awk '{ if($6<0) {print $4" "0} else {print $4" "1} }' ) | uniq | bgzip > data/DNASeq.ColonTransverse.txt.gz
cat <(echo 'SNP DNASeqKidney_d') <(zcat ../data/ENCODE5/converted/ENCODE5.Kidney.DNASeq.bed.gz | bedtools intersect -loj -a data/gtex.WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.bed.gz -b stdin | awk '{ if($6<0) {print $4" "0} else {print $4" "1} }' ) | uniq | bgzip > data/DNASeq.KidneyCortex.txt.gz
cat <(echo 'SNP DNASeqLung_d') <(zcat ../data/ENCODE5/converted/ENCODE5.Lung.DNASeq.bed.gz | bedtools intersect -loj -a data/gtex.WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.bed.gz -b stdin | awk '{ if($6<0) {print $4" "0} else {print $4" "1} }' ) | uniq | bgzip > data/DNASeq.Lung.txt.gz
cat <(echo 'SNP DNASeqMuscle_d') <(zcat ../data/ENCODE5/converted/ENCODE5.Muscle.DNASeq.bed.gz | bedtools intersect -loj -a data/gtex.WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.bed.gz -b stdin | awk '{ if($6<0) {print $4" "0} else {print $4" "1} }' ) | uniq | bgzip > data/DNASeq.MuscleSkeletal.txt.gz
cat <(echo 'SNP DNASeqOvary_d') <(zcat ../data/ENCODE5/converted/ENCODE5.Ovary.DNASeq.bed.gz | bedtools intersect -loj -a data/gtex.WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.bed.gz -b stdin | awk '{ if($6<0) {print $4" "0} else {print $4" "1} }' ) | uniq | bgzip > data/DNASeq.Ovary.txt.gz
cat <(echo 'SNP DNASeqProstate_d') <(zcat ../data/ENCODE5/converted/ENCODE5.Prostate.DNASeq.bed.gz | bedtools intersect -loj -a data/gtex.WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.bed.gz -b stdin | awk '{ if($6<0) {print $4" "0} else {print $4" "1} }' ) | uniq | bgzip > data/DNASeq.Prostate.txt.gz
cat <(echo 'SNP DNASeqTestis_d') <(zcat ../data/ENCODE5/converted/ENCODE5.Testis.DNASeq.bed.gz | bedtools intersect -loj -a data/gtex.WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.bed.gz -b stdin | awk '{ if($6<0) {print $4" "0} else {print $4" "1} }' ) | uniq | bgzip > data/DNASeq.Testis.txt.gz

# Generate TORUS-format ENCODE TF binding cluster anotation
# encRegTfbsClustered.gz  is extracted from UCSC (Jan 2021).

for i in $(seq 1 17); do bash scripts/generate_bed_to_TORUS_TFs.cmd $i;done

echo 'Merge gene regulatory cCRE with VEP annotations"
paste <(zcat data/cpgIslandExt.txt.gz) <(zcat data/cpgIslandShoreExt.txt.gz | cut -d' ' -f 2) <(zcat data/cpgIslandShelfExt.txt.gz | cut -d' ' -f 2) <(zcat data/encodeCcreCombined.CTCF.txt.gz | cut -d' ' -f 2 ) <(zcat data/encodeCcreCombined.enhD.txt.gz | cut -d' ' -f 2 ) <(zcat data/encodeCcreCombined.enhP.txt.gz | cut -d' ' -f 2 ) <(zcat data/encodeCcreCombined.prom.txt.gz | cut -d' ' -f 2 ) <(zcat data/VEP.splice.txt.gz | cut -d' ' -f 2 )  <(zcat data/VEP.5_prime_UTR_variant.txt.gz | cut -d' ' -f 2 ) <(zcat data/VEP.3_prime_UTR_variant.txt.gz | cut -d' ' -f 2) <(zcat data/VEP.intron_variant.txt.gz | cut -d' ' -f 2) <(zcat data/VEP.non_coding_transcript_exon_variant.txt.gz | cut -d' ' -f 2) <(zcat data/VEP.CDS.txt.gz | cut -d' ' -f 2) | tr ' ' '\t' | bgzip > torus_annots/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.complemented.txt.gz

echo 'Merge DNASeq annotations"
paste <(zcat data/DNASeq.BreastMammaryTissue.txt.gz) <(zcat data/DNASeq.ColonTransverse.txt.gz | cut -d' ' -f 2) <(zcat data/DNASeq.KidneyCortex.txt.gz | cut -d' ' -f 2) <(zcat data/DNASeq.Lung.txt.gz | cut -d' ' -f 2) <(zcat data/DNASeq.MuscleSkeletal.txt.gz | cut -d' ' -f 2) <(zcat data/DNASeq.Ovary.txt.gz | cut -d' ' -f 2) <(zcat data/DNASeq.Prostate.txt.gz | cut -d' ' -f 2) <(zcat data/DNASeq.Testis.txt.gz | cut -d' ' -f 2) | tr ' ' '\t' | bgzip > torus_annots/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.DNASeq.txt.gz


#for tissue in $(echo "Breast Colon Kidney Lung Muscle Ovary Prostate Testis"); do
        echo $tissue;
        bedtools merge -i <(zcat data/ENCODE5/original/$tissue.1.DNASeq.bed.gz | awk -v tissue=$tissue '{print $1"\t"$2"\t"$3"\t"tissue }' | sort -k1,1 -k2,2n) -c 4 -o collapse | bgzip > data/ENCODE5/converted/ENCODE5.$tissue.DNASeq.bed.gz
done
