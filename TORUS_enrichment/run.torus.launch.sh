# Beware of changing memory allocations for eQTL tissues. Def. 20GB, Muscle, Ovary, Prostate 50GB, Testis 70 GB
for tissue in $(cat ../../coloc/data/tissues.txt | sed 's/_//g'); do echo $tissue; bash scripts/run.torus.motifs.sh mqtls All ENSEMBL_REGULATORY_BUILD $tissue;done
for tissue in $(cat ../../coloc/data/tissues.txt | sed 's/_//g'); do echo $tissue; bash scripts/run.torus.motifs.sh eqtls All ENSEMBL_REGULATORY_BUILD $tissue;done

for tissue in $(cat ../../coloc/data/tissues.txt | sed 's/_//g'); do echo $tissue; bash scripts/run.torus.DNASeq.sh mqtls All DNASeq $tissue;done
for tissue in $(cat ../../coloc/data/tissues.txt | sed 's/_//g'); do echo $tissue; bash scripts/run.torus.DNASeq.sh eqtls All DNASeq $tissue;done

for block in $(seq 1 17); do
echo TF block $block
for tissue in $(cat ../../coloc/data/tissues.txt | sed 's/_//g'); do echo $tissue; bash scripts/run.torus.TF.sh mqtls TF $block $tissue;done
for tissue in $(cat ../../coloc/data/tissues.txt | sed 's/_//g' | grep -A 1000 Musc | head -n 3); do echo $tissue; bash scripts/run.torus.TF.sh eqtls TF $block $tissue;done
done
