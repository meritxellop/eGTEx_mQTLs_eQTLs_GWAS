for tissue in $(cat ../../coloc/data/tissues.txt | sed 's/_//g'); do echo $tissue; bash scripts/generate_torus_inputfile.fastqtlgtex.lite.sh $tissue;done
