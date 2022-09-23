# Select random 21k mQTL assocs to run TORUS. For eQTLs, the whole summary stats will be used.
for tissue in $(cat ../../coloc/data/tissues.txt | sed 's/_//g'); do echo $tissue; bash scripts/generate_torus_inputfile.fastqtlgtex.lite.sh $tissue;done

## Run TORUS for different annotations across tissues and QTL types 
## For that, one needs to have the annotations in place. Some scripts are provided here to generate annotations
bash run.torus.launch.sh

# Use summarize*sh scripts to generate the input for the meta-analysis. Additional scripts to extract the overlap of m/QTLs with annotations (summarize*R) are provided.
bash summarize.results.Chmm.active_states.1.sh
R --vanilla < summarize.results.Chmm.active_states.2.R
bash summarize.results.Chmm.active_states.sh
bash summarize.results.DNASeq.1.sh
bash summarize.results.TF.1.sh
R --vanilla < summarize.results.TF.2.R
bash summarize.results.vep.1.sh
R --vanilla < summarize.results.vep.2.R


## Perform cross-tissue meta-analysis
R --vanilla < scripts/perform.metaanalysis.random_effects.model.Chmm.active_states.R
R --vanilla < scripts/perform.metaanalysis.random_effects.model.TF.R
R --vanilla < scripts/perform.metaanalysis.random_effects.model.vep.R
