## Select data for mash
bash pipeline.eQTLs.sh

## Generate mashable files
R --vanilla < generate.mashable.file.eQTLs.R

## Run mash
## Comment the jobs sequentially, and allocate resources accordingly
bash mash.eQTLs.cmd

## Extract LFSR < 0.05 hits
R --vanilla < extract.signif.eQTLs.R
