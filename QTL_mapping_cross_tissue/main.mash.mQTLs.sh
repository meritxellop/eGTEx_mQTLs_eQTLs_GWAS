## Select data for mash
bash pipeline.mQTLs.sh

## Generate mashable files
R --vanilla < generate.mashable.file.mQTLs.R

## Run mash
## Comment the jobs sequentially, and allocate resources accordingly
bash mash.mQTLs.cmd

## Extract LFSR < 0.05 hits
R --vanilla < extract.signif.eQTLs.R
