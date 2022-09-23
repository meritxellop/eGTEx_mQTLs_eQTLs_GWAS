#### Perform an eQTM across-tissue multivariate analysis with MASHR, to leverage shared signal across tissues and increase the eQTM discovery set

# Generate PEER-residualized expression and methylation values
bash scripts/generate.residualized.methylation.and.expression.sh 

# Test for eQTMs by spearman correlation on PEER-residualized expression and methylation values. Exhaustive testing in trans: in each tissue, for each CpG, all expressed genes are tested.
bash scripts/run.CpG_gene_correlation.launch.sh

# Calculate overlaps of CpG cis (-/+ 1.5Mb) regions with genes, annotate CpG-TSS distance
bash scripts/CpG_gene_nearby_15Mb.discover.sh

# Extract CpG-gene cis eQTM tests, calculate significance after adjusting for multiple testing at within and across CpG level, extract ststs corresponding to eQTMs significant (qvalue<0.05) in at leas one tissue
bash scripts/run.extract.cor.stats.per.block.sh
bash scripts/joincor.stats.per.block.sh
bash scripts/calculate.qval.sh 
bash scripts/signif.cpg.gene.pairs.qvalue.sh
bash scripts/extract.acrosstissues.for_mash.signif.cpg.gene.pairs.stats.sh

# Select N=100k random CpG-gene correlations, from cis (-/+ 1.5Mb) CpG-gene correlation tests. 
# Importantly, the distribution of the corresponding p-values is inflated, indicating, that there is non-null CpG-gene correlation in CpG cis (-/+ 1.5Mb) regions. 
# If we select this background as oposed to selecting CpG-gene correlations throughout the genome, we are being more strict.
# This option has been selected to mimic MASHR eQTL analysis (GTEx main paper), where cis eQTL tests were selected. 
bash scripts/extract.acrosstissues.for_mash.random.cpg.gene.pairs.stats.sh

# Alternatively, one can select CpG-gene correlations throughout the genome
#bash scripts/run.CpG_gene_random.discover.sh
#bash scripts/run.extract.cor.stats.per.block.random.sh
#bash scripts/joincor.stats.per.block.random.sh
#bash scripts/extract.acrosstissues.for_mash.random.trans.cpg.gene.pairs.stats.sh

# Generate MASHR input. Spearman correlation values (rhos) are z-transformed
bash scripts/generate.mash.input.strong.sh
bash scripts/generate.mash.input.random.sh

R --vanilla < scripts/generate.mashable.file.Fishers_z_transformation.R

# Run MASHR. Perform each step separately, by commenting the other ones subsequently. Follow indications to alocate time and memory
bash scripts/mash.CpG.gene.cors.cmd

# Wrap up stats for significant cases, defined as the union of a) eQTMs significant ( bonferroni-corrected P<0.05 + qvalue<0.05 ) or b) LFSR<0.05 in at least one tissue
bash scripts/run.CpG_gene_correlation.launch.sh
