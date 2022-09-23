#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
args = commandArgs(trailingOnly=TRUE)
tissue=args[1]

library(psych,lib.loc='/gpfs/apps/haswell/software/gcc-6.2.0/R/3.5.0/lib64/R/library')
library(reshape2)
library(readr)

# Load significant eQTMs univariate
fdr005=read.table('results/all.signif.cpg.gene.pairs.qvalue.005.txt',header=T)
fdr005$eQTMtissue=paste0(fdr005$cpg,':',fdr005$gene,':',fdr005$Tissue)
fdr005=subset(fdr005,Tissue%in%tissue)

# Load mash results
load(file='mash_output/CpG.gene.cors.mash.output.Robj')

# Load mash input
dat = readRDS('results/CpG.gene.cors.mash.rds')

# Fill with NAs tests with no input
CpG.gene.cors.m2$result$lfsr[is.na(dat$strong.b)]<-NA

# Format results
# Extract MASHR stats for eQTMs univariate (FDR<0.05) + multivariate (LFSR<0.05) hits
stats=na.omit(cbind(melt(CpG.gene.cors.m2$result$lfsr),fisherz2r(melt(CpG.gene.cors.m2$result$PosteriorMean)[,3]),melt(CpG.gene.cors.m2$result$PosteriorMean)[,3],melt(CpG.gene.cors.m2$result$PosteriorSD)[,3]))
colnames(stats)=c('eQTM','Tissue','lfsr','smoothed_rho','PosteriorMean','PosteriorSD')
stats$eQTMtissue=paste0(stats$eQTM,':',stats$Tissue)
stats=subset(stats,Tissue%in%tissue)
allsig=union(as.character(fdr005$eQTMtissue),as.character(stats$eQTMtissue[stats$lfsr<0.05]))
print('FDR<0.05 hits')
length(fdr005$eQTMtissue)
print('LFSR<0.05 hits')
length(stats$eQTMtissue[stats$lfsr<0.05])
print('FDR<0.05 + LFSR<0.05 union hits')
length(allsig)
stats=stats[stats$eQTMtissue%in%allsig,]
stats$Gene=gsub('.*:','',stats$eQTM)
stats$cpg=gsub(':.*','',stats$eQTM)

# Get q-values
minpval=read.table(paste0('results/',tissue,'.cor.stats.min.pval.txt'),header=T)

# Get full stats: HEAVY!
fullstats=as.data.frame(read_delim(paste0('results/',tissue,'.cor.stats.txt'),delim='\t'))
fullstats$Tissue=tissue

# Merge stats
merged_stats=merge(fullstats,stats,by.x=c('Tissue','cpg','gene'),by.y=c('Tissue','cpg','Gene'),all.y=T)
merged_stats=merge(merged_stats,minpval[,c('cpg','qvalue')],by.x='cpg',by.y='cpg',all.x=T)
dim(merged_stats)
merged_stats=merged_stats[,c('Tissue','cpg','gene','CpG_to_TSS_dist','rho_lower','rho','row_upper','pvalue','padj','qvalue','ttest','ttest_se','lfsr','smoothed_rho','PosteriorMean','PosteriorSD')]
write.table(file=paste0('results/',tissue,'.cor.stats.complete.txt'),quote=F,sep='\t',row.names=F,merged_stats)
quit(save='no')
