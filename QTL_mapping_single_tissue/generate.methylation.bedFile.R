#!/usr/bin/env /apps/software/gcc-6.2.0/R/3.5.0/bin/Rscript
library(argparse)
library(peer)
library(readr)
library(stringr)

parser=ArgumentParser()
parser$add_argument("-t", "--tissue", type="character", default=FALSE,
    help="Tissue with Methylation data. Path: wdir/data/Phenotypes/normalized/TISSUE.noob_final_BMIQ_9-24-2019.RData")
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false",
    dest="verbose", help="Print NO output")
args <- parser$parse_args()

wdir='/gpfs/data/gtex-group/mpavia/methylation'

# Read CpG metadata
filename=paste0(wdir,'/data/Annotations/MethylationEPIC_v-1-0_B4.filtered.hg38.bed')
cmd=paste('cut -d";" -f 1 ',filename)
cpg_annots=read.table(pipe(cmd),header=F,row.names=1)
colnames(cpg_annots)=c('#chr','start','end','gene_id')

# Read methylation data
filename=paste0(wdir,'/data/Phenotypes/normalized/',args$tissue,'.noob_final_BMIQ_9-24-2019.RData')
load(filename)
metobject=ls()[5];
metdata=eval(as.name(metobject));
metobject
dim(metdata) # CpGs x samples
colnames(metdata)=str_extract(colnames(metdata),'GTEX-\\w+')
selected=rownames(cpg_annots)[rownames(cpg_annots)%in%rownames(metdata)]
metdata=metdata[selected,];

# Rank-transform and quantile normalize
print('Rank transform');
metdata=qnorm(t(apply(metdata, 1, rank, ties.method = "average"))/ (ncol(metdata)+1));
metdata=cbind(cpg_annots[rownames(metdata),],metdata);
dim(metdata) # CpGs x samples

print('Write file');
write.table(metdata,file=paste(c(wdir,"/data/Phenotypes/bedfiles/",args$tissue,".bed"),collapse=''),quote=F,row.names=F,col.names=T,sep='\t');
quit(save='no')
