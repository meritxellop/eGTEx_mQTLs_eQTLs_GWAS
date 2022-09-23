#!/usr/bin/env /apps/software/gcc-6.2.0/R/3.5.0/bin/Rscript
library(argparse)
library(qvalue)
library(readr)

parser=ArgumentParser()
parser$add_argument("-d", "--dir", type="character", default="/gpfs/data/gtex-group/mpavia/methylation/results/permuted/",
    help="Directory (WDIR) where the results are stored: $WDIR/$TISSUE.regular.perm.txt",
    metavar="integer")
parser$add_argument("-t", "--tissue", type="character", default=FALSE,
    help="Tissue with mQTL fastqtl N=1000 permutations results. Colon_Transverse	Lung	Ovary	Prostate	Breast_Mammary_Tissue	Kidney_Cortex	Muscle_Skeletal	Testis	Whole_Blood")
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false",
    dest="verbose", help="Print NO output")
args <- parser$parse_args()

# Read mQTL fastqtl N=1000 permutations results file
filename=paste0(args$dir,args$tissue,'.regular.perm.txt')
mqtls=as.data.frame(read_delim(filename,col_names=F,del='\t'))
print(paste0('# CpGs: ',nrow(mqtls)))
q=qvalue(mqtls$X17)
print(paste0('pi1: ',1-q$pi0))
print(paste0('# mGenes ( q-value<0.05 ): ',table(q$qvalues<0.05)[['TRUE']]))
mqtls=mqtls[,c(1,7,11,14,15,13,17)]
mqtls=cbind(mqtls,q$qvalues)
colnames(mqtls)=c('cpg_id','variant_id','maf','slope','slope_se','pval_nominal','pval_permuted','qval')
filename=paste0(args$dir,args$tissue,'.regular.perm.fdr.txt')
write.table(file=filename,mqtls,quote=F,row.names=F,col.names=T,sep='\t')
head(mqtls)
quit(save='no')
