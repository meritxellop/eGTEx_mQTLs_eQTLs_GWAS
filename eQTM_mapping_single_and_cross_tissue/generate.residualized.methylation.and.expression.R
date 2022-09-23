library(psych)
args = commandArgs(trailingOnly=TRUE)
tissue=args[1]
print(tissue)
TISSUE=args[2]

methylation=read.table(paste0('/gpfs/data/gtex-group/mpavia/methylation/data/Phenotypes/bedfiles/',TISSUE,'.bed.gz'),row.names=4,comment.char='!',header=T,check.names=F)
methylation=methylation[,4:ncol(methylation)]
orim=colnames(methylation)
print(length(orim))
#covs_methylation=read.table(paste0('/gpfs/data/gtex-group/mpavia/methylation/data/Covariates/',TISSUE,'.covariates.txt'),header=T,row.names=1,check.names=F)
covs_methylation=read.table(pipe(paste0('grep -P "GTEX|PEER" /gpfs/data/gtex-group/mpavia/methylation/data/Covariates/',TISSUE,'.covariates.txt')),header=T,row.names=1,check.names=F)
common_methylation=colnames(methylation)[colnames(methylation)%in%colnames(covs_methylation)]
methylation=t(resid(lm(t(methylation[,common_methylation])~t(covs_methylation[,common_methylation]))))
print(length(common_methylation))
write.table(file=paste0('/scratch/mpavia/fastQTL/results/',TISSUE,'/mqtls/NonPermutedResults/tmpFiles/',TISSUE,'.methylation.residuals.txt'),quote=F,sep='\t',row.names=T,col.names=T,methylation)
cmd=paste0("system('module load htslib; bgzip -f /scratch/mpavia/fastQTL/results/",TISSUE,"/mqtls/NonPermutedResults/tmpFiles/",TISSUE,".methylation.residuals.txt')")
eval(parse(text=cmd))

geneexp=read.table(paste0('/gpfs/data/gtex-group/v8/62552/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices/',tissue,'.v8.normalized_expression.bed.gz'),row.names=4,comment.char='!',header=T,check.names=F)
geneexp=geneexp[,4:ncol(geneexp)]
orie=colnames(geneexp)
print(length(orie))
#covs_expression=read.table(paste0('/gpfs/data/gtex-group/sex_biased_regulation_v8/data/Covariates_Data/GTEx_Analysis_v8_eQTL_covariates/',tissue,'.v8.covariates.txt'),header=T,row.names=1,check.names=F)
covs_expression=read.table(pipe(paste0('grep -P "GTEX|Inferred" /gpfs/data/gtex-group/sex_biased_regulation_v8/data/Covariates_Data/GTEx_Analysis_v8_eQTL_covariates/',tissue,'.v8.covariates.txt')),header=T,row.names=1,check.names=F)
common_expression=colnames(geneexp)[colnames(geneexp)%in%colnames(covs_expression)]
geneexp=t(resid(lm(t(geneexp[,common_expression])~t(covs_expression[,common_expression]))))
print(length(common_expression))
print(length(orim[orim%in%orie]))
common=common_expression[common_expression%in%common_methylation]
print(length(common))

write.table(file=paste0('/scratch/mpavia/fastQTL/results/',TISSUE,'/mqtls/NonPermutedResults/tmpFiles/',TISSUE,'.expression.residuals.txt'),quote=F,sep='\t',row.names=T,col.names=T,geneexp)
cmd=paste0("system('module load htslib; bgzip -f /scratch/mpavia/fastQTL/results/",TISSUE,"/mqtls/NonPermutedResults/tmpFiles/",TISSUE,".expression.residuals.txt')")
eval(parse(text=cmd))
quit(save='no')
