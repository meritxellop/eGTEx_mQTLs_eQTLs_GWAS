#!/usr/bin/env /apps/software/gcc-6.2.0/R/3.4.1/bin/Rscript
library(coloc)
library(argparse)
library(readr)
library(stringr)

parser=ArgumentParser()
parser$add_argument("-t", "--tissue", type="character", default="Lung",
    help="Tissue with Methylation and eQTL data. Path: TO-BE-COMPLETED")
parser$add_argument("-m", "--mode", type="character", default="regular",
    help="mode:regular")
parser$add_argument("-g", "--gwas", type="character", default="UKB_GIANT_2018_WHRadjBMI_COMBINED_AllAncestries",
    help="GWAS")
parser$add_argument("-p", "--gwasp", type="character", default="5-e08",
    help="GWAS threshold")
parser$add_argument("-I", "--gwasid", type="character", default="region1.chr22_7777_7778",
    help="GWAS id")
parser$add_argument("--p1", type="double", default="1e-04",
    help="p1")
parser$add_argument("--p2", type="double", default="1e-04",
    help="p2")
parser$add_argument("--p12", type="double", default="1e-05",
    help="p12")
parser$add_argument("--sample_size", type="integer",
    help="GWAS Sample size")
parser$add_argument("--cases", type="character",
    help="GWAS Cases")
parser$add_argument("--type", type="character", default="quant",
    help="GWAS type: cc or quant")
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false",
    dest="verbose", help="Print NO output")
args <- parser$parse_args()


wdir='/gpfs/data/gtex-group/mpavia/methylation/coloc'
tmpdir='/scratch/mpavia/fastQTL'

file=paste0(tmpdir,'/results/',args$tissue,'/mqtls/NonPermutedResults/tmpFiles/',args$tissue,'.',args$mode,'.',args$gwas,'.',args$gwasp,'.',args$gwasid,'.all_eQTLs.txt')
file_genes=paste0(tmpdir,'/results/',args$tissue,'/mqtls/NonPermutedResults/tmpFiles/',args$tissue,'.',args$mode,'.',args$gwas,'.',args$gwasp,'.',args$gwasid,'.genes.all_eQTLs.txt')
print(file)
statsall=read.table(file,header=F)
genes=read.table(file_genes,header=F)[,1]

header=c('variant_id','gene_id','tss','ma1','ma2','maf','pvalueeQTL','slope','slope_se','panel_variant_id','chromosome','position','effect_allele','non_effect_allele','current_build','frequency','sample_size','zscore','pvalue','effect_size','standard_error','imputation_status','n_cases');
colnames(statsall)=header

# Eliminate eQTL NA cases
statsall=statsall[!is.na(statsall$slope_se),]

# Fill GWAS effect_size and standard_error by z,N,maf
#in_effect_size_na=is.na(statsall$effect_size)
#if(sum(in_effect_size_na) > 0) {
#        statsall$sample_size=args$sample_size
#        statsall$effect_size[in_effect_size_na]<-statsall$zscore[in_effect_size_na]/((args$sample_size*statsall$frequency[in_effect_size_na]*(1-statsall$frequency[in_effect_size_na]))^(1/2))
#        in_standard_error_na=is.na(statsall$standard_error)
#        statsall$standard_error[in_standard_error_na]<-statsall$effect_size[in_standard_error_na]/statsall$zscore[in_standard_error_na]
#}
#statsall=statsall[!is.na(statsall$standard_error),]

statsall$sample_size=args$sample_size
statsall$effect_size<-statsall$zscore/((args$sample_size*statsall$frequency*(1-statsall$frequency))^(1/2))
statsall$standard_error<-statsall$effect_size/statsall$zscore
statsall=statsall[!is.na(statsall$standard_error),]

# Obtain eQTL sample size
cmd=paste0('head -n 1 /gpfs/data/gtex-group/mpavia/methylation/data/Covariates/eQTLs.',args$tissue,'.v8.covariates.txt | tr \'\t\' \'\n\' | grep -c \'GTEX\'')
qtl_sample_size=as.numeric(system(cmd,intern=T))

my.res.all=as.data.frame(t(c(NA,NA,NA,NA,NA,NA)))
colnames(my.res.all)=c('nsnps','PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf')
my.res.all.top=as.data.frame(t(c(NA,NA)))
colnames(my.res.all.top)=c('topColoc_variant','SNP.PP.H4')
genes_tested=c()
for (gene in as.character(genes)) {
	print(gene)
	stats=subset(statsall,gene_id%in%gene)
	# Minimum to perform coloc set to 50 GWAS-eQTL variants
	if (nrow(stats) < 50) {next}

	# Impute GWAS sample size:mean
	#stats$sample_size[stats$sample_size%in%NA]=mean(stats$sample_size,na.rm=T)
	# Use general GWAS sample size

	eQTL=stats[,1:9]
	GWAS=stats[,10:23]
	if( args$type=="quant" ) {
                try(my.res <- coloc.abf(dataset1=list(beta=GWAS$effect_size, varbeta=GWAS$standard_error^2,type="quant",MAF=GWAS$frequency,N=args$sample_size),
                    dataset2=list(beta=eQTL$slope, varbeta=eQTL$slope_se^2,type="quant",MAF=eQTL$maf,N=qtl_sample_size),
                    p1=args$p1,p2=args$p2,p12=args$p12))
        } else if( args$type=="cc" ) {
                args$cases=as.numeric(args$cases)
                try(my.res <- coloc.abf(dataset1=list(beta=GWAS$effect_size, varbeta=GWAS$standard_error^2,type="cc",MAF=GWAS$frequency,N=args$sample_size,s=args$cases/args$sample_size),
                    dataset2=list(beta=eQTL$slope, varbeta=eQTL$slope_se^2,type="quant",MAF=eQTL$maf,N=qtl_sample_size),
                    p1=args$p1,p2=args$p2,p12=args$p12))
        }
        if (!exists(deparse(substitute(my.res)))) {print(paste0(gene," failed"));next}
        genes_tested=c(genes_tested,gene)
#	file=paste0(tmpdir,'/results/',args$tissue,'/mqtls/NonPermutedResults/tmpFiles/',args$tissue,'.',args$mode,'.',args$gwas,'.',gene,'.Robj')
#	save(file=file,my.res)

#	file=paste0(tmpdir,'/results/',args$tissue,'/mqtls/NonPermutedResults/tmpFiles/',args$tissue,'.',args$mode,'.',args$gwas,'.',gene,'.all_eQTLs.txt')
#	write.table(file=file,as.data.frame(my.res$summary),quote=F,sep='\t')
        my.res$results$snp=gsub('SNP.','',my.res$results$snp)
        my.res$results=my.res$results[order(as.numeric(as.character(my.res$results$snp))),]
        rownames(my.res$results)=stats$variant_id
        my.res$results=my.res$results[order(my.res$results$SNP.PP.H4,decreasing=T),]
        variant=rownames(my.res$results)[1]
        names(variant)='topColoc_variant'
        SNP.PP.H4=my.res$results$SNP.PP.H4[1]
        names(SNP.PP.H4)='SNP.PP.H4'
        my.res.all.top=rbind(my.res.all.top,c(variant,SNP.PP.H4))
#	print(cbind(stats,my.res$results)[order(my.res$results$SNP.PP.H4),c('panel_variant_id','panel_variant_id','SNP.PP.H4','pvalueeQTL','pvalue')])
	
	my.res.all=rbind(my.res.all,my.res$summary)	
	rm(my.res)
}

#summary is a vector giving the number of SNPs analysed, and the posterior probabilities of H0 (no causal variant), H1 (causal variant for trait 1 only), H2 (causal variant for trait 2 only), H3 (two distinct causal variants) and H4 (one common causal variant)

#results is an annotated version of the input data containing log Approximate Bayes Factors and intermediate calculations, and the posterior probability SNP.PP.H4 of the SNP being causal for the shared signal
	my.res.all=my.res.all[-1,]
	rownames(my.res.all)=as.character(genes_tested)
	my.res.all.top=my.res.all.top[-1,]
        rownames(my.res.all.top)=as.character(genes_tested)
        my.res.all=cbind(my.res.all,my.res.all.top)
        file=paste0(wdir,'/results/',args$tissue,'/',args$gwas,'/colocsummary/',args$tissue,'.standard_eQTLs.',args$mode,'.',args$gwas,'.',args$gwasp,'.',args$gwasid,'.standard_GWAS.all_eQTLs.LD.txt')

	print(file)
	write.table(file=file,my.res.all,quote=F,sep='\t')
quit(save="no")
