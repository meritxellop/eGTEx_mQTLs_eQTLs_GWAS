#!/usr/bin/env /apps/software/gcc-6.2.0/R/3.6.3/bin/Rscript
library(argparse)
library(readr)
library(stringr)
library(hyprcoloc)
library(reshape2)
library(permute)

parser=ArgumentParser()
parser$add_argument("-g", "--gwas", type="character", default="Astle_et_al_2016_Eosinophil_counts",
    help="GWAS")
parser$add_argument("-I", "--gwasr", type="character", default="region1004",
    help="GWAS id")
parser$add_argument("--p1", type="double", default="1e-04",
    help="p1")
parser$add_argument("--p2", type="double", default="0.98",
    help="p2")
parser$add_argument("--sample_size", type="integer",
    help="GWAS Sample size")
parser$add_argument("--cases", type="character",
    help="GWAS Cases")
parser$add_argument("--type", type="character", default="quant",
    help="GWAS type: cc or quant")
parser$add_argument("--cpg", type="character", default="cg12622986",
    help="cpg id, from top colocalizing cpg-tissue pair")
parser$add_argument("--tissue", type="character", default="MuscleSkeletal",
    help="tissue id, from top colocalizing cpg-tissue pair")
parser$add_argument("--coloc_category", type="character", default="mQTL_specific",
    help="m/eQTL-GWAS colocalization category, either mQTL specific (no eQTL colocalizing with GWAS hit) or mQTL shared (eQTL(s) and mQTL(s) colocalizing with GWAS hit)")
parser$add_argument("--eQTLtrait", type="character", default="TwinsUK_ge_blood",
    help="eQTL mapping set")
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false",
    dest="verbose", help="Print NO output")
args <- parser$parse_args()


wdir='/gpfs/data/gtex-group/mpavia/methylation/HyPrColoc'
tmpdir='/scratch/mpavia/fastQTL'

header=c('variant_id','cpg_id','tss','ma1','ma2','maf','pvaluemQTL','slope','slope_se','panel_variant_id','chromosome','position','effect_allele','non_effect_allele','current_build','frequency','sample_size','zscore','pvalue','effect_size','standard_error','imputation_status','n_cases','Molecular_trait_id','Chromosome','Position','Ref','Alt','Ma_samples','Maf','Pvalue','Beta','Se','Type','Ac','An','R2','Molecular_trait_object_id','Gene_id','Median_tpm');
file=paste0(tmpdir,'/results/',args$tissue,'/mqtls/NonPermutedResults/tmpFiles/HyPrColoc.',args$tissue,'.',args$cpg,'.',args$gwas,'.',args$gwasr,'.mQTLs.',args$eQTLtrait,'.txt')
statsall=read.table(file,header=F,fill=T,col.names=header)

colnames(statsall)=header
statsall$gwas=args$gwas
statsall$gwasr=args$gwasr

# Eliminate mQTL NA cases
statsall=statsall[!is.na(statsall$slope_se),]

# Fill GWAS effect_size and standard_error by z,N,maf
statsall$sample_size=args$sample_size
statsall$effect_size<-statsall$zscore/((args$sample_size*statsall$frequency*(1-statsall$frequency))^(1/2))
statsall$standard_error<-statsall$effect_size/statsall$zscore
statsall=statsall[!is.na(statsall$standard_error),]

# Obtain mQTL sample size
cmd=paste0('head -n 1 /gpfs/data/gtex-group/mpavia/methylation/data/Covariates/',args$tissue,'.covariates.txt | tr \'\t\' \'\n\' | grep -c \'GTEX\'')
qtl_sample_size=as.numeric(system(cmd,intern=T))

# Keep only those genes with eQTL signal compatible with GWAS signal: genes with GWAS lead variant corresponding to eQTL nominal p-value threshold
lead_gwas_variant=as.character(statsall[order(statsall$pvalue),'variant_id'][1])
print(lead_gwas_variant)
threshold=1e-02
genes=as.character(unique(statsall$Molecular_trait_id[statsall$variant_id == lead_gwas_variant & statsall$Pvalue<threshold]))
genes=genes[!is.na(genes)]

if (length(genes)<1) {
	print(paste0(file," does not contain any gene with minimal eQTL signal, GWAS lead variant P<",threshold,", to test"))
	quit(save='no')
} else {
	print(paste0(file," contains ",length(genes)," genes with eQTL signal to test: ",genes))
}

# Format data to HyprColoc standards
stats=subset(statsall,Molecular_trait_id%in%genes)
header=c('Trait','variant','effect_size')
gwas_eff_sizes=unique(stats[,c('gwas','variant_id','effect_size')]);colnames(gwas_eff_sizes)=header;
mQTL_eff_sizes=unique(stats[,c('cpg_id','variant_id','slope')]);colnames(mQTL_eff_sizes)=header
eQTL_eff_sizes=stats[,c('Molecular_trait_id','variant_id','Beta')];colnames(eQTL_eff_sizes)=header
eff_sizes=rbind(rbind(gwas_eff_sizes,mQTL_eff_sizes),eQTL_eff_sizes)
traits=unique(eff_sizes$Trait)
eff_sizes=dcast(eff_sizes,Trait~variant,value.var='effect_size')
rownames(eff_sizes)=eff_sizes[,1]
eff_sizes=eff_sizes[,-1]
eff_sizes=eff_sizes[traits,]

header=c('Trait','variant','effect_size_se')
gwas_eff_sizes_se=unique(stats[,c('gwas','variant_id','standard_error')]);colnames(gwas_eff_sizes_se)=header;
mQTL_eff_sizes_se=unique(stats[,c('cpg_id','variant_id','slope_se')]);colnames(mQTL_eff_sizes_se)=header
eQTL_eff_sizes_se=stats[,c('Molecular_trait_id','variant_id','Se')];colnames(eQTL_eff_sizes_se)=header
eff_sizes_se=rbind(rbind(gwas_eff_sizes_se,mQTL_eff_sizes_se),eQTL_eff_sizes_se)
traits=unique(eff_sizes_se$Trait)
eff_sizes_se=dcast(eff_sizes_se,Trait~variant,value.var='effect_size_se')
rownames(eff_sizes_se)=eff_sizes_se[,1]
eff_sizes_se=eff_sizes_se[,-1]
eff_sizes_se=eff_sizes_se[traits,]

# Keep genes with eQTL cis region overlapping >50% the mQTL-GWAS region
# Keep variants with GWAS, mQTL and eQTL signal for all genes selected
overlap=1-apply(eff_sizes,1,function(x) sum(is.na(x)))/ncol(eff_sizes)
eff_sizes=eff_sizes[overlap>0.50,]
eff_sizes=eff_sizes[,!apply(is.na(eff_sizes), 2, any)]
eff_sizes_se=eff_sizes_se[overlap>0.50,]
eff_sizes_se=eff_sizes_se[,!apply(is.na(eff_sizes_se), 2, any)]
if (nrow(eff_sizes)<3) {
        print(paste0(file," does not contain eQTL cis regions collectivelly overlapping >50% the mQTL-GWAS region"))
        quit(save='no')
} else { genes=paste(rownames(eff_sizes)[-1][-1],collapse=','); print(paste0("Genes to be tested ",genes)) }

# HyPrColoc
# binary.outcomes param is a binary vector of dimension the number of traits: 1 represents a binary trait 0 otherwise
Binary.outcomes=rep(0,nrow(eff_sizes))
if( args$type=="cc" ) { Binary.outcomes[1]=1 }

# The HyPrColoc variant specific priors require the specification of two parameters:
# “prior.1”, denoting the probability that a snp is associated with a single trait
# “prior.2”, such that (i) 1−prior.2 is the prior probaility that a snp is associated with an additonal trait given that it is associated with one trait; (ii) 1−(prior.2)2 is the prior probability that the snp is associated with a third trait given it is associated with two other traits; (iii) 1−(prior.2)3 is the prior probability that the snp is associated with a fourth trait given three and so on… By default prior.1=10−4, hence the prior probability that any snp is a ssociated with a single trait is 1 in 10000 (matching that used in COLOC). The default for prior.2=0.98, meaning that the prior probability that the snp is associated with a second trait, conditional on association with one trait, is 1 in 50; the prior probaility that it is associated with a third trait given association with two other traits is around 1 in 25; and so on…

# We will assess stability of clusters via a sensitivity analysis
# We choose the prior specification corresponding to allowing a  variant to have specific priors, which, for each variant, focuses on tuning the probability that a variant is colocalized with a subset of traits, for all possible subsets
# When using the variant specific prior, results tend to be most sensitive to the choice of prior.2 as this controls the prior probability of each additional colocalized trait. Hence, we focus our sensitivity analysis on varying this parameter. We consider the range of values prior.2∈{0.95,0.98,0.99,0.999}, meaning that the probability of a snp being associated with a second trait given its association with one trait can be: (i) 1 in 20, i.e. prior.2=0.95; (ii) 1 in 50, i.e. prior.2=0.98; (iii) 1 in 100, i.e. prior.2=0.99; (ii) 1 in 1000, i.e. prior.2=0.999

p1 = args$p1 # 1e-04
p2 = args$p2 # Prior of choice: 0.98

my.res.all.top=as.data.frame(t(c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)))
colnames(my.res.all.top)=c('iteration','traits','posterior_prob','regional_prob','candidate_snp','posterior_explained_by_snp','dropped_trait','nsnps','eQTL_trait','genes_tested','prior.1','prior.2','prior.c','prior.2.choice','shuftime')
for(shuftime in 1:1000){
	num=shuftime
	print(num)
	set.seed(num)
	# Randomize eQTL signal
	shuffled_ind=shuffle(1:ncol(eff_sizes))
	eff_sizes[3:nrow(eff_sizes),]=eff_sizes[3:nrow(eff_sizes),shuffled_ind]
	eff_sizes_se[3:nrow(eff_sizes_se),]=eff_sizes_se[3:nrow(eff_sizes_se),shuffled_ind]
	pc=1-p2
 	results.priorc <- hyprcoloc(as.matrix(t(eff_sizes)), as.matrix(t(eff_sizes_se)), trait.names=rownames(eff_sizes), snp.id=colnames(eff_sizes), uniform.priors = FALSE, prior.1=p1, prior.c=pc, binary.outcomes=Binary.outcomes)[[1]]
	results.priorc$nsnps=ncol(eff_sizes)
	results.priorc$eQTL_trait=args$eQTLtrait
	results.priorc$genes_tested=genes
	results.priorc$prior.1=p1
	results.priorc$prior.2=p2
	results.priorc$prior.c=pc
	results.priorc$prior.2.choice='Yes'
	results.priorc$shuftime=shuftime
	my.res.all.top=rbind(my.res.all.top,results.priorc)
}

my.res.all.top=my.res.all.top[-1,]
#save.image(file='image.test.Rimage')
file=paste0(wdir,'/results_random_eQTL/',args$coloc_category,'/',args$gwas,'.',args$gwasr,'.',args$cpg,'.',args$tissue,'.',args$eQTLtrait,'.HyPrColoc.random_eQTL.txt')
print(file)
write.table(file=file,my.res.all.top,quote=F,sep='\t',row.names=F)
quit(save="no")
