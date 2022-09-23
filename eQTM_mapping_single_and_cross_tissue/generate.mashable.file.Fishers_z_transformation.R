library(psych)
dat=list()
strong.rho=read.table('results/rho.fdr005.txt')
random.rho=read.table('results/rho.random.txt')

# Apply Fisher's z transformation
dat$strong.b=fisherz(strong.rho)
dat$random.b=fisherz(random.rho)

# Calculate z standard error
n=c(34,75,131,32,112,44,25,37)
names(n)=colnames(dat$strong.b)

dat$strong.s=dat$strong.b
dat$random.s=dat$random.b
for (tissueindex in 1:length(n)) {
	dat$strong.s[,tissueindex][!is.na(dat$strong.s[,tissueindex])]<-(1/(n[tissueindex]-3))^(1/2)
	dat$random.s[,tissueindex][!is.na(dat$random.s[,tissueindex])]<-(1/(n[tissueindex]-3))^(1/2)
}

names=c('strong.b','strong.s','random.b','random.s');

for (name in names) {
# convert to matrix
dat[[name]]=as.matrix(dat[[name]]);
}

# set NaN se to 0
dat[['strong.s']][which(is.nan(dat[['strong.s']]))] = 0
dat[['random.s']][which(is.nan(dat[['random.s']]))] = 0

saveRDS(file='results/CpG.gene.cors.mash.rds',dat)
