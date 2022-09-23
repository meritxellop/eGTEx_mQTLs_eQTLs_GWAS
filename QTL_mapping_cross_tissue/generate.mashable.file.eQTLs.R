dat=list()
dat$strong.b=read.table('eQTLs/beta.eQTL_top.fdr005.all.txt')
dat$strong.s=read.table('eQTLs/beta_se.eQTL_top.fdr005.all.txt')
dat$random.b=read.table('eQTLs/beta.eQTL_top.fdr005.all.random.txt')
dat$random.s=read.table('eQTLs/beta_se.eQTL_top.fdr005.all.random.txt')

names=c('strong.b','strong.s','random.b','random.s');

for (name in names) {
# convert to matrix
dat[[name]]=as.matrix(dat[[name]]);
}

# set NaN se to 0
dat[['strong.s']][which(is.nan(dat[['strong.s']]))] = 0
dat[['random.s']][which(is.nan(dat[['random.s']]))] = 0

saveRDS(file='eQTLs/eQTLs.FastQTLSumStats.mash.rds',dat)
