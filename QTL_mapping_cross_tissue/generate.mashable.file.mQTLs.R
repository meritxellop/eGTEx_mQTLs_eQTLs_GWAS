dat=list()
dat$strong.b=read.table('mQTLs/beta.mQTL_top.fdr005.txt')
dat$strong.s=read.table('mQTLs/beta_se.mQTL_top.fdr005.txt')
dat$random.b=read.table('mQTLs/beta.mQTL_top.fdr005.random.txt')
dat$random.s=read.table('mQTLs/beta_se.mQTL_top.fdr005.random.txt')

names=c('strong.b','strong.s','random.b','random.s');

for (name in names) {
# convert to matrix
dat[[name]]=as.matrix(dat[[name]]);
}

# set NaN se to 0
dat[['strong.s']][which(is.nan(dat[['strong.s']]))] = 0
dat[['random.s']][which(is.nan(dat[['random.s']]))] = 0

saveRDS(file='mQTLs/mQTLs.FastQTLSumStats.mash.rds',dat)
