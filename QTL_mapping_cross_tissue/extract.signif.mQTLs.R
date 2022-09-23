library(reshape2)

load('./mash_output/mQTLs.mash.output.Robj')
rownames(mQTLs.m2$result$lfsr)=gsub(':.*','',rownames(mQTLs.m2$result$lfsr))
mQTLs.m2$result$lfsr[mQTLs.m2$result$lfsr<0]<-0
m=melt(mQTLs.m2$result$lfsr)
colnames(m)=c('gene_or_cpg','tissue','lfsr')
m=m[m$lfsr<0.05,]
write.table(file='mQTLs/lfsr.005.txt',sep='\t',row.names=F,col.names=T,quote=F,m)


