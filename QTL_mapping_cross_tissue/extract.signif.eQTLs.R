library(reshape2)

load('./mash_output/eQTLs.mash.output.Robj')
rownames(eQTLs.m2$result$lfsr)=gsub(':.*','',rownames(eQTLs.m2$result$lfsr))
eQTLs.m2$result$lfsr[eQTLs.m2$result$lfsr<0]<-0
m=melt(eQTLs.m2$result$lfsr)
colnames(m)=c('gene_or_cpg','tissue','lfsr')
m=m[m$lfsr<0.05,]
write.table(file='eQTLs/lfsr.005.txt',sep='\t',row.names=F,col.names=T,quote=F,m)


