library(ggplot2)
source('code/themes.R')
library(reshape2)

mQTLs.lfsr.long=read.table('data/mQTLs.lfsr.005.txt',header=T)
mQTLs.fdr.long=read.table('data/mQTLs.fdr.005.txt',header=T)
discovery=read.table('data/discovery.set.counts.txt',header=T,comment.char = '!')

mQTLs.fdr.m=dcast(mQTLs.fdr.long,gene_or_cpg~tissue)
rownames(mQTLs.fdr.m)=mQTLs.fdr.m[,1]
mQTLs.fdr.m=mQTLs.fdr.m[,-1]
mQTLs.fdr.m[!is.na(mQTLs.fdr.m)]<-1
mQTLs.fdr.m[is.na(mQTLs.fdr.m)]<-0

mQTLs.lfsr.m=dcast(mQTLs.lfsr.long,gene_or_cpg~tissue)
rownames(mQTLs.lfsr.m)=mQTLs.lfsr.m[,1]
mQTLs.lfsr.m=mQTLs.lfsr.m[,-1]
mQTLs.lfsr.m[!is.na(mQTLs.lfsr.m)]<-1
mQTLs.lfsr.m[is.na(mQTLs.lfsr.m)]<-0

cpgs=unique(c(rownames(mQTLs.lfsr.m),rownames(mQTLs.fdr.m)))
mQTLs.fdr.m=merge(cpgs,mQTLs.fdr.m,by.x=1,by.y=0,all.x=T)[,-1]
mQTLs.fdr.m[is.na(mQTLs.fdr.m)]<-0
mQTLs.lfsr.m=merge(cpgs,mQTLs.lfsr.m,by.x=1,by.y=0,all.x=T)[,-1]
mQTLs.lfsr.m[is.na(mQTLs.lfsr.m)]<-0
mQTLs=as.matrix(mQTLs.lfsr.m)+as.matrix(mQTLs.fdr.m)
mQTLs[mQTLs>1]<-1
rownames(mQTLs)=cpgs

cpgc=rowSums(mQTLs)
names(cpgc)=cpgs

m=as.data.frame(c(table(genec)/sum(table(genec)),table(cpgc)/sum(table(cpgc))))
m$QTL=c(rep('eQTLs',9),rep('mQTLs',9))
m$num=rep(seq(1,9),2)
colnames(m)=c('Percent of QTLs','QTL','# Tissues')
m$QTL=factor(m$QTL,levels=c('mQTLs','eQTLs'))
ggplot(data = m, aes(y = (`Percent of QTLs`)*100, x = factor(`# Tissues`))) +
    geom_col(position='dodge',aes(fill = QTL),color='black',size=0.4) +
    scale_fill_manual(values=c('red','blue')) +
    theme_Publication() +
    labs(tag='a') +
    theme(plot.tag = element_text(face='bold')) +
    ylab("Percent of QTLs") +
    xlab("# Tissues")
#600x400


# Are mQTLs more tissue-shared that eQTLs?
wilcox.test(cpgc,genec) 


