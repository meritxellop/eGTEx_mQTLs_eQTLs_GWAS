library(ggplot2)
source('themes.R')
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

cpgsl=list()
cpgsm=data.frame(t(rep(NA,3)))
colnames(cpgsm)=c('Tissue','# Tissues\n(Avg. %\nmCpGs)','counts')
for(num in seq(1,9)) {
cpgsl[[num]]=cpgs[mQTLs[,num]>0]
m=cbind(rep(discovery$Tissue[num],9),as.data.frame(table(cpgc[as.character(cpgsl[[num]])])))
colnames(m)=c('Tissue','# Tissues\n(Avg. %\nmCpGs)','counts')
cpgsm=rbind(cpgsm,m)
}
cpgsm=cpgsm[-1,]
cpgsm$Tissue=factor(cpgsm$Tissue,levels=discovery$Tissue[order(discovery$FDR005LFSR005,decreasing=F)])
percents=unlist(lapply(seq(1,9),function(x) round(100*mean(subset(cpgsm,`# Tissues\n(Avg. %\nmCpGs)`%in%x)$counts/unlist(lapply(seq(1,9),function(x) length(cpgsl[[x]])))),2)))

p1 <- ggplot(data = cpgsm, aes(y = counts, x = Tissue)) +
    geom_col(position="fill",aes(fill = Tissue, alpha = `# Tissues\n(Avg. %\nmCpGs)`),color='black',size=0.4) +
    coord_flip() +
    theme_Publication() +
    scale_fill_manual(values = as.character(discovery$Color[order(discovery$FDR005LFSR005,decreasing=F)]),labels=levels(cpgsm$Tissue),guide = guide_legend(reverse = TRUE)) +
    labs(tag='') +
    theme(plot.tag = element_text(face='bold')) +
    ylab("% mCpGs\n(FDR < 0.05 U LFSR < 0.05)\n") +
    scale_alpha_manual(values = c(0.01,0.075,0.1,0.2,0.3,0.4,0.45,0.65,1),labels = paste0(seq(1,9),'  (',percents,')')) +
    xlab(NULL)+theme(axis.text.y  = element_blank())


#400x450

## Correlation of tissue-sharing with sample size

cor.test(100*(subset(cpgsm,`# Tissues\n(Avg. %\nmCpGs)`%in%1)$counts/unlist(lapply(seq(1,9),function(x) length(cpgsl[[x]])))),discovery$SampleS,method='spearman')

eQTLs.lfsr.long=read.table('eQTLs.lfsr.005.txt',header=T)
eQTLs.fdr.long=read.table('eQTLs.fdr.005.txt',header=T)
discovery=read.table('discovery.set.counts.txt',header=T,comment.char = '!')

eQTLs.fdr.m=dcast(eQTLs.fdr.long,gene_or_cpg~tissue)
rownames(eQTLs.fdr.m)=eQTLs.fdr.m[,1]
eQTLs.fdr.m=eQTLs.fdr.m[,-1]
eQTLs.fdr.m[!is.na(eQTLs.fdr.m)]<-1
eQTLs.fdr.m[is.na(eQTLs.fdr.m)]<-0

eQTLs.lfsr.m=dcast(eQTLs.lfsr.long,gene_or_cpg~tissue)
rownames(eQTLs.lfsr.m)=eQTLs.lfsr.m[,1]
eQTLs.lfsr.m=eQTLs.lfsr.m[,-1]
eQTLs.lfsr.m[!is.na(eQTLs.lfsr.m)]<-1
eQTLs.lfsr.m[is.na(eQTLs.lfsr.m)]<-0

genes=unique(c(rownames(eQTLs.lfsr.m),rownames(eQTLs.fdr.m)))
eQTLs.fdr.m=merge(genes,eQTLs.fdr.m,by.x=1,by.y=0,all.x=T)[,-1]
eQTLs.fdr.m[is.na(eQTLs.fdr.m)]<-0
eQTLs.lfsr.m=merge(genes,eQTLs.lfsr.m,by.x=1,by.y=0,all.x=T)[,-1]
eQTLs.lfsr.m[is.na(eQTLs.lfsr.m)]<-0
eQTLs=as.matrix(eQTLs.lfsr.m)+as.matrix(eQTLs.fdr.m)
eQTLs[eQTLs>1]<-1
rownames(eQTLs)=genes

genec=rowSums(eQTLs)
names(genec)=genes

genesl=list()
genesm=data.frame(t(rep(NA,3)))
colnames(genesm)=c('Tissue','# Tissues\n   (Avg. % eGenes)','counts')
for(num in seq(1,9)) {
  genesl[[num]]=genes[eQTLs[,num]>0]
  m=cbind(rep(discovery$Tissue[num],9),as.data.frame(table(genec[as.character(genesl[[num]])])))
  colnames(m)=c('Tissue','# Tissues\n   (Avg. % eGenes)','counts')
  genesm=rbind(genesm,m)
}
genesm=genesm[-1,]
genesm$Tissue=factor(genesm$Tissue,levels=discovery$Tissue[order(discovery$FDR005LFSR005)])
percents=unlist(lapply(seq(1,9),function(x) round(100*mean(subset(genesm,`# Tissues\n   (Avg. % eGenes)`%in%x)$counts/unlist(lapply(seq(1,9),function(x) length(genesl[[x]])))),2)))

p2 <- ggplot(data = genesm, aes(y = counts, x = Tissue)) +
  geom_col(position="fill",aes(fill = Tissue, alpha = `# Tissues\n   (Avg. % eGenes)`),color='black',size=0.4) +
  coord_flip() +
  theme_Publication() +
  scale_fill_manual(values = as.character(discovery$Color[order(discovery$FDR005LFSR005)])) +
  #scale_fill_manual(values = c('red','blue')) +
  labs(tag='b') +
  theme(plot.tag = element_text(face='bold')) +
  ylab("% eGenes\n(FDR < 0.05 U LFSR < 0.05)") +
  scale_alpha_manual(values = c(0.01,0.075,0.1,0.2,0.3,0.4,0.45,0.65,1),labels = paste0(seq(1,9),'  (',percents,')'))+xlab(NULL)+theme(axis.text.y  = element_blank())
#450x450

## Correlation of tissue-sharing with sample size

genesm=cbind(genesm,unlist(lapply(seq(1,9),function(x) round(100*(subset(genesm,`# Tissues\n   (Avg. % eGenes)`%in%x)$counts/unlist(lapply(seq(1,9),function(x) length(genesl[[x]])))),2))))
colnames(genesm)[4]='mean_percent_sharing'
sampleS=discovery$SampleS
names(sampleS)=discovery$Tissue
cor.test(genesm$mean_percent_sharing,sampleS[genesm$Tissue],method='spearman')