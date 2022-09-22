library(ggplot2)
source('code/themes.R')
library(reshape2)
library(cowplot)

mCpGHowManyeGenes=read.table('data/Pleiotropy.mCpGHowManyeGenes.txt',header=T)
mCpGHowManyeGenes$tissue=gsub('WholeBlood','Blood',gsub('MuscleSkeletal','Muscle',gsub('KidneyCortex','Kidney',gsub('ColonTransverse','Colon',gsub('BreastMammaryTissue','Breast',mCpGHowManyeGenes$tissue)))))
discovery=read.table('data/discovery.set.counts.txt',header=T,comment.char = '!')
mCpGHowManyeGenes.df=data.frame(t(c(NA,NA,NA,NA,NA,NA)))
colnames(mCpGHowManyeGenes.df)=c(1:5,'>5')

for (Tissue in as.character(discovery$Tissue)) {
m=table(subset(mCpGHowManyeGenes,tissue%in%Tissue,'mCpGHowManyeGenes')[,1])
m=c(m[1:5],sum(m[as.numeric(names(m))>5]))
names(m)=c(1:5,'>5')
mCpGHowManyeGenes.df=rbind(mCpGHowManyeGenes.df,m)
}
mCpGHowManyeGenes.df=mCpGHowManyeGenes.df[-1,]
mCpGHowManyeGenes.df=round(100*(mCpGHowManyeGenes.df/rowSums(mCpGHowManyeGenes.df)),2)
mCpGHowManyeGenes.df$Tissue=discovery$Tissue
mCpGHowManyeGenes.df=mCpGHowManyeGenes.df[order(discovery$FDR005LFSR005),]
mCpGHowManyeGenes.df$Tissue=factor(mCpGHowManyeGenes.df$Tissue,levels=unique(mCpGHowManyeGenes.df$Tissue))
m=melt(mCpGHowManyeGenes.df)
colnames(m)=c('Tissue','# eGenes/mCpG\n(Avg. % mCpGs)','% mCpGs')

p1 <- ggplot(data = m, aes(y = `% mCpGs`, x = Tissue)) +
    geom_col(position="fill",aes(fill = Tissue, alpha = `# eGenes/mCpG\n(Avg. % mCpGs)`),color='black',size=0.4) +
    coord_flip() +
    theme_Publication() +
    scale_fill_manual(values = as.character(discovery$Color[order(discovery$FDR005LFSR005,decreasing=F)]),guide = FALSE) +
    labs(tag='C') +
    theme(legend.position='bottom') +
    theme(plot.tag = element_text(face='bold')) +
    ylab("% mCpGs") +
    scale_alpha_manual(values = c(0.1,0.35,0.50,0.65,0.80,1),labels = paste0(unique(m$`# eGenes/mCpG\n(Avg. % mCpGs)`),'  (',round(colMeans(mCpGHowManyeGenes.df[,1:6]),2),')')) +
    #xlab(NULL)+theme(axis.text.y  = element_blank()) +
    xlab(NULL) +
    scale_y_reverse()
#250x400

eGeneHowManymCpGs=read.table('data/Pleiotropy.eGeneHowManymCpGs.txt',header=T)
eGeneHowManymCpGs$tissue=gsub('WholeBlood','Blood',gsub('MuscleSkeletal','Muscle',gsub('KidneyCortex','Kidney',gsub('ColonTransverse','Colon',gsub('BreastMammaryTissue','Breast',eGeneHowManymCpGs$tissue)))))
eGeneHowManymCpGs.df=data.frame(t(c(NA,NA,NA,NA,NA,NA)))
colnames(eGeneHowManymCpGs.df)=c(1:5,'>5')

for (Tissue in as.character(discovery$Tissue)) {
m=table(subset(eGeneHowManymCpGs,tissue%in%Tissue,'eGeneHowManymCpGs')[,1])
m=c(m[1:5],sum(m[as.numeric(names(m))>5]))
names(m)=c(1:5,'>5')
eGeneHowManymCpGs.df=rbind(eGeneHowManymCpGs.df,m)
}
eGeneHowManymCpGs.df=eGeneHowManymCpGs.df[-1,]
eGeneHowManymCpGs.df=round(100*(eGeneHowManymCpGs.df/rowSums(eGeneHowManymCpGs.df)),2)
eGeneHowManymCpGs.df$Tissue=discovery$Tissue
eGeneHowManymCpGs.df=eGeneHowManymCpGs.df[order(discovery$FDR005LFSR005),]
eGeneHowManymCpGs.df$Tissue=factor(eGeneHowManymCpGs.df$Tissue,levels=unique(eGeneHowManymCpGs.df$Tissue))
m=melt(eGeneHowManymCpGs.df)
colnames(m)=c('Tissue','# mCpGs/eGene\n(Avg. % eGenes)','% eGenes')

p2 <- ggplot(data = m, aes(y = `% eGenes`, x = Tissue)) +
    geom_col(position="fill",aes(fill = Tissue, alpha = `# mCpGs/eGene\n(Avg. % eGenes)`),color='black',size=0.4) +
    coord_flip() +
    theme_Publication() +
    scale_fill_manual(values = as.character(discovery$Color[order(discovery$FDR005LFSR005,decreasing=F)]),guide = FALSE) +
    labs(tag='') +
    theme(legend.position='bottom') +
    theme(plot.tag = element_text(face='bold',color='white')) +
    ylab("% eGenes") +
    scale_alpha_manual(values = c(0.1,0.35,0.50,0.65,0.80,1),labels = paste0(unique(m$`# mCpGs/eGene\n(Avg. % eGenes)`),'  (',round(colMeans(eGeneHowManymCpGs.df[,1:6]),2),')')) +
    xlab(NULL)+theme(axis.text.y  = element_blank()) +
    scale_y_reverse()
#175x400

plot_grid(nrow = 1,p1,p2)

subm=subset(m,`# mCpGs/eGene\n(Avg. % eGenes)`%in%1)
rownames(discovery)=discovery$Tissue
cor.test(subm$`% eGenes`,discovery[as.character(subm$Tissue),'SampleS'],method='spearman')
