library(eulerr)
library(ggplot2)
source('code/themes.R')


eQTLs=read.table(pipe('cut -f 1,2 standard_GWAS.all_eQTLs.pp4_030.rcp_030.txt | tail -n+2 | sort | uniq'),header=F)
mQTLs=read.table(pipe('cut -f 1,2 standard_GWAS.all_mQTLs.pp4_030.rcp_030.txt | tail -n+2 | sort | uniq'),header=F)
QTLs=unique(rbind(eQTLs,mQTLs))

# This definition does not take into accoiunt tissue-matching colocs or independent signals per region
eQTLs$eQTLs=TRUE
QTLs.m=merge(QTLs,eQTLs,all.x=T,by.x=c(1,2),by.y=c(1,2))
mQTLs$mQTLs=TRUE
QTLs.m=merge(QTLs.m,mQTLs,all.x=T,by.x=c(1,2),by.y=c(1,2))
QTLs.m$eQTLs[is.na(QTLs.m$eQTLs)]<-FALSE
QTLs.m$mQTLs[is.na(QTLs.m$mQTLs)]<-FALSE

plot(euler(QTLs.m[, 3:4]), quantities = TRUE,fills = c('blue','red'))
# Save in pdf as is

#https://ggplot2-book.org/facet.html

results=read.table('data/GWAS.abbv.txt',header=T,comment.char = '!')
colnames(results)[4]="# Colocalized GWAS hits"
tags=read.table('data/GWAS.groups.txt',sep='\t',row.names=2)
results=unique(merge(rbind(cbind(as.character(results$GWASid),as.character(results$GWAS),'mQTLs'),cbind(as.character(results$GWASid),as.character(results$GWAS),'eQTLs')),results,all.x=T,by.x=c(1,2,3),by.y=c(1,2,3)))
results$`# Colocalized GWAS hits`[is.na(results$`# Colocalized GWAS hits`)]<-0
colnames(results)[1:3]=c('GWASid','GWAS','QTL')
results$QTL=factor(results$QTL,levels=c('mQTLs','eQTLs'))
results=results[order(results$`# Colocalized GWAS hits`,decreasing=T),]
results$GWAS=factor(results$GWAS,levels=unique(results$GWAS))
results$Group=tags[as.character(results$GWASid),]
results2=subset(results,!Group%in%c('Anthropometric','Blood Cell Counts','Cardiometabolic'))
results=subset(results,Group%in%c('Anthropometric','Blood Cell Counts','Cardiometabolic'))
results$Group=factor(results$Group,levels=unique(results$Group))
results2$Group=factor(results2$Group,levels=unique(results2$Group))

p=ggplot(data = results, aes(y = `# Colocalized GWAS hits`, x = GWAS)) +
    geom_col(position='dodge',aes(fill = QTL),color='black',size=0.4) +
    theme_Publication() +
    scale_fill_manual(values = c('red','blue')) +
    xlab(NULL)+
theme(axis.text.x = element_text(angle = 50, size=7, hjust=1), axis.text.y = element_text(size=7), axis.title.y  = element_text(size=7))+
    theme(legend.direction = "horizontal", legend.position = "top",legend.title=element_blank()) +
    labs(tag='A') +
    theme(plot.tag = element_text(face='bold')) + 
    facet_grid(.~Group,scales='free_x',space="free") +
    scale_y_sqrt() +
    theme(strip.background = element_rect(fill='white',color='white'),strip.text = element_text(size = 12, face = "plain")) 
p
#900x375

results2$Group=as.character(results2$Group)
results2$Group[results2$Group%in%"Endocrine system disease"]<-"Other"
results2$Group[results2$Group%in%"Psychiatric_neurologic"]<-"Psychiatric and Neurologic"
results2$Group[results2$Group%in%"Allergy"]<-"Immune and Allergy"
results2$Group[results2$Group%in%"Immune"]<-"Immune and Allergy"
results2$Group[results2$Group%in%"Morphology"]<-"Other"
results2$Group[results2$Group%in%"Skeletal system disease"]<-"Other"
results2$Group[results2$Group%in%"Aging"]<-"Other"
results2$Group=factor(results2$Group,levels=c('Immune and Allergy','Psychiatric and Neurologic','Cancer','Other'))

p=ggplot(data = results2, aes(y = `# Colocalized GWAS hits`, x = GWAS)) +
    geom_col(position='dodge',aes(fill = QTL),color='black',size=0.4) +
    theme_Publication() +
    scale_fill_manual(values = c('red','blue'),guide=FALSE) +
    xlab(NULL)+
    ylab("\n# Colocalized GWAS hits")+
theme(axis.text.x = element_text(angle = 50, size=7, hjust=1), axis.text.y = element_text(size=7), axis.title.y  = element_text(size=7))+
    facet_grid(.~Group,scales='free_x',space="free") +
    scale_y_sqrt() +
    theme(strip.background = element_rect(fill='white',color='white'),strip.text = element_text(size = 12, face = "plain"))
p
#900x350
