library(ggplot2)

eQTLs=read.table(pipe('cut -f 1,2 data/standard_GWAS.all_eQTLs.pp4_030.rcp_030.txt | tail -n+2 | sort | uniq'),header=F)
mQTLs=read.table(pipe('cut -f 1,2 data/standard_GWAS.all_mQTLs.pp4_030.rcp_030.txt | tail -n+2 | sort | uniq'),header=F)
QTLs=unique(rbind(eQTLs,mQTLs))

# This definition does not take into account tissue-matching colocs or independent signals per region
eQTLs$eQTLs=TRUE
QTLs.m=merge(QTLs,eQTLs,all.x=T,by.x=c(1,2),by.y=c(1,2))
mQTLs$mQTLs=TRUE
QTLs.m=merge(QTLs.m,mQTLs,all.x=T,by.x=c(1,2),by.y=c(1,2))
QTLs.m$eQTLs[is.na(QTLs.m$eQTLs)]<-FALSE
QTLs.m$mQTLs[is.na(QTLs.m$mQTLs)]<-FALSE

percent_shared <-c()
percent_mQTL_specific <-c()
percent_eQTL_specific <-c()
for (trait in unique(QTLs.m$V1)) {
t.QTLs.m <- QTLs.m[grep(trait,QTLs.m$V1),]
percent_shared <- c(percent_shared,100*(nrow(subset(t.QTLs.m,eQTLs & mQTLs))/nrow(t.QTLs.m)))
percent_mQTL_specific <- c(percent_mQTL_specific,100*(nrow(subset(t.QTLs.m,!eQTLs & mQTLs))/nrow(t.QTLs.m)))
percent_eQTL_specific <- c(percent_eQTL_specific,100*(nrow(subset(t.QTLs.m,eQTLs & !mQTLs))/nrow(t.QTLs.m)))
}
percents_QTL_sharing <- cbind(unique(QTLs.m$V1),percent_shared,percent_mQTL_specific,percent_eQTL_specific)
polygenic_traits <- names(table(QTLs.m$V1))[table(QTLs.m$V1)>=10]
rownames(percents_QTL_sharing) <- percents_QTL_sharing[,1]
percents_QTL_sharing <- cbind(table(QTLs.m$V1)[polygenic_traits],percents_QTL_sharing[polygenic_traits,c(2,3,4)])
colnames(percents_QTL_sharing)[1] <- 'number_colocalized_GWAS_hits'
percents_QTL_sharing <- cbind(g[rownames(percents_QTL_sharing),'Category'],percents_QTL_sharing)
g<- read.table('data/GWAS.groups.txt',header=T,sep='\t')
abbv <- read.table('data/GWAS.abbv.txt',header=T,sep=' ')
abbv <- unique(abbv[,c(1,2)])
rownames(abbv) <- abbv$GWASid
rownames(g) <- g$Tag
percents_QTL_sharing<-as.data.frame(percents_QTL_sharing)
colnames(percents_QTL_sharing)[1] <- 'GWASGroup'
percents_QTL_sharing$GWAS <- abbv[polygenic_traits,'GWAS']
percents_QTL_sharing[,c(2,3,4,5)] <- apply(percents_QTL_sharing[,c(2,3,4,5)],2,function(x) as.numeric(x))
percents_QTL_sharing.m<-melt(percents_QTL_sharing,id.vars = c('GWAS','GWASGroup','number_colocalized_GWAS_hits'))

percents_QTL_sharing.m$GWASGroup[percents_QTL_sharing.m$GWASGroup%in%"Endocrine system disease"]<-"Other"
percents_QTL_sharing.m$GWASGroup[percents_QTL_sharing.m$GWASGroup%in%"Psychiatric_neurologic"]<-"Psychiatric and Neurologic"
percents_QTL_sharing.m$GWASGroup[percents_QTL_sharing.m$GWASGroup%in%"Allergy"]<-"Immune and Allergy"
percents_QTL_sharing.m$GWASGroup[percents_QTL_sharing.m$GWASGroup%in%"Immune"]<-"Immune and Allergy"
percents_QTL_sharing.m$GWASGroup[percents_QTL_sharing.m$GWASGroup%in%"Morphology"]<-"Other"
percents_QTL_sharing.m$GWASGroup[percents_QTL_sharing.m$GWASGroup%in%"Skeletal system disease"]<-"Other"
percents_QTL_sharing.m$GWASGroup[percents_QTL_sharing.m$GWASGroup%in%"Aging"]<-"Other"
percents_QTL_sharing.m$variable <- as.character(percents_QTL_sharing.m$variable)
percents_QTL_sharing.m$variable[percents_QTL_sharing.m$variable%in%"percent_mQTL_specific"]<-"mQTL-specific"
percents_QTL_sharing.m$variable[percents_QTL_sharing.m$variable%in%"percent_eQTL_specific"]<-"eQTL-specific"
percents_QTL_sharing.m$variable[percents_QTL_sharing.m$variable%in%"percent_shared"]<-"e/mQTL-shared"


p=ggplot(data = percents_QTL_sharing.m, aes(y = value, x = GWAS)) +
    geom_col(position='dodge',aes(fill = GWASGroup),color='black',size=0.4) +
    theme_Publication() +
    xlab(NULL)+
    ylab("\n# Colocalized GWAS hits")+
    theme(axis.text.x = element_text(angle = 50, size=7, hjust=1), axis.text.y = element_text(size=7), axis.title.y  = element_text(size=7))+
    facet_grid(variable~GWASGroup,scales='free_x',space="free") +
    scale_y_sqrt() +
    theme(strip.background = element_rect(fill='white',color='white'),strip.text = element_text(size = 8, face = "plain"))+labs(tag='a') + theme(plot.tag = element_text(face='bold'))
p
# 1275 x 450
