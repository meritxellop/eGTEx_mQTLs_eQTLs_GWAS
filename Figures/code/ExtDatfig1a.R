library(ggplot2)
source('code/themes.R')

discovery=read.table('data/discovery.set.counts.txt',header=T,comment.char = '!')
# Methylome sample sizes, alphabetically sorted tissues
ss=read.table('data/tSNE.ss.txt')
samples=read.table('data/tSNE.samples.txt')[,1]

plotdf=read.table('data/tSNE.all.CpGs.pcaT.txt',header=T,sep=" ")

# Exclude flagged samples
plotdf$Tissue=unlist(lapply(seq(1,9),function(x) rep(discovery$Tissue[x],ss$V1[x])))
plotdf$Tissue=factor(plotdf$Tissue,levels=unique(plotdf$Tissue))
plotdf$Sample_Id=samples

toflag=c()
m=cbind(abs(scale(subset(plotdf,Tissue%in%'Breast')$V1))>2.5,abs(scale(subset(plotdf,Tissue%in%'Breast')$V2))>2.5)
rownames(m)=subset(plotdf,Tissue%in%'Breast')$Sample_Id
toflag=c(toflag,rownames(m)[m[,1] | m[,2]])

m=cbind(abs(scale(subset(plotdf,Tissue%in%'Colon')$V1))>2.5,abs(scale(subset(plotdf,Tissue%in%'Colon')$V2))>2.5)
rownames(m)=subset(plotdf,Tissue%in%'Colon')$Sample_Id
toflag=c(toflag,rownames(m)[m[,1] | m[,2]])

m=cbind(abs(scale(subset(plotdf,Tissue%in%'Kidney')$V1))>2.5,abs(scale(subset(plotdf,Tissue%in%'Kidney')$V2))>2.5)
rownames(m)=subset(plotdf,Tissue%in%'Kidney')$Sample_Id
toflag=c(toflag,rownames(m)[m[,1] | m[,2]])

m=cbind(abs(scale(subset(plotdf,Tissue%in%'Lung')$V1))>2.5,abs(scale(subset(plotdf,Tissue%in%'Lung')$V2))>2.5)
rownames(m)=subset(plotdf,Tissue%in%'Lung')$Sample_Id
toflag=c(toflag,rownames(m)[m[,1] | m[,2]])

m=cbind(abs(scale(subset(plotdf,Tissue%in%'Muscle')$V1))>2.5,abs(scale(subset(plotdf,Tissue%in%'Muscle')$V2))>2.5)
rownames(m)=subset(plotdf,Tissue%in%'Muscle')$Sample_Id
toflag=c(toflag,rownames(m)[m[,1] | m[,2]])

m=cbind(abs(scale(subset(plotdf,Tissue%in%'Ovary')$V1))>2.5,abs(scale(subset(plotdf,Tissue%in%'Ovary')$V2))>2.5)
rownames(m)=subset(plotdf,Tissue%in%'Ovary')$Sample_Id
toflag=c(toflag,rownames(m)[m[,1] | m[,2]])

m=cbind(abs(scale(subset(plotdf,Tissue%in%'Prostate')$V1))>2.5,abs(scale(subset(plotdf,Tissue%in%'Prostate')$V2))>2.5)
rownames(m)=subset(plotdf,Tissue%in%'Prostate')$Sample_Id
toflag=c(toflag,rownames(m)[m[,1] | m[,2]])

m=cbind(abs(scale(subset(plotdf,Tissue%in%'Testis')$V1))>2.5,abs(scale(subset(plotdf,Tissue%in%'Testis')$V2))>2.5)
rownames(m)=subset(plotdf,Tissue%in%'Testis')$Sample_Id
toflag=c(toflag,rownames(m)[m[,1] | m[,2]])

m=cbind(abs(scale(subset(plotdf,Tissue%in%'Blood')$V1))>2.5,abs(scale(subset(plotdf,Tissue%in%'Blood')$V2))>2.5)
rownames(m)=subset(plotdf,Tissue%in%'Blood')$Sample_Id
toflag=c(toflag,rownames(m)[m[,1] | m[,2]])

ggplot(subset(plotdf,!Sample_Id%in%toflag), aes(x = V1, y = V2)) +
geom_point(shape=21,colour = "black", aes(fill=Tissue), size = 4, stroke = 0.5) +
scale_fill_manual(values = as.character(discovery$Color),guide=F) +
theme_Publication() +
xlab('t-SNE 1') +
ylab('t-SNE 2') +
annotate(size=5,geom="text", x=43, y=67, label="Whole\nBlood") +
annotate(size=5,geom="text", x=-88, y=20, label="Muscle\nSkeletal") +
annotate(size=5,geom="text", x=20, y=30, label="Lung") +
annotate(size=5,geom="text", x=-35, y=30, label="Colon\nTransverse") +
annotate(size=5,geom="text", x=75, y=20, label="Kidney\nCortex") +
annotate(size=5,geom="text", x=5, y=-20, label="Ovary") +
annotate(size=5,geom="text", x=-15, y=65, label="Prostate") +
annotate(size=5,geom="text", x=10, y=95, label="Testis") +
annotate(size=5,geom="text", x=-48, y=-45, label="Breast\nMammary") +
labs(tag='a') +
theme(plot.tag = element_text(face='bold'))
#1000x375
