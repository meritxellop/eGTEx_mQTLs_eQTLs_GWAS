library(pheatmap)

hoxdcpgs=read.table('data/HOXD4.ENSG00000170166.5.mCpGs.chr2_176174850_G_T_b38.mQTLs.txt',header=F)
hoxdgenes=read.table('data/HOXD.eGenes.chr2_176174850_G_T_b38.eQTLs.txt',header=F)

paletteLength <- 500
myColor <- colorRampPalette(c("navy", "white", "firebrick3"))(paletteLength)
myColor2 <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

m=as.data.frame(cbind(as.character(hoxdgenes$V1),as.character(hoxdgenes$V2),hoxdgenes$V9))
colnames(m)=c('Tissue','Gene','EffectSize')
m$EffectSize=-1*as.numeric(as.character(m$EffectSize))
mm=dcast(m,Tissue~Gene)
mm=mm[,-1]
rownames(mm)=c('Breast','Colon','Kidney','Lung','Muscle','Ovary','Prostate','Testis','Blood')
mina=abs(min(mm,na.rm=T))
maxa=abs(max(mm,na.rm=T))
maxa=max(c(mina,maxa))
mina=max(c(mina,maxa))
mina=(-1)*mina
myBreaks <- c(seq(mina, 0, length.out=ceiling(paletteLength/2) + 1),
seq(max(mm,na.rm=T)/paletteLength, maxa, length.out=floor(paletteLength/2)))
pheatmap(mm, color=myColor2, breaks=myBreaks)
#450x450

n=as.data.frame(cbind(as.character(hoxdcpgs$V1),as.character(hoxdcpgs$V2),hoxdcpgs$V9))
colnames(n)=c('Tissue','CpG','EffectSize')
n$EffectSize=-1*as.numeric(as.character(n$EffectSize))
nn=dcast(n,Tissue~CpG)
rownames(nn)=c('Breast','Colon','Kidney','Lung','Muscle','Ovary','Prostate','Testis','Blood')
nn=nn[,-1]
mina=abs(min(nn,na.rm=T))
maxa=abs(max(nn,na.rm=T))
maxa=max(c(mina,maxa))
mina=max(c(mina,maxa))
mina=(-1)*mina
myBreaks <- c(seq(mina, 0, length.out=ceiling(paletteLength/2) + 1),
seq(max(nn,na.rm=T)/paletteLength, maxa, length.out=floor(paletteLength/2)))
pheatmap(nn, color=myColor2, breaks=myBreaks,show_rownames=T,fontsize_col =7)
#1275x350
