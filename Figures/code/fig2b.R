library(corrplot)

mQTL=read.table('data/cluster.masrh.mQTL.magnitude.txt')
eQTL=read.table('data/cluster.masrh.eQTL.magnitude.txt')
rownames(eQTL)=rownames(mQTL)=colnames(eQTL)=colnames(mQTL)=rownames(mQTL)=colnames(eQTL)=colnames(mQTL)=c('Breast','Colon','Kidney','Lung','Muscle','Ovary','Prostate','Testis','Blood')

colmat <- colorRampPalette(c("white","white","white","lightyellow","yellow","orange", "red"))
eQTL[is.na(eQTL)] <- 0
mQTL[is.na(mQTL)] <- 0

kk=corrplot(as.matrix(eQTL), method = 'square', type = 'lower', diag = FALSE,addCoef.col = 'black',number.cex = 0.8,tl.srt = 30,tl.col='black',cl.lim=c(0,1),col=colmat(200),order='hclust',hclust.method = 'complete',ylab='eQTLs')
text(4,-1, "Proportion of shared efects")
text(6,7, "eQTLs",cex=1.2,font=2)
kk

kk=corrplot(as.matrix(mQTL), method = 'square', type = 'lower', diag = FALSE,addCoef.col = 'black',number.cex = 0.8,tl.srt = 30,tl.col='black',cl.lim=c(0,1),col=colmat(200),order='hclust',hclust.method = 'complete',ylab='mQTLs')
text(4,-1, "Proportion of shared efects")
text(6,7, "mQTLs",cex=1.2,font=2)
text(-2,9.5, "b",font=2,cex=1.5)
kk