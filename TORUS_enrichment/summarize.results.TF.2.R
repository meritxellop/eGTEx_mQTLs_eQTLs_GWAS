library(readr)
for (TFblock in seq(1,17)) {
r=as.data.frame(read_delim(paste0('torus_annots/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.TF.',TFblock,'.txt.gz'),delim='\t'))
rownames(r)=r$SNP
r=r[,-1]

top=as.character(read.table(file='data/mqtls/signif.top.mQTL.variants.txt',header=F)[,1])
percentsm=c()

fts=colnames(r)
fts=fts[!fts%in%'RBM17_d']
for (ft in fts) {
	print(ft)
	t=table(r[top,ft])
	percentsm=c(percentsm,100*(t[['1']]/length(top)))
}
names(percentsm)=fts
percentsm=percentsm[order(names(percentsm))]

top=as.character(read.table(file='data/eqtls/signif.top.eQTL.variants.txt',header=F)[,1])
percentse=c()

for (ft in fts) {
        print(ft)
        t=table(r[top,ft])
        percentse=c(percentse,100*(t[['1']]/length(top)))
}
names(percentse)=fts
percentse=percentse[order(names(percentse))]

meta=read.table('output/TF/meta.summary.TF.txt',header=T)
nam=as.character(subset(meta,QTL%in%'mQTLs','ft')[,1])
nam=gsub('\\.1','_d',nam)
nam=nam[nam%in%names(percentsm)]
meta=meta[meta$ft%in%gsub('_d','.1',nam),]
meta$percent_top_QTL_variants=c(percentsm[nam],percentse[nam])
write.table(file=paste0('output/TF/meta.summary.TF.',TFblock,'.full.txt'),meta,row.names=F,col.names=T,sep='\t',quote=F)
}
