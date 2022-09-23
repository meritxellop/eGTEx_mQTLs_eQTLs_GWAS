library(readr)
r=as.data.frame(read_delim('torus_annots/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.complemented.txt.gz',delim='\t'))
rownames(r)=r$SNP
r=r[,-1]

top=as.character(read.table(file='data/mqtls/signif.top.mQTL.variants.txt',header=F)[,1])
percentsm=c()

fts=colnames(r)
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

meta=read.table('output/VEP/meta.summary.vep.txt',header=T)
nam=as.character(subset(meta,QTL%in%'mQTLs','ft')[,1])
nam=gsub('\\.1','_d',nam)
meta$percent_top_QTL_variants=c(percentsm[nam],percentse[nam])
write.table(file='output/VEP/meta.summary.vep.full.txt',meta,row.names=F,col.names=T,sep='\t',quote=F)
