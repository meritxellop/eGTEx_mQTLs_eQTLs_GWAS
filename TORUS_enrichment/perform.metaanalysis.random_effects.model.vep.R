library(metafor)

res=read.table('output/VEP/summary.vep.txt',header=T)
res$yi=res$logOR
res$vi=NA
res$sei=NA
res$sei=replmiss(res$sei,with(res,(logORciu-logORcil)/(2*1.96)))
res$vi=res$sei^2

fts=unique(res$feature)
res_summary_all=data.frame(rep(NA,10))
res$qtl=gsub('qtl','QTL',res$qtl)
qtls=unique(res$qtl)

for (QTL in qtls) {
	for (ft in fts) {
		print(QTL)
		print(ft)
		res_s=subset(res,qtl%in%QTL & feature%in%ft)
		res_summary=summary(rma(yi,vi,data=res_s))
		names(QTL)='QTL'
		res_summary_all=cbind(res_summary_all,t(cbind(QTL,ft,coef(summary(res_summary)),summary(res_summary)[23],summary(res_summary)[20])))
	}
}
res_summary_all=t(res_summary_all)[-1,]
rownames(res_summary_all)=seq(1,nrow(res_summary_all))
res_summary_all=as.data.frame(res_summary_all)
res_summary_all$padj=p.adjust(as.numeric(as.character(res_summary_all$pval)),method='bonferroni')
write.table(file='output/VEP/meta.summary.vep.txt',res_summary_all,quote=F,row.names=F,col.names=T,sep='\t')
