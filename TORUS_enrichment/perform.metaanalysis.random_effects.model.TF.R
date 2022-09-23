library(metafor)

res=read.table('output/TF/summary.TF.txt',header=T)
res$yi=res$logOR
res$vi=NA
res$sei=NA
res$sei=replmiss(res$sei,with(res,(logORciu-logORcil)/(2*1.96)))
res$vi=res$sei^2

fts=unique(res$feature)
fts=fts[!fts%in%'RBM17.1'] # not quantified in mQTLs: too little signal for TORUS to handle?
res_summary_all=data.frame(rep(NA,10))
res$qtl=gsub('qtl','QTL',res$qtl)
qtls=unique(res$qtl)

for (QTL in qtls) {
	for (ft in fts) {
		print(QTL)
		print(ft)
		res_s=subset(res,qtl%in%QTL & feature%in%ft)
		res_summary=summary(rma(yi,vi,data=res_s,control=list(maxiter=1000)))
		names(QTL)='QTL'
		res_summary_all=cbind(res_summary_all,t(cbind(QTL,ft,coef(summary(res_summary)),summary(res_summary)[23],summary(res_summary)[20])))
	}
}
res_summary_all=t(res_summary_all)[-1,]
rownames(res_summary_all)=seq(1,nrow(res_summary_all))
res_summary_all=as.data.frame(res_summary_all)
res_summary_all$padj=p.adjust(as.numeric(as.character(res_summary_all$pval)),method='bonferroni')
write.table(file='output/TF/meta.summary.TF.txt',res_summary_all,quote=F,row.names=F,col.names=T,sep='\t')
