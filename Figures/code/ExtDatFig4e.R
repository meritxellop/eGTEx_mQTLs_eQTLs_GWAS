source('code/themes.R')
library(reshape)
library(ggplot2)
library(psych)

colors=read.table('data/gtex_tissue_colors2.csv',sep=',',header=T,comment.char = '!')
rownames(colors)=gsub(' ' ,'',gsub('-','',colors$tissue_site_detail))

stats2=read.table('data/tissue.specific.mCpGs.GTEx.stats.txt')
stats.m=merge(stats,stats2[,c(1,2,5)],by.x=c(2,3),by.y=c(1,2))
stats.m$V12[stats.m$Cohort%in%'FUSION (Muscle)']=-1*(stats.m$V12[stats.m$Cohort%in%'FUSION (Muscle)'])
stats_lite$cor=NA
stats_lite$cor_CIL=NA
stats_lite$cor_CIU=NA
for (tis in unique(stats.m$tissue)) {print(tis); 
	for (coh in unique(stats.m$Cohort)) { 
		print(tis);
		print(coh);
		print(cor(subset(stats.m,tissue%in%tis & Cohort%in%coh,c('V12','V5.y')),method='spearman')[2,1]) ;
		rho=cor(subset(stats.m,tissue%in%tis & Cohort%in%coh,c('V12','V5.y')),method='spearman')[2,1]
		stats_lite$cor[stats_lite$`GTEx Tissue`%in%tis & stats_lite$Cohort%in%coh]=rho
		n=table(stats$tissue%in%tis & stats$Cohort%in%coh)[['TRUE']]
		c=r.con(rho,n,p=.95,twotailed=TRUE)
		stats_lite$cor_CIL[stats_lite$`GTEx Tissue`%in%tis & stats_lite$Cohort%in%coh]=c[1]
		stats_lite$cor_CIU[stats_lite$`GTEx Tissue`%in%tis & stats_lite$Cohort%in%coh]=c[2]
	}
}
ggplot(stats_lite,aes(x = `GTEx Tissue`, y = cor, ymin = cor_CIL, ymax = cor_CIU )) +
    geom_pointrange(aes(col = `GTEx Tissue`), size = 1) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_errorbar(aes(ymin = cor_CIL, ymax = cor_CIU, col = `GTEx Tissue`),
                  width = 0.25, cex = 1) +
    facet_wrap(.~Cohort, strip.position = "top", nrow = 1)+ scale_color_manual(values = as.character(colors[levels(stats$Tissue),'tissue_color_hex'])) + theme_Publication()+ theme(
        axis.text.x = element_blank(),
        axis.ticks = element_blank())+labs(tag = "e")+theme(plot.tag = element_text(face="bold"))+ theme(legend.position = "none")+labs(y = "Spearman Rho\n(95% Confidence Interval) ")+ scale_y_continuous(breaks=seq(-0.4,1,0.1))+ theme(legend.position = "none")+theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white") )


# 500 x 300

