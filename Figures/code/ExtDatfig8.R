library(ggplot2)
library(viridis)
library(ggalt)
library(ggthemes)
library(grid)
library(cowplot)
source('code/themes.R')

m=read.table('data/UKB_50_Standing_height.colocsummary.defPrior.comparison.all.txt',header=F)
colnames(m)=c('MP','coloc_PP4_fastenloc','coloc_PP4_defPriors')

rho=as.numeric(round(cor.test(m$coloc_PP4_fastenloc, m$coloc_PP4_defPriors,method = 'spearman')[4]$estimate,2))
p1 <- ggplot(m, aes(coloc_PP4_fastenloc, coloc_PP4_defPriors)) + geom_point(shape=16, size=0.1, show.legend = FALSE,alpha=0.2,col="#440154FF") + stat_bkde2d(aes(fill=..level..), geom="polygon",bandwidth=c(0.10, 0.05)) + scale_fill_viridis()+ annotate(geom="text", x=0.9, y=0.075, label=bquote(rho==.(rho)), color="black")+ guides(fill=guide_legend(title="Density\nlevel"))+xlab("coloc PP4\nfastenloc-derived priors")+ ylab("coloc PP4\ndefault priors")+theme_Publication()+labs(tag='a') + theme(plot.tag = element_text(face='bold'))
#500x400

all_mQTLs.pp4_010.rcp_010=read.table('data/standard_GWAS.all_mQTLs.pp4_010.rcp_010.txt',header=T)

rho=as.numeric(round(cor.test(all_mQTLs.pp4_010.rcp_010$coloc_PP4, all_mQTLs.pp4_010.rcp_010$fastenloc_RCP,method = 'spearman')[4]$estimate,2))
p2 <- ggplot(all_mQTLs.pp4_010.rcp_010, aes(coloc_PP4, fastenloc_RCP)) + geom_point(col="#440154FF",shape=16, size=0.1, alpha=0.2, show.legend = FALSE) + stat_bkde2d(aes(fill=..level..), geom="polygon") + scale_fill_viridis()+ annotate(geom="text", x=0.9, y=0.075, label=bquote(rho==.(rho)), color="black")+ guides(fill=guide_legend(title="Density\nlevel"))+xlab("coloc PP4\nfastenloc-derived priors")+ ylab("fastenloc RCP")+ geom_hline(yintercept=0.3, linetype="dashed",alpha=0.8)+ geom_vline(xintercept=0.3, linetype="dashed",alpha=0.8)+theme_Publication()+labs(tag='b')+scale_x_continuous(breaks=seq(0,1,0.2))+scale_y_continuous(breaks=seq(0,1,0.2))+theme(plot.tag = element_text(face='bold'))
#500x400

plot_grid(nrow = 1,p1,p2)


#all_eQTLs.pp4_010.rcp_010=read.table('data/standard_GWAS.all_eQTLs.pp4_010.rcp_010.txt',header=T)
#eQTLs=unique(subset(all_eQTLs.pp4_010.rcp_010,coloc_PP4>0.3 & fastenloc_RCP > 0.3)[,c(1,2)])
#inter=merge(eQTLs,mQTLs,by.x=c(1,2),by.y=c(1,2))
#eQTLs.all=subset(all_eQTLs.pp4_010.rcp_010,coloc_PP4>0.3 & fastenloc_RCP > 0.3)
#mQTLs.all=subset(all_mQTLs.pp4_010.rcp_010,coloc_PP4>0.3 & fastenloc_RCP > 0.3)
#m=merge(merge(inter,mQTLs.all,by.x=c(1,2),by.y=c(1,2)),eQTLs.all,by.x=c(1,2),by.y=c(1,2))
#fastenloc.i=tapply(m$fastenloc_RCP.x,paste0(m$GWAS,':',m$GWAS_region),max)>tapply(m$fastenloc_RCP.y,paste0(m$GWAS,':',m$GWAS_region),max)
#coloc.i=tapply(m$coloc_PP4.x,paste0(m$GWAS,':',m$GWAS_region),max)>tapply(m$coloc_PP4.y,paste0(m$GWAS,':',m$GWAS_region),max)
#table(coloc.i & fastenloc.i)

