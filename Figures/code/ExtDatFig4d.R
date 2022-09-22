source('code/themes.R')
library(reshape)
library(ggplot2)
library(psych)

stats=read.table('data/tissue.specific.mCpGs.all.cohorts.txt',header=F)
colors=read.table('data/gtex_tissue_colors2.csv',sep=',',header=T,comment.char = '!')
rownames(colors)=gsub(' ' ,'',gsub('-','',colors$tissue_site_detail))
colnames(stats)[2]='Tissue'
stats$tissue=gsub('BreastMammaryTissue','Breast',stats$Tissue)
stats$tissue=gsub('KidneyCortex','Kidney',stats$tissue)
stats$tissue=gsub('MuscleSkeletal','Muscle',stats$tissue)
stats$tissue=gsub('WholeBlood','Blood',stats$tissue)
stats$Cohort=gsub('GoDMC','GoDMC (Blood)',gsub('FUSION','FUSION (Muscle)',gsub('ROSMAP','ROSMAP (Brain)',stats$V1)))
stats_lite=cbind(melt(tapply(abs(stats$V12),list(stats$Cohort,stats$tissue),mean)),melt(tapply(abs(stats$V12),list(stats$Cohort,stats$tissue),function(x) sd(x)/sqrt(length(x))))[,3])
colnames(stats_lite)=c('Cohort','GTEx Tissue','eff','sd')
stats_lite$`GTEx Tissue`=factor(stats_lite$`GTEx Tissue`,levels=unique(stats$tissue))
stats$Tissue <- as.factor(stats$Tissue)

p<- ggplot(stats_lite, aes(x=`GTEx Tissue`, y=abs(eff), fill=`GTEx Tissue`)) + geom_bar(stat="identity", color="black", position=position_dodge()) +geom_errorbar(aes(ymin=eff-sd, ymax=eff+sd), width=.2,
                                                                                                                                    position=position_dodge(.9))
p+labs(y = "| mQTL effect size |")+ theme_Publication() + scale_fill_manual(values = as.character(colors[levels(stats$Tissue),'tissue_color_hex'])) +facet_wrap(.~Cohort)+ theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank())+labs(tag = "d")+theme(plot.tag = element_text(face="bold"))+ theme(legend.position = "none") 

# Alternative in boxplot format.

p <- ggplot(stats, aes(x=Tissue, y=abs(V12)))  +geom_boxplot(aes(fill=Tissue))+labs(y = "| mQTL effect size |")+labs(x = "GTEx Tissue") + theme_Publication() + scale_fill_manual(values = as.character(colors[levels(stats$Tissue),'tissue_color_hex'])) +facet_wrap(.~Cohort)+ theme(axis.text.x = element_blank(),axis.ticks = element_blank())+labs(tag = "d")+theme(plot.tag = element_text(face="bold"))+ theme(legend.position = "none")+theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white") )
p

# 500 x 300
