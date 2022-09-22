library(ggplot2)
source('code/themes.R')
library(reshape2)

tiers=read.table('data/pleiotropy.tiers.txt',header=T,check.names=F)
discovery=read.table('data/discovery.set.counts.txt',header=T,comment.char = '!')
tiers$color=discovery$Color
tiers$SampleS=discovery$SampleS
tiers=tiers[order(discovery$FDR005LFSR005),]
tiers$Tissue=factor(tiers$Tissue,levels=unique(tiers$Tissue))
colnames(tiers)[2:5]=c('Tier 1: Single mCpG, Single eGene','Tier 2: Single mCpG, Multiple eGenes','Tier 3: Multiple mCpGs, Single eGenes','Tier 4: Multiple mCpGs, Multiple eGenes')
tiers.m=melt(tiers[,1:5])
colnames(tiers.m)=c('Tissue','Pleiotropy Tier (Avg. % mCpGs)','% mCpGs')
percents=unlist(lapply(unique(tiers.m$'Pleiotropy Tier (Avg. % mCpGs)'),function(x) round(mean(subset(tiers.m,`Pleiotropy Tier (Avg. % mCpGs)`%in%x)$`% mCpGs`),2)))

ggplot(data = tiers.m, aes(y = `% mCpGs`, x = Tissue)) +
    geom_col(position="fill",aes(fill = Tissue, alpha = `Pleiotropy Tier (Avg. % mCpGs)`),color='black',size=0.4) +
    coord_flip() +
    theme_Publication() +
    scale_fill_manual(values = as.character(discovery$Color[order(discovery$FDR005LFSR005,decreasing=F)]),guide = FALSE) +
    labs(tag='b') +
    theme(legend.position='bottom') +
    theme(plot.tag = element_text(face='bold')) +
    ylab("% mCpGs\n") +
    scale_alpha_manual(values = c(0.1,0.25,0.65,1),labels = paste0(unique(tiers.m$`Pleiotropy Tier (Avg. % mCpGs)`),'  (',percents,')')) +
    xlab(NULL)+
    scale_y_reverse()

#450x450


