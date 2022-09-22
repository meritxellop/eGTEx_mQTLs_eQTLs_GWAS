library(ggplot2)
source('code/themes.R')

discoveryeQTMs=read.table('data/eQTM.discovery.txt',header=T,comment.char = '!')
discovery=read.table('data/discovery.set.counts.txt',header=T,comment.char = '!')[-3,] # Kidney is excluded from eQTM analysis
discovery=discovery[order(discoveryeQTMs$SampleS),]
discovery$Tissue=factor(discovery$Tissue,levels=discovery$Tissue)
colors=as.character(discovery$Color)

discoveryeQTMs$SampleSL=paste0('(N = ',discoveryeQTMs$SampleS,') ',discoveryeQTMs$FDR005LFSR005)
discoveryeQTMs=discoveryeQTMs[order(discoveryeQTMs$SampleS),]
discoveryeQTMs$Tissue=factor(discoveryeQTMs$Tissue,levels=discoveryeQTMs$Tissue)

ggplot(data = discoveryeQTMs, aes(y = FDR005LFSR005, x = Tissue)) +
    geom_col(aes(fill = Tissue),color='black',size=0.4) +
    geom_text(aes(label = SampleSL), nudge_y = -1250, size = 3.5) +
    coord_flip() +
    theme_Publication() +
    scale_fill_manual(values = as.character(colors)) +
    labs(tag='a') +
    theme(plot.tag = element_text(face='bold')) + 
    guides(fill=FALSE) +
    ylab("# eQTMs\n(FDR < 0.05 U LFSR < 0.05)")# +
#492x459
