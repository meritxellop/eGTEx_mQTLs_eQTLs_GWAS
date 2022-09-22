library(ggplot2)
source('code/themes.R')

discovery=read.table('data/discovery.set.counts.txt',header=T,comment.char = '!')
discovery=discovery[order(discovery$FDR005LFSR005),]
discovery$Tissue=factor(discovery$Tissue,levels=discovery$Tissue)
colors=as.character(discovery$Color)
discovery$SampleSL=paste0('(N = ',discovery$SampleS,') ',discovery$FDR005LFSR005)
ggplot(data = discovery, aes(y = FDR005LFSR005, x = Tissue)) +
    geom_col(aes(fill = Tissue),color='black',size=0.4) +
    geom_text(aes(label = SampleSL), nudge_y = -55000, size = 3.5) +
    coord_flip() +
    theme_Publication() +
    scale_fill_manual(values = as.character(colors)) +
    labs(tag='A') +
    theme(plot.tag = element_text(face='bold')) + guides(fill=FALSE) + ylab("# mCpGs / 1e03\n(FDR < 0.05 U LFSR < 0.05)")+scale_y_continuous(labels=c(0,50,100,150,200),breaks=c(0,50000,100000,150000,200000))
#375x450