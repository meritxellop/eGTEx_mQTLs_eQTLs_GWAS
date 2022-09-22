library(ggplot2)
source('code/themes.R')
library(reshape2)

discoveryeQTMs=read.table('data/eQTM.tissuesharing.txt',header=T,comment.char = '!',sep='\t',check.names=F)

discoveryeQTMs$FDR005LFSR005eCpGs_per=discoveryeQTMs$FDR005LFSR005eCpGs/sum(discoveryeQTMs$FDR005LFSR005eCpGs)*100
discoveryeQTMs$FDR005LFSR005eQTMs_per=discoveryeQTMs$FDR005LFSR005eQTMs/sum(discoveryeQTMs$FDR005LFSR005eQTMs)*100
ggplot(data = discoveryeQTMs, aes(y = FDR005LFSR005eQTMs_per, x = factor(`# Tissues`))) +
    geom_col(position='dodge',color='black',size=0.4,fill='lightgrey') +
    theme_Publication() +
    labs(tag='d') +
    theme(plot.tag = element_text(face='bold')) +
    ylab("Percent of eQTMs") +
    xlab("# Tissues")
#450x459
