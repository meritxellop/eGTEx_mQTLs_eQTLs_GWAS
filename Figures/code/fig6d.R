source('code/themes.R')
counts=read.table('data/counts.summary.0.98.fpr.txt',header=T)
counts$class=gsub('mQTL_specific','GWAS hits colocalized\nw/ mQTL but no eQTL\nin the same GTEx tissue',gsub('mQTL_shared','GWAS hits colocalized\nw/ mQTL and eQTL\nin the same GTEx tissue',counts$class))
p <- ggplot(counts, aes(x=class, y=counts, fill=class)) + geom_violin(position=position_dodge(1)) +geom_boxplot(width=0.1,fill='white')+ scale_y_sqrt()
p + ylab('# colocalized eQTL contexts\nper GWAS hit') + xlab('') + scale_fill_manual(values=c('purple','red')) + theme_Publication()+theme(legend.position = 'none')+labs(tag='d')   + theme(plot.tag = element_text(face='bold')) 

# 500 x 400
