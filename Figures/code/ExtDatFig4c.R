#http://www.sthda.com/english/articles/32-r-graphics-essentials/132-plot-grouped-data-box-plot-bar-plot-and-more/

library(ggpubr)
res.melt<-read.table('data/magnitude_beta_by_tissue_sharing_cCRE.txt',header=T,sep='\t')
colnames(res.melt)<-c('T','Gene_Regulatory_region','Proportion_of_tissue_shared_effects')
my_comparisons <- list( c("Promoter", "Distal Enh."), c("Promoter", "Proximal Enh."), c("Proximal Enh.", "Distal Enh.") )
ggboxplot(res.melt, x = "Gene_Regulatory_region", y = "Proportion_of_tissue_shared_effects", color = "Gene_Regulatory_region", palette = "jco")+ stat_compare_means(comparisons = my_comparisons,paired = T,method = 'wilcox.test')+ stat_compare_means(label.y = 1) + xlab('') +labs(tag='c',size=20)+ theme(plot.tag = element_text(face='bold')) + font("ylab", size = 20)+ font("xy.text", size = 20) 
# 425 x 643
