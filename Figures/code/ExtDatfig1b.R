library(pheatmap)
discovery=read.table('data/discovery.set.counts.txt',header=T,comment.char = '!')

load('data/clustering_inputs.methylation.Rdata')
colnames(meth_mean_sub)=discovery$Tissue

p=pheatmap(t(meth_mean_sub[topdivergent,]),cutree_rows=2,scale='column',show_colnames=F, color = colorRampPalette(c("blue", "white", "red"))(500),clustering_method = 'complete',clustering_distance_rows=cor.methylation,fontsize_col=7, fontsize_row=7, legend = FALSE,border_color = NA)
p

load('data/clustering_inputs.expression.Rdata')
colnames(exp_mean_sub)=discovery$Tissue

p=pheatmap(t(exp_mean_sub[topdivergent,]),cutree_rows=2,scale='column',show_colnames=F, color = colorRampPalette(c("blue", "white", "red"))(500),clustering_method = 'complete',clustering_distance_rows=cor.expression,fontsize_col=7, fontsize_row=7, legend = FALSE,border_color = NA)
p


