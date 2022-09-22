library(eulerr)
library(cowplot)

# P< 1e-05
mcpgs=read.table('data/all.mcpgs.matrix.txt',header=T,row.names=1)
p1 <- plot(euler(mcpgs), quantities = TRUE,fills = c('#E69F00','#56B4E9','#009E73','#D55E00'),main = 'All mCpGs, P< 1e-05' ,legend=T)
p2 <- plot(euler(mcpgs[rownames(mcpgs)%in%as.character(read.table('data/HM450.hg38.manifest.cpgs.txt')[,1]),]), quantities = TRUE,fills = c('#E69F00','#56B4E9','#009E73','#D55E00'),main = '450K mCpGs, P< 1e-05' ,legend=T)

print('How many GTEx exclusive mCpGs are identified at P< 1e-05?')
nrow(subset(mcpgs,GTEx & !ROSMAP & !FUSION & !GoDMC))
100*(nrow(subset(mcpgs,GTEx & !ROSMAP & !FUSION & !GoDMC))/nrow(subset(mcpgs,GTEx)))

print('How many GTEx exclusive mCpGs are identified at P< 1e-05?')
nrow(subset(mcpgs[rownames(mcpgs)%in%as.character(read.table('data/HM450.hg38.manifest.cpgs.txt')[,1]),],GTEx & !ROSMAP & !FUSION & !GoDMC))
100*(nrow(subset(mcpgs[rownames(mcpgs)%in%as.character(read.table('data/HM450.hg38.manifest.cpgs.txt')[,1]),],GTEx & !ROSMAP & !FUSION & !GoDMC))/nrow(subset(mcpgs[rownames(mcpgs)%in%as.character(read.table('data/HM450.hg38.manifest.cpgs.txt')[,1]),],GTEx)))

print('Considering GTEx mCpGs included in both arrays, only 36% (21% when considering only blood-derived mCpGs) are not included in the top-powered blood mQTL catalog (Min et al. 2021) at P < 1e-05')
HM450.cpgs <- as.character(read.table('data/HM450.hg38.manifest.cpgs.txt')[,1])
wb.GTEx.mCpGs <- read.table(pipe('grep WholeBlood data/fdr.005.txt'))[,1]
mcpgs.HM450 <- mcpgs[rownames(mcpgs)%in%HM450.cpgs,]
GoDMC.mcpgs.both_arrays <- rownames(subset(mcpgs.HM450,GoDMC))
GTEx.mcpgs.both_arrays <- rownames(subset(mcpgs.HM450,GTEx))
GTEx.wb.mcpgs.both_arrays <- GTEx.mcpgs.both_arrays[GTEx.mcpgs.both_arrays%in%wb.GTEx.mCpGs]
table(GTEx.mcpgs.both_arrays%in%GoDMC.mcpgs.both_arrays)[['FALSE']]/length(GTEx.mcpgs.both_arrays) 
table(GTEx.wb.mcpgs.both_arrays%in%GoDMC.mcpgs.both_arrays)[['FALSE']]/length(GTEx.wb.mcpgs.both_arrays)

# P< 1e-03
mcpgs.0001=read.table('data/all.mcpgs.matrix.0001.txt',header=T,row.names=1)
p3 <- plot(euler(mcpgs.0001), quantities = TRUE,fills = c('#E69F00','#56B4E9','#009E73','#D55E00'),main = 'All mCpGs, P< 1e-03' ,legend=T)
p4 <- plot(euler(mcpgs.0001[rownames(mcpgs.0001)%in%as.character(read.table('data/HM450.hg38.manifest.cpgs.txt')[,1]),]), quantities = TRUE,fills = c('#E69F00','#56B4E9','#009E73','#D55E00'),main = '450K mCpGs, P< 1e-03' ,legend=T)

plot_grid(nrow = 2,p1,p2,p3,p4)
# 700 x 450
