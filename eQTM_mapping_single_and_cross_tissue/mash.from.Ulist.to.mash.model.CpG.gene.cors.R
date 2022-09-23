# https://stephenslab.github.io/mashr/articles/mQTL_outline.html
# https://github.com/stephenslab/gtexresults/blob/master/workflows/mashr_flashr_workflow.ipynb

library(mashr)
library(mixsqp)
library(flashr)

load("mash_output/CpG.gene.cors.FLASH.Robj")
load(file='mash_output/CpG.gene.cors.mash.simple.vhat.Robj')
load(file='mash_output/CpG.gene.cors.Ulist.Robj')

# Re-estimate vhat
#random.subset = sample(1:nrow(data.random$Bhat), min(6000, nrow(data.random$Bhat)))
#random.subset = mash_set_data(data.random$Bhat[random.subset,], data.random$Shat[random.subset,], alpha=1, zero_Bhat_Shat_reset = 1E3)
#vhat = estimate_null_correlation(random.subset, Ulist, max_iter = 6)
#save(file='mash_output/CpG.gene.cors.vhat.Robj',vhat)
load(file='mash_output/CpG.gene.cors.vhat.Robj')

# Re-estimate QTL signal with re-estimated vhat
dat = readRDS('results/CpG.gene.cors.mash.rds')
data.random = mash_set_data(dat$random.b,dat$random.s,V=vhat, alpha=1, zero_Bhat_Shat_reset = 1E3)
data.strong = mash_set_data(dat$strong.b,dat$strong.s,V=vhat, alpha=1, zero_Bhat_Shat_reset = 1E3)
mash_data <- data.strong

# mashr mixture model fitting
mash.model = mash(data.random, Ulist = Ulist , optmethod = 'mixSQP', outputlevel = 1, algorithm.version = "Rcpp")
save(file='mash_output/CpG.gene.cors.mash.mixture.model.Robj',mash.model)
quit(save='no')
