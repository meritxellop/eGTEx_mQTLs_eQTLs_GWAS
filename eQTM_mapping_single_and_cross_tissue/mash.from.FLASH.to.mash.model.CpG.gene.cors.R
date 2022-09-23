# https://stephenslab.github.io/mashr/articles/eQTL_outline.html
# https://github.com/stephenslab/gtexresults/blob/master/workflows/mashr_flashr_workflow.ipynb

library(mashr)
library(mixsqp)
library(flashr)

load("mash_output/CpG.gene.cors.FLASH.Robj")
load(file='mash_output/CpG.gene.cors.mash.simple.vhat.Robj')
mash_data <- data.strong

# SVD matrices
U.pca = cov_pca(mash_data,5)

# Empirical cov matrix
X.center = apply(mash_data$Bhat, 2, function(x) x - mean(x))

# Denoised data-driven matrices
U.ed = cov_ed(mash_data, c(U.flash, U.pca, list("XX" = t(X.center) %*% X.center / nrow(X.center))), logfile='mash_output/CpG.gene.cors.U.ed.log')

# Canonical matrices
U.can = cov_canonical(mash_data)
Ulist = c(U.ed, U.can)
save(file='mash_output/CpG.gene.cors.Ulist.Robj',Ulist)

# Re-estimate vhat
random.subset = sample(1:nrow(data.random$Bhat), min(6000, nrow(data.random$Bhat)))
random.subset = mash_set_data(data.random$Bhat[random.subset,], data.random$Shat[random.subset,], alpha=1, zero_Bhat_Shat_reset = 1E3)
vhat = estimate_null_correlation(random.subset, Ulist, max_iter = 6)
save(file='mash_output/CpG.gene.cors.vhat.Robj',vhat)

# Generate mash data with re-estimated vhat
dat = readRDS('results/CpG.gene.cors.mash.rds')
data.random = mash_set_data(dat$random.b,dat$random.s,V=vhat, alpha=1, zero_Bhat_Shat_reset = 1E3)

# mashr mixture model fitting
mash.model = mash(data.random, Ulist = Ulist , optmethod = 'mixSQP', outputlevel = 1, algorithm.version = "Rcpp")
save(file='mash_output/CpG.gene.cors.mash.mixture.model.Robj',mash.model)
quit(save='no')
