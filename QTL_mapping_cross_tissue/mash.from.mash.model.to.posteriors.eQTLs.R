# https://stephenslab.github.io/mashr/articles/eQTL_outline.html
# https://github.com/stephenslab/gtexresults/blob/master/workflows/mashr_flashr_workflow.ipynb

library(mashr)
library(mixsqp)

load('mash_output/eQTLs.mash.mixture.model.Robj')
load(file='mash_output/eQTLs.vhat.Robj')

# Generate mash data with re-estimated vhat
dat = readRDS('eQTLs/eQTLs.FastQTLSumStats.mash.rds')
data.strong = mash_set_data(dat$strong.b,dat$strong.s,V=vhat, alpha=1, zero_Bhat_Shat_reset = 1E3)
mash_data <- data.strong

# Compute posterior summaries
eQTLs.m2 = mash(mash_data, g=get_fitted_g(mash.model), fixg=TRUE)
#eQTLs.m2 = mash_compute_posterior_matrices(mash.model, mash_data)
save(file='mash_output/eQTLs.mash.output.Robj',eQTLs.m2)
quit(save='no')
