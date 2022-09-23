# https://stephenslab.github.io/mashr/articles/eQTL_outline.html
# https://github.com/stephenslab/gtexresults/blob/master/workflows/mashr_flashr_workflow.ipynb

library(mashr)
library(flashr)
library(mixsqp)

###### FLASH code #######

    my_init_fn <- function(Y, K = 1) {
      ret = flashr:::udv_si(Y, K)
      pos_sum = sum(ret$v[ret$v > 0])
      neg_sum = -sum(ret$v[ret$v < 0])
      if (neg_sum > pos_sum) {
        return(list(u = -ret$u, d = ret$d, v = -ret$v))
      } else
      return(ret)
    }

    flash_pipeline = function(data, ...) {
      ## current state-of-the art
      ## suggested by Jason Willwerscheid
      ## cf: discussion section of
      ## https://willwerscheid.github.io/MASHvFLASH/MASHvFLASHnn2.html
      ebnm_fn = "ebnm_ash"
      ebnm_param = list(l = list(mixcompdist = "normal",
                               optmethod = "mixSQP"),
                        f = list(mixcompdist = "+uniform",
                               optmethod = "mixSQP"))
      ##
      fl_g <- flashr:::flash_greedy_workhorse(data,
                    var_type = "constant",
                    ebnm_fn = ebnm_fn,
                    ebnm_param = ebnm_param,
                    init_fn = "my_init_fn",
                    stopping_rule = "factors",
                    tol = 1e-3,
                    verbose_output = "odF")
      fl_b <- flashr:::flash_backfit_workhorse(data,
                    f_init = fl_g,
                    var_type = "constant",
                    ebnm_fn = ebnm_fn,
                    ebnm_param = ebnm_param,
                    stopping_rule = "factors",
                    tol = 1e-3,
                    verbose_output = "odF")
      return(fl_b)
    }

    cov_flash = function(data, subset = NULL, non_singleton = FALSE, save_model = NULL) {
      if(is.null(subset)) subset = 1:mashr:::n_effects(data)
      b.center = apply(data$Bhat[subset,], 2, function(x) x - mean(x))
      ## Only keep factors with at least two values greater than 1 / sqrt(n)
      find_nonunique_effects <- function(fl) {
        thresh <- 1/sqrt(ncol(fl$fitted_values))
        vals_above_avg <- colSums(fl$ldf$f > thresh)
        nonuniq_effects <- which(vals_above_avg > 1)
        return(fl$ldf$f[, nonuniq_effects, drop = FALSE])
      }

      fmodel = flash_pipeline(b.center)
      if (non_singleton)
          flash_f = find_nonunique_effects(fmodel)
      else 
          flash_f = fmodel$ldf$f
      ## row.names(flash_f) = colnames(b)
      if (!is.null(save_model)) saveRDS(list(model=fmodel, factors=flash_f), save_model)
      if(ncol(flash_f) == 0){
        U.flash = list("tFLASH" = t(fmodel$fitted_values) %*% fmodel$fitted_values / nrow(fmodel$fitted_values))
      } else{
        U.flash = c(mashr:::cov_from_factors(t(as.matrix(flash_f)), "FLASH"),
                    list("tFLASH" = t(fmodel$fitted_values) %*% fmodel$fitted_values / nrow(fmodel$fitted_values)))
      }
      return(U.flash)
    }
###### FLASH code #######

dat = readRDS('results/CpG.gene.cors.mash.rds')

# Estimate vhat: "simple" method (using null z-scores)
vhat = estimate_null_correlation_simple(mash_set_data(dat$random.b, Shat=dat$random.s, alpha=1, zero_Bhat_Shat_reset = 1E3))

# Split data: strong vs random
data.random = mash_set_data(dat$random.b,dat$random.s,V=vhat, alpha=1, zero_Bhat_Shat_reset = 1E3)
data.strong = mash_set_data(dat$strong.b,dat$strong.s,V=vhat, alpha=1, zero_Bhat_Shat_reset = 1E3)
mash_data <- data.strong
save(file='mash_output/CpG.gene.cors.mash.simple.vhat.Robj',list=c('data.random','data.strong'))

# Compute prior matrices
# FLASH prior covariances
U.flash = cov_flash(mash_data, non_singleton = TRUE, save_model = "mash_output/CpG.gene.cors.mash.FLASH.model.rds")
save(file="mash_output/CpG.gene.cors.FLASH.Robj",U.flash)
quit(save='no')
