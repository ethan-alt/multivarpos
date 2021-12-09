
set.seed(123)
n <- 981
library(rstanarm)
options(mc.cores = parallel::detectCores())

##
## First, generate covariates
##
   ## RACECAT1
   p <- .8838
   repeat{
     RACECAT1 = rbinom(n, 1, p)
     phat = mean(RACECAT1)
     if ( abs((p - phat) / p) < .01 )
       break
   }
   data <- data.frame(RACECAT1 = RACECAT1)

   ## AGE
   beta    = c(61.842, 8.199)
   sigma   = sqrt(156.8618)
   parms   = c(beta, sigma)
   repeat {
     data$AGE = beta[1] + beta[2] * data$RACECAT1 + rnorm(n, sd = sigma )
     fit = lm(AGE ~ RACECAT1, data = data)
     betahat  = coef(fit)
     sigmahat = summary(fit)$sigma
     parmshat = c(betahat, sigma)
     if ( max( abs( (parms - parmshat) / parms ) ) < .01 )
       break
   }

   ## SELFPAY
   beta = c(3.87712, -0.52450, -0.09903)
   repeat {
     eta = beta[1] + beta[2] * data$RACECAT1 + beta[3] * data$AGE
     p   = binomial()$linkinv(eta)
     data$SELFPAY = rbinom(n, 1, p)
     fit = glm(SELFPAY ~ RACECAT1 + AGE, data = data, family = binomial())
     betahat = coef(fit)
     if ( max( abs((beta - betahat) / beta) ) < .02 )
       break
    }

   ## HisStk
   beta = c(-2.236529, -0.299556, 0.011494, -0.862221)
   repeat {
     eta         = beta[1] + beta[2] * data$RACECAT1 + beta[3] * data$AGE + beta[4] * data$SELFPAY
     p           = binomial()$linkinv(eta)
     data$HisStk = rbinom(n, 1, p)
     fit = glm(HisStk ~ RACECAT1 + AGE + SELFPAY, data = data, family = binomial())
     betahat = coef(fit)
     if ( max( abs((beta - betahat) / beta) ) < .03 )
       break
   }

   ## SSCAT3
   beta = c(-.945750, -.307482, 0.014162, 0.215857, 0.172128)
   repeat {
     eta         = beta[1] + beta[2] * data$RACECAT1 + beta[3] * data$AGE + beta[4] * data$SELFPAY + beta[5] * data$HisStk
     p           = binomial()$linkinv(eta)
     data$SSCAT3 = rbinom(n, 1, p)
     fit = glm(SSCAT3 ~ RACECAT1 + AGE + SELFPAY + HisStk, data = data, family = binomial())
     betahat = coef(fit)
     if ( max( abs((beta - betahat) / beta) ) < .15)
       break
   }

   ## STKDIAG
   beta = c(3.186419, -0.616160, -0.029901, 0.353613, -0.300814, 1.560311)
   repeat {
     eta         = beta[1] + beta[2] * data$RACECAT1 + beta[3] * data$AGE + beta[4] * data$SELFPAY + beta[5] * data$HisStk + beta[6] * data$SSCAT3
     p           = binomial()$linkinv(eta)
     data$STKDIAG = rbinom(n, 1, p)
     fit = glm(STKDIAG ~ RACECAT1 + AGE + SELFPAY + HisStk + SSCAT3, data = data, family = binomial())
     betahat = coef(fit)
     if ( max( abs((beta - betahat) / beta) ) < .10)
       break
   }
setwd('~/Projects/PoS/compass_sims/Redo/Data')
saveRDS(data, 'x_histdata.rds')






## sample from covariate distribution
remove(list = ls())
setwd('~/Projects/PoS/compass_sims/Redo/Data')
data = readRDS('x_histdata.rds')
data$AGE_C = scale(data$AGE)
X.formula.list = list(
  RACECAT1 ~ 1,
  AGE      ~ RACECAT1,
  SELFPAY  ~ RACECAT1 + AGE_C,
  HisStk   ~ RACECAT1 + AGE_C + SELFPAY,
  SSCAT3   ~ RACECAT1 + AGE_C + SELFPAY + HisStk,
  STKDIAG  ~ RACECAT1 + AGE_C + SELFPAY + HisStk + SSCAT3
)
X.family.list = list(
  binomial(),
  gaussian(),
  binomial(),
  binomial(),
  binomial(),
  binomial()
)



fitlist = vector('list', length(X.formula.list))
for ( i in seq_along(X.formula.list) ) {
  fitlist[[i]] = stan_glm(
    X.formula.list[[i]], X.family.list[[i]], data, prior = NULL, prior_aux = NULL,
    chains = 4, iter = 12500 + 2000, warmup = 2000, thin = 5
  )
}

setwd('~/Projects/PoS/compass_sims/Redo/Samples')
saveRDS(fitlist, 'covar_samples_stan.rds')








## generate covariates
remove(list = ls())
setwd('~/Projects/PoS/compass_sims/Redo/Samples')
fitlist = readRDS('covar_samples_stan.rds')
mcmc = lapply(fitlist, function(x) as.mcmc(as.matrix(x)))
lapply(mcmc, function(x) autocorr.diag(x))

X.formula.list = list(
  RACECAT1 ~ 1,
  AGE      ~ RACECAT1,
  SELFPAY  ~ RACECAT1 + AGE_C,
  HisStk   ~ RACECAT1 + AGE_C + SELFPAY,
  SSCAT3   ~ RACECAT1 + AGE_C + SELFPAY + HisStk,
  STKDIAG  ~ RACECAT1 + AGE_C + SELFPAY + HisStk + SSCAT3
)
X.family.list = list(
  binomial(),
  gaussian(),
  binomial(),
  binomial(),
  binomial(),
  binomial()
)


X.formula.list.rhs = lapply(X.formula.list, function(f) f[-2])
nfut = 2000


## generate from prior predictive distribution and store sampled data sets
gen_priorpred = function(formula, family, smpl, data = NULL, n = 2000) {
  rhs.formula = formula[-2]
  yname       = all.vars(formula)[1]
  beta = smpl
  
  if ( family$family == 'gaussian' ) {
    beta  = smpl[-length(smpl)]
    sigma = smpl[length(smpl)]
  }
  if ( is.null(data) ) {
    X    = matrix(rep(1, n), ncol = 1)
    data = data.frame(rep(NA, times = n))
    names(data) = yname
  } else {
    X    = model.matrix(rhs.formula, data)
  }
  mu   = family$linkinv(X %*% beta)
  if ( family$family == 'binomial' ) {
    data[[yname]] = rbinom(n, 1, prob = mu)
    return(data)
  }
  data[[yname]] = rnorm(n, mean = mu, sd = sigma)
  data[[paste0(yname, '_C')]] = (data[[yname]] - mean(data[[yname]])) / sd(data[[yname]])
  return(data)
}


setwd('~/Projects/PoS/compass_sims/Redo/Data/xdata')
for ( i in 1:nrow(mcmc[[2]]) ) {
  data = NULL
  for ( j in seq_along(X.formula.list) ) {
    data = gen_priorpred(X.formula.list[[j]], X.family.list[[j]], mcmc[[j]][i, ], data = data, n = nfut)
  }
  data$AGESQ     = data$AGE^2
  data$AGESQ_C   = ( data$AGESQ - mean(data$AGESQ) ) / sd(data$AGESQ)
  data$ECAREPLAN = rbinom(nfut, size = 1, prob = 0.5)
  saveRDS(data, paste0('xdata_', i, '.rds'))
}


