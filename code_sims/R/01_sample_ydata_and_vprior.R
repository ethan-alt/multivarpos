library(tidyverse)
library(surbayes)

remove(list = ls())
set.seed(123)

savedir = '~/Projects/PoS/compass_sims/Redo/Samples/y_vprior'
datadir = '~/Projects/PoS/compass_sims/Redo/Data'


## formula lists for posterior samples
formula.list = list(
  SCORE_SIS16  ~ ECAREPLAN + SSCAT3 + RACECAT1 + AGE_C + SELFPAY + HisStk + STKDIAG,
  RateHltn     ~ ECAREPLAN + SSCAT3 + RACECAT1 + AGE_C + SELFPAY + HisStk + STKDIAG,
  SCORE_PROMIS ~ ECAREPLAN + SSCAT3 + RACECAT1 + AGE_C + SELFPAY + HisStk + STKDIAG
)
formula.list.h0 = list(
  SCORE_SIS16  ~ SSCAT3 + RACECAT1 + AGE_C + SELFPAY + HisStk + STKDIAG,
  RateHltn     ~ SSCAT3 + RACECAT1 + AGE_C + SELFPAY + HisStk + STKDIAG,
  SCORE_PROMIS ~ SSCAT3 + RACECAT1 + AGE_C + SELFPAY + HisStk + STKDIAG
)

## parameters for Y's
beta1  = c(4.0388, 0.0333, -0.0437, 0.1123, 0.0114, -0.0001057, -0.103568, -0.095856, 0.0015638)
beta2  = c(0.6737, 0.1667, -0.0888, 0.2173, 0.0501, -0.0003, -.3489, -.2644, 0.1347)
beta3  = c(64.7912, -0.5980, 2.3036, 1.5166, -0.2990, 0.0012, 1.2206, 1.2789, -0.6305)
beta3  = -beta3
beta   = c(beta1, beta2, beta3)
vars   = c(0.07451509, 1.1203780, 110.187863)
sds    = sqrt(vars / 2)
parms  = c(beta1[2], beta2[2], beta3[2], sds)
names(parms) = c('beta_11', 'beta_12', 'beta_13', 'sd_1', 'sd_2', 'sd_3')
ynames = c('SCORE_SIS16', 'RateHltn', 'SCORE_PROMIS')


## compute mean of Y's (doesn't change between correlations)
histdata           = readRDS(file.path(datadir, 'x_histdata.rds'))
histdata$ECAREPLAN = rbinom(nrow(histdata), 1, 0.50)
histdata$AGESQ     = histdata$AGE^2
histdata$AGE_C     = as.numeric(scale(histdata$AGE))

n  = nrow(histdata)
X  = model.matrix(~ ECAREPLAN + SSCAT3 + RACECAT1 + AGE + AGESQ + SELFPAY + HisStk + STKDIAG, histdata)
Mu = cbind(X %*% beta1, X %*% beta2, X %*% beta3)
colnames(mu) = ynames



## SIMULATION PARAMETERS
nulls = list(
  c(1,2), c(1,3), c(2,3)
)
corrveclist = list(
  c(-.7, -.4, -.3),
  c(-.3, -.2, -.1),
  c(0, 0, 0),
  c(.3, .2, .1),
  c(.7, .4, .3)
)
corrindx = seq_len( length(corrveclist) )
corrindx.lbl = c('high.neg', 'low.neg', 'ind', 'low.pos', 'high.pos')

hyp = c('h0', 'h1')
grid = expand.grid(nulls = nulls, hyp = hyp, corrindx = corrindx)

corrvec2corrmat = function(corrvec) {
  res = diag(3)
  res[1,2] = res[2,1] = corrvec[1]
  res[1,3] = res[3,1] = corrvec[2]
  res[2,3] = res[3,2] = corrvec[3]
  return(res)
}

corrmat2covmat = function(corrmat, sds) {
  D = diag(sds)
  D %*% corrmat %*% D
}


corrmatlist = lapply(corrveclist, corrvec2corrmat)
Sigmalist = lapply(corrmatlist, corrmat2covmat, sds = sds)


## Generate Y's based on correlation index
Ylist = list()
for ( i in seq_along(Sigmalist) ) {
  print(paste("i =", i))
  repeat {
    Y      = Mu + mvtnorm::rmvnorm(n = n, sigma = Sigmalist[[i]])
    colnames(Y) = ynames
    histdata    = cbind(Y, histdata)
    fitlist = lapply(formula.list, function(x, data) summary(lm(formula = x, data = histdata)), data = histdata)
    betatrt.hat = sapply(fitlist, function(x) x$coefficients['ECAREPLAN', 'Estimate'])
    sigma.hat = sapply(fitlist, function(x) x$sigma)
    param.hat = c(betatrt.hat, sigma.hat)

    if ( max(abs((parms - param.hat) / parms)) < .08 ){
      Ylist[[i]] = Y
      break
    }
  }
}
saveRDS(Ylist, file.path(datadir, 'ysamples.rds'))


for ( id in seq_len(nrow(grid)) ) {
  print(paste('id = ', id))
  grid.id = grid[id, ]
  null.id = grid.id$nulls[[1]]
  hyp.id = grid.id$hyp
  corrindx.id = grid.id$corrindx
  
  data = cbind(Ylist[[corrindx.id]], histdata)
  
  fitnames = sapply(fitlist, function(x) rownames(coef(x)))
  fitnames = sapply(1:3, function(i) paste0(i, '.', fitnames[, i]))
  fitnames = c(fitnames)
  
  ## validation prior sample for power
  if(hyp.id == 'h1') {
    beta.smpl  = c()
    sigma.smpl = list()
    repeat {
      smpl       = sur_sample(formula.list, data, M = 10000)
      bad.names  = paste0(null.id, '.ECAREPLAN')
      bad.indx   = do.call(pmax, as.data.frame(smpl$betadraw[, bad.names])) < 0
      beta.smpl  = rbind(beta.smpl, smpl$betadraw[!bad.indx, ])
      sigma.smpl = c(sigma.smpl, smpl$Sigmalist[!bad.indx])
      if ( nrow(beta.smpl) >= 10000 )
        break
    }
    beta.smpl = beta.smpl[1:10000, ]
    sigma.smpl = sigma.smpl[1:10000]
  }
  
  ## validation prior sample for type 1 error
  if(hyp.id == 'h0') {
    fmla.list = formula.list.h0
    for ( j in seq_along(1:3) ) {
      if ( !(j %in% null.id) )
        fmla.list[[j]] = formula.list[[j]]
    }
    smpl = sur_sample(fmla.list, data, M = 10000)
    sigma.smpl = smpl$Sigmalist
    beta.smpl = smpl$betadraw
    beta.smpl = as.data.frame(beta.smpl)
    beta.smpl[, paste0(null.id, '.ECAREPLAN')] = 0
    beta.smpl = as.matrix(beta.smpl)
    beta.smpl = beta.smpl[, fitnames]          ## reorder columns
  }
  
  res = list(
    'beta.smpl'    = beta.smpl,
    'sigma.smpl'   = sigma.smpl,
    'hyp'          = hyp.id,
    'null'         = null.id,
    'corr'         = corrindx.lbl[corrindx.id],
    'beta1'        = beta1,
    'beta2'        = beta2,
    'beta3'        = beta3,
    'Sigma'        = Sigmalist[[corrindx.id]],
    'formula.list' = formula.list
  )
  
  saveRDS(res, file.path(savedir, paste0('y_vprior_', id, '.rds')))
}









