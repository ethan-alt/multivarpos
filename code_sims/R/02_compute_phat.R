remove(list = ls())
library(tidyverse)
library(surbayes)
library(mvPOS)
library(mvtnorm)


B = 10000
M = 60000


save.dir      = '~/Projects/PoS/Redo/phat'
y.vprior.dir  = '~/Projects/PoS/Redo/Samples/y_vprior'
x.data.dir    = '~/Projects/PoS/Redo/Data/xdata'


###############################
## Things that DO NOT change 
## between batches
###############################

## sim scenarios
nfut        = seq(200, 400, by = 20)
vprior.indx = seq(1, length(list.files(y.vprior.dir, pattern = '.rds')))
grid = expand.grid('vprior.indx' = vprior.indx, 'nfut' = nfut)

vpriordir = '~/Projects/PoS/compass_sims/corr_indcovar/vpriors'
trt.prob = 0.5


id = as.numeric( Sys.getenv('SLURM_ARRAY_TASK_ID') )
if (is.na(id)) id = 1

grid.id = grid[id, ]
y.vprior.indx = grid.id['vprior.indx']
nfut          = grid.id$nfut


## read in validation prior corresponding to indx
y.vprior = readRDS(file.path(y.vprior.dir, paste0('y_vprior_', y.vprior.indx, '.rds')))
beta.smpl        = y.vprior$beta.smpl
sigma.smpl       = y.vprior$sigma.smpl
null.smpl        = y.vprior$hyp
null.endpt       = paste0(y.vprior$null, collapse = '&')
corr.lbl         = y.vprior$corr
formula.list     = y.vprior$formula.list
rhs.formula.list = lapply(formula.list, function(f) f[-2])
ynames           = sapply(formula.list, function(f) all.vars(f)[1])
remove(list = 'y.vprior')


## store treatment names
trtnames = paste0(1:3, '.ECAREPLAN')
newnames = paste0('trt_', 1:3)


## formula list to generate data
formula.list.gen = list(
  SCORE_SIS16  ~ ECAREPLAN + SSCAT3 + RACECAT1 + AGE + SELFPAY + HisStk + STKDIAG,
  RateHltn     ~ ECAREPLAN + SSCAT3 + RACECAT1 + AGE + SELFPAY + HisStk + STKDIAG,
  SCORE_PROMIS ~ ECAREPLAN + SSCAT3 + RACECAT1 + AGE + SELFPAY + HisStk + STKDIAG
)
rhs.formula.list.gen = lapply(formula.list.gen, function(f) f[-2])

for ( b in 1:B ) {
  ## get data for covariates
  data = readRDS( file.path(x.data.dir, paste0('xdata_', b, '.rds') ) )[1:nfut, ]
  data$AGE_C = as.numeric(scale(data$AGE))
  
  ## compute means
  Mu = Matrix::bdiag( lapply( rhs.formula.list, function(x) model.matrix(x, data = data) ) )
  Mu = Mu %*% beta.smpl[b, ]
  Mu = matrix(Mu, nrow = nfut)
  Y  = Mu + rmvnorm(nfut, sigma = sigma.smpl[[b]])
  colnames(Y) = ynames
  data = cbind(Y, data)
  remove(list = c('Y', 'Mu'))
  
  ## conduct SUR sample; compute POS
  smpl = suppressMessages(
    sur_sample(formula.list, data, M = M)$betadraw[, trtnames]
  )
  smpl = data.frame(smpl)
  names(smpl) = newnames
  smpl = smpl %>%
    summarize(
      'bayes_1'   = mean(trt_1 > 0),
      'bayes_2'   = mean(trt_2 > 0),
      'bayes_3'   = mean(trt_3 > 0),
      'bayes_1.2' = mean( pmax(trt_1, trt_2) > 0 ),
      'bayes_1.3' = mean( pmax(trt_1, trt_3) > 0 ),
      'bayes_2.3' = mean( pmax(trt_2, trt_3) > 0 ),
      'corr_12'   = cor(trt_1, trt_2),
      'corr_13'   = cor(trt_1, trt_3),
      'corr_23'   = cor(trt_2, trt_3)
    )
  
  ## conduct frequentist analysis
  fitlist = lapply(formula.list, function(f) summary(lm(f, data = data)))
  tval    = sapply(fitlist, function(x) x$coefficients['ECAREPLAN', 't value'])
  df      = sapply(fitlist, function(x) x$df[2])
  pval    = pt(tval, df = df, lower.tail = F)
  names(pval) = paste0('pval_', 1:length(pval))
  
  smpl = cbind( 'sim_id' = id, 'data_id' = b, 'nfut' = nfut, 'null_smpl' = null.smpl, 'null_endpt' = null.endpt, 'corr' = corr.lbl, smpl, t(pval) )
  if ( b == 1 ) {
    res = smpl
    next
  }
  res = rbind(res, smpl)
  if ( b %% 100 == 0)
    saveRDS(res, file.path(save.dir, paste0('phat_', id, '.rds')))
}
saveRDS(res, file.path(save.dir, paste0('phat_', id, '.rds')))


