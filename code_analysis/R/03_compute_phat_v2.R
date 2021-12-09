remove(list = ls())
library(tidyverse)
library(mvtnorm)
library(Rcpp)

sourceCpp('~/Projects/PoS/gibbs/gibbs_funs.cpp')
source('~/Projects/PoS/gibbs/gibbs_funs_r.R')


savedir = '~/Projects/PoS/Redo_DataAnalysis/phat'

b02grid     = c(0.50, 0.75, 1.00)
b01grid     = c(0.00, 0.50, 0.75, 1.00)
b0grid      = expand.grid('b02' = b02grid, 'b01' = b01grid)
a0grid   = c(0, 0.50, 1.00)
nfutgrid = seq(1000, 2000, by = 50)
grid     = expand.grid('a0' = a0grid, 'nfut' = nfutgrid, 'x_id' = 1:nrow(b0grid))
grid[, c('b02', 'b01')] = b0grid[grid$x_id, ]
grid     = grid %>%
  filter( 
      (a0 == 0.00 & b02 == 1.00 & b01 == 0.00)
    | (a0 == 0.50 & b02 == 1.00 & b01 == 0.00)
    | (a0 == 1.00 & b02 == 1.00 & b01 == 0.00)
    | (a0 == 0.50 & b02 == 1.00 & b01 == 1.00)
    | (a0 == 1.00 & b02 == 1.00 & b01 == 1.00)
    )

B    = 9000
# B    = 2
# M    = 60000
burn = 1000
thin = 2

grid = grid[order(grid$a0, grid$x_id, grid$nfut), ]
rownames(grid) = NULL

id = as.numeric( Sys.getenv('SLURM_ARRAY_TASK_ID') )
if(is.na(id)){
  id = nrow(grid)
  # id = 21
}

a0   = grid[id, 'a0']
nfut = grid[id, 'nfut']
x.id = grid[id, 'x_id']
b02  = b0grid[x.id, 'b02']
b01  = b0grid[x.id, 'b01']

remove(list = c('b02grid', 'b01grid', 'b0grid', 'a0grid', 'nfutgrid', 'grid') )

y.formula.list = list(
  SCORE_SIS16  ~ ECAREPLAN + AGE_C + RACECAT1 + HisStk + SSCAT3 + STKDIAG,
  RateHltn     ~ ECAREPLAN + AGE_C + RACECAT1 + HisStk + SSCAT3 + STKDIAG,
  SCORE_PROMIS ~ ECAREPLAN + AGE_C + RACECAT1 + HisStk + SSCAT3 + STKDIAG
)

histdata = NULL
X0 = NULL
Y0 = NULL
ynames = sapply(y.formula.list, function(f) all.vars(f)[1])
if ( a0 > 0 ) {
  setwd('/proj/psiodalab/users/ethanalt/compass')
  histdata <- readRDS('pilot.rds')
  histdata$AGE_C = as.numeric(scale(histdata$AGE))
  histdata = data.frame( lapply(histdata, function(x) { x = ifelse(is.infinite(x), NA, x); return(x) }) )
  histdata$SCORE_PROMIS = -histdata$SCORE_PROMIS
  histdata = histdata[complete.cases(histdata), ]
  
  X0 = model.matrix(y.formula.list[[1]], histdata)
  Y0 = as.matrix(histdata[, ynames])
}



main    = '~/Projects/PoS/Redo_DataAnalysis'
setwd(main)
dat.dir = file.path(main, 'x_sample', x.id)
for ( i in 1:B ) {
  ## read in new data from validation prior sample
  newdat = readRDS( file.path(dat.dir, paste0(i, '.rds') ) )
  if ( i == 1 ) {
    b01 = newdat$b01
    b02 = newdat$b02
  }
  newdat = newdat$data[1:nfut, ]
  Y = as.matrix(newdat[, ynames])
  X = model.matrix(y.formula.list[[1]], data = newdat)
  ## since all X's are the same, posterior is matrix t;
  phat = compute_phat_sameX (
    Y = Y, X = X, Y0 = Y0, X0 = X0, a0 = a0, trt.name = 'ECAREPLAN'
  )
  
  ## compute frequentist pvals
  fitlist = lapply(y.formula.list, function(f) lm(f, data = newdat) )
  df      = sapply(fitlist, function(f) f$df.residual)
  tval    = sapply(fitlist, function(f) summary(f)$coefficients['ECAREPLAN', 't value'])
  pval    = pt(tval, df = df, lower.tail = F)
  names(pval) = paste0('pval_', 1:length(pval))
  
  ## create result vector
  phat_i = data.frame(sim_id = id, xparam_id = x.id, data_id = i, b01 = b01, b02 = b02, a0 = a0, nfut = nfut, t(phat), t(pval))
  if ( i == 1 ) {
    res = phat_i
    next
  }
  res = rbind(res, phat_i)
  if ( i %% 100 == 0)
    saveRDS(res, file.path( savedir, paste0('phat_', id, '.rds') ) )
}
saveRDS(res, file.path( savedir, paste0('phat_', id, '.rds') ) )

