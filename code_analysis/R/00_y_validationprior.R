remove(list = ls())
library(tidyverse)
library(mvtnorm)
library(Rcpp)



B    = 10000
thin = 5
burn = 1000


sourceCpp('~/Projects/PoS/gibbs/gibbs_funs.cpp')
source('~/Projects/PoS/gibbs/gibbs_funs_r.R')


setwd('/proj/psiodalab/users/ethanalt/compass')
data = readRDS('compass.rds')
data$SCORE_PROMIS = -data$SCORE_PROMIS   ## take negative so positive values better
data$AGE_C = scale(data$AGE)             ## center/scale age

y.formula.list = list(
  SCORE_SIS16  ~ ECAREPLAN + AGE_C + RACECAT1 + HisStk + SSCAT3 + STKDIAG,
  RateHltn     ~ ECAREPLAN + AGE_C + RACECAT1 + HisStk + SSCAT3 + STKDIAG,
  SCORE_PROMIS ~ ECAREPLAN + AGE_C + RACECAT1 + HisStk + SSCAT3 + STKDIAG
)


smpl = sur_sample_gibbs(y.formula.list, data, nsmpl = B, burn = burn, thin = thin)
smpl$y.formula.list = y.formula.list
saveRDS(smpl, file = '~/Projects/PoS/Redo_DataAnalysis/y_sample.rds')




