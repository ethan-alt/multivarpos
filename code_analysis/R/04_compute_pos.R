library(tidyverse)


filedir   = '~/Projects/PoS/Redo_DataAnalysis/phat'
filelist  = list.files(filedir, pattern = '.rds')

destdir  = '~/Projects/PoS/Redo_DataAnalysis/Results'
filename = 'pos_dataanalysis.rds'



## get gammastar for composite outcomes
gammastar = as.data.frame( readRDS('~/Projects/PoS/compass_sims/POSfix_trtonly/Data/gammastar_sm.rds')$dim2 )
loess.fit = loess(gammastar ~ corr, data = gammastar, span = 0.05)
corr.min      = min(gammastar$corr)
corr.max      = max(gammastar$corr)
gammastar.min = as.numeric(predict(loess.fit, newdata = data.frame(corr = corr.min)))
gammastar.max = as.numeric(predict(loess.fit, newdata = data.frame(corr = corr.max)))
get_gammastar = function(x) {
  res = predict(loess.fit, newdata = data.frame(corr = x))
  res[x < corr.min] = gammastar.min
  res[x > corr.max] = gammastar.max
  res
}


## set threshold for success
gamma = 0.95
alpha = 1 - gamma



for ( i in 1:length(filelist) ) {
  file  = readRDS( file.path( filedir, paste0('phat_', i, '.rds') ) )
  file = file %>% mutate(
    gammastar_12 = get_gammastar(cor_12),
    gammastar_13 = get_gammastar(cor_13),
    gammastar_23 = get_gammastar(cor_23)
  )
  file$succ_1_2.3 = with(file, (bayes_1 > gamma) & (bayes_2.3 > gammastar_23) )
  pos_i = file %>%
    group_by(sim_id, xparam_id, b01, b02, a0, nfut) %>%
    summarize(
        'nsims'       = length(sim_id)
      , 'bayes_1'     = mean(bayes_1 > gamma)
      , 'bayes_2'     = mean(bayes_2 > gamma)
      , 'bayes_3'     = mean(bayes_3 > gamma)
      , 'bayes_1.2'   = mean(bayes_1.2 > gammastar_12)
      , 'bayes_1.3'   = mean(bayes_1.3 > gammastar_13)
      , 'bayes_2.3'   = mean(bayes_2.3 > gammastar_23)
      , 'bayes_1_2.3' = mean( succ_1_2.3 )
      , 'freq_1'      = mean(pval_1 < alpha)
      , 'freq_2'      = mean(pval_2 < alpha)
      , 'freq_3'      = mean(pval_3 < alpha)
      , 'freq_1.2'    = mean( pmin(pval_1, pval_2) < alpha / 2 )
      , 'freq_1.3'    = mean( pmin(pval_1, pval_3) < alpha / 2 )
      , 'freq_2.3'    = mean( pmin(pval_2, pval_3) < alpha / 2 )
      , 'gate_1_2.3'  = mean( (pval_1 < alpha) & (pmin(pval_2, pval_3) < alpha / 2) )
      , .groups       = 'keep'
    )
  if ( i == 1 ) {
    res = pos_i
    next
  }
  res = rbind(res, pos_i)
}

saveRDS(res, file.path(destdir, filename))

sub = res %>% filter(
  (a0 == 0 & b01 == 0 & b02 == 1) | (a0 == 1 & b01 == 1 & b02 == 1)
)





