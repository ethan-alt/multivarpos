library(tidyverse)
library(ggthemes)
library(ggpubr)

phat.dir = '~/Projects/PoS/Redo/phat'
save.dir = '~/Projects/PoS/Redo'

filelist = list.files(phat.dir, pattern = '.rds')
nfiles   = length(filelist)




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



for ( i in 1:nfiles ) {
  phat_i = readRDS( file.path(phat.dir, paste0('phat_', i, '.rds') ) )
  phat_i = phat_i %>% mutate(
    gammastar_12 = get_gammastar(corr_12),
    gammastar_13 = get_gammastar(corr_13),
    gammastar_23 = get_gammastar(corr_23)
  )
  pos_i  = phat_i %>%
    group_by(sim_id, nfut, null_smpl, null_endpt, corr) %>%
    summarize(
      'ndatasets'   = length(bayes_1),
      'bayes_1_2.3' = mean( (bayes_1 > gamma) & ( bayes_2.3 > gammastar_23 ) ),
      'bayes_2_1.3' = mean( (bayes_2 > gamma) & ( bayes_1.3 > gammastar_13 ) ),
      'bayes_3_1.2' = mean( (bayes_3 > gamma)& ( bayes_1.2 > gammastar_12 ) ),
      'bayes_1'     = mean(bayes_1 > gamma),
      'bayes_2'     = mean(bayes_2 > gamma),
      'bayes_3'     = mean(bayes_3 > gamma),
      'bayes_1.2'   = mean(bayes_1.2 > gammastar_12),
      'bayes_1.3'   = mean(bayes_1.3 > gammastar_13),
      'bayes_2.3'   = mean(bayes_2.3 > gammastar_23),
      'freq_1'      = mean(pval_1 < alpha),
      'freq_2'      = mean(pval_2 < alpha),
      'freq_3'      = mean(pval_3 < alpha),
      'holm_1.2'    = mean( pmin(pval_1, pval_2) < alpha / 2 ),
      'holm_1.3'    = mean( pmin(pval_1, pval_3) < alpha / 2 ),
      'holm_2.3'    = mean( pmin(pval_2, pval_3) < alpha / 2 ),
      .groups = 'keep'
    )
  if ( i == 1 ) {
    res = pos_i
    next
  }
  res = rbind(res, pos_i)
}


res.long = res %>% pivot_longer(starts_with(c('bayes_', 'freq_', 'holm_')))
endpt = str_split(res.long$name, pattern = '_', n = 2)
res.long$method = sapply(endpt, function(x) x[1])
res.long$endpoint  = sapply(endpt, function(x) x[2])


endpt.mtx = rbind(c(1,2), c(1,3), c(2,3))
plotlist  = list()
for ( i in 1:nrow(endpt.mtx ) ) {
  endpt.num = endpt.mtx[i, ]
  endpt      = paste0(endpt.num, collapse = '.')
  null_endpt = paste0(endpt.num, collapse = '&')
  plotlist[[i]] = ggplot(
    data = res.long %>% filter(null_smpl == 'h1', endpoint == endpt, null_endpt == null_endpt),
    aes(x = nfut, y = value, color = method)
  ) + 
    geom_smooth(se = F) +
    facet_wrap(~corr, ncol = 1) + 
    scale_color_tableau() + 
    ggtitle(paste0(endpt.num[1], ' or ', endpt.num[2] ) )
}

ggarrange(plotlist = plotlist, ncol = 3, common.legend = T, legend = 'bottom')


saveRDS(
  list(
    'pos.wide' = res,
    'pos.long' = res.long
  ),
  file = file.path( save.dir, 'pos_v2.rds' )
)




