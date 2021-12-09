
main = '~/Projects/PoS/Redo_DataAnalysis/x_sample'
setwd(main)
xdata.list = list.files(main, pattern = '.rds')

nfut.max = 5000
id = as.numeric( Sys.getenv('SLURM_ARRAY_TASK_ID') )
if(is.na(id))
  id = 1


b0grid = expand.grid('b02' = c(0.5, 0.75, 1), 'b01' = c(0, 0.5, 0.75, 1))


x.vprior       = readRDS(paste0('x_sample_', id, '.rds'))
y.vprior       = readRDS('~/Projects/PoS/Redo_DataAnalysis/y_sample.rds')


y.formula.list = y.vprior$y.formula.list
  
  
x.formula.list = list(
  RACECAT1 ~ 1,
  AGE      ~ RACECAT1,
  SELFPAY  ~ RACECAT1 + AGE,
  HisStk   ~ RACECAT1 + AGE + SELFPAY,
  SSCAT3   ~ RACECAT1 + AGE + SELFPAY + HisStk,
  STKDIAG  ~ RACECAT1 + AGE + SELFPAY + HisStk + SSCAT3
)
x.family.list = list(
  binomial(), gaussian(), binomial(), binomial(), binomial(), binomial()
)

x.pp.id = x.vprior$x_pp_id
b01     = b0grid[x.pp.id, 'b01']
b02     = b0grid[x.pp.id, 'b02']
x.vprior = x.vprior[seq_along(x.formula.list)]

## rbind the list to simplify sampling process
x.vprior = lapply(x.vprior, function(x) do.call(rbind, x))
nsmpl    = nrow(x.vprior[[1]])

gen_glm = function(formula, family, parms, data = NULL, nfut.max = nfut.max) {
  if ( is.null(data) ) {
    data = data.frame('temp' = rep(1, nfut.max))
  }
  parms       = parms[-length(parms)]     ## remove lp__ column
  varname     = all.vars(formula)[1]
  rhs.formula = formula[-2]
  X           = model.matrix(rhs.formula, data)
  if ( family$family == 'binomial' ) {
    beta   = parms
    mean   = family$linkinv(X %*% beta)
    data[[varname]] = rbinom(nfut.max, 1, mean)
  }
  if ( family$family == 'gaussian' ) {
    sigmasq.indx = length(parms)
    sigma = sqrt(parms[sigmasq.indx])
    beta  = parms[-sigmasq.indx]
    mean  = family$linkinv(X %*% beta)
    data[[varname]] = rnorm(nfut.max, mean = mean, sd = sigma)
  }
  if ( 'temp' %in% names(data) ) {
    data = data.frame(data[[varname]])
    names(data) = varname
  }
  return(data)
}


gen_y = function(y.formula.list, beta, sigma, data) {
  y.rhs.formula.list = lapply(y.formula.list, function(f) f[-2])
  ynames             = sapply(y.formula.list, function(f) all.vars(f)[1])
  X = Matrix::bdiag( lapply(y.rhs.formula.list, function(f) model.matrix(f, data) ) )
  y = matrix(X %*% beta, ncol = 3) + mvtnorm::rmvnorm(n = nfut.max, sigma = sigma)
  colnames(y) = ynames
  # y[, 3] = -y[, 3]   ## take negative so positive values are better
  data.frame(y, data)
}



dir.create(file.path(main, id), showWarnings = FALSE)
setwd(file.path(main, id))
for ( i in 1:nsmpl ) {
  data = NULL
  for ( j in seq_along(x.formula.list) ) {
    data = gen_glm(x.formula.list[[j]], x.family.list[[j]], parms = x.vprior[[j]][i, ], data = data, nfut.max = nfut.max)
  }
  data$AGE_C = as.numeric(scale(data$AGE))
  data$ECAREPLAN = rbinom(nfut.max, 1, 0.5)
  data = gen_y(y.formula.list, y.vprior$beta_sample[i, ], y.vprior$sigma_sample[,,i], data = data)
  saveRDS(
    list(
      'data' = data,
      'b01'  = b01,
      'b02'  = b02,
      'x.id' = x.pp.id
    ),
    file = paste0(i, '.rds')
  )
}




