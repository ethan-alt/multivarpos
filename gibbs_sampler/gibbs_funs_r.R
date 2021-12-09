
library(mvtnorm)
library(Rcpp)

sur_sample_gibbs = function(
  formula.list, data, histdata = NULL, a0 = 0, nsmpl = 60000, burn = 0, thin = 1, keep.sigma = T, d0 = 2 * length(formula.list)
) {
  J       = length(formula.list)
  Xlist   = lapply(formula.list, model.matrix, data = data)
  Y       = as.matrix( data[, sapply(formula.list, function(f) all.vars(f)[1])] )
  Xwide   = do.call(cbind, Xlist)
  XtX     = crossprod(Xwide)
  XtYlist = lapply(Xlist, function(x) crossprod(x, Y))
  pj      = sapply(Xlist, ncol)
  
  ## get starting value for beta
  beta = unlist( sapply(formula.list, function(f) coef(lm(f, data))) )
  if ( ! is.null(a0) ) {
    if ( a0 < 0 | a0 > 1 )
      stop('a0 must be between 0 and 1')
  }
  
  if ( a0 > 0 ) {
    X0list    = lapply(formula.list, model.matrix, data = histdata)
    Y0        = as.matrix( histdata[, sapply(formula.list, function(f) all.vars(f)[1])] )
    X0wide    = do.call(cbind, X0list)
    X0tX0     = crossprod(X0wide)
    X0tY0list = lapply(X0list, function(x) crossprod(x, Y0))
    
    smpl = sur_sample_gibbs_pp_rcpp(
      Xwide, X0wide, Y, Y0, XtX, X0tX0, XtYlist, X0tY0list, 
      d0, beta, nsmpl, pj, burn = burn, thin = thin, keep_sigma = keep.sigma, a0 = a0
    )
  } else {
    smpl = sur_sample_gibbs_rcpp(
      Xwide, Y, XtX, XtYlist, d0, beta, nsmpl, pj, burn = burn, 
      thin = thin, keep_sigma = keep.sigma
    )
  }
  beta.names = lapply(Xlist, function(x) colnames(x))
  beta.names = sapply(1:J, function(i) paste0('eq', i, '_', beta.names[[i]]))
  colnames(smpl$beta_sample) = unlist( beta.names )
  
  smpl
}



mvlm_sample = function(
  formula.list, data, histdata = NULL, a0 = NULL, nsmpl = 60000, keep.sigma = T, d0 = 2 * length(formula.list)
) {
  Y = data[, sapply(formula.list, function(f) all.vars(f)[1])]
  Y = as.matrix(Y)
  X = model.matrix(formula.list[[1]], data = data)
  
  ## check to make sure indep vars are the same
  rhs.formula.list = as.character( sapply(formula.list, function(f) f[-2]) )
  if ( length(unique(rhs.formula.list)) > 1 )
    stop("independent variables in formula.list must be the same; try sur")
  if(is.null(a0) | is.null(histdata)) {
    smpl = mvlm_sample_rcpp(Y, X, d0, nsmpl, keep_sigma = keep.sigma)
  } else {
    if ( a0 > 0 & is.null(histdata) )
      stop("when a0 is positive, must specify histdata")
    if ( !(is.null(histdata)) & (a0 <= 0 | a0 > 1) )
      stop("when histdata is specified, a0 must be larger than 0 and no larger than 1")
    if ( a0 > 0 ) {
      Y0 = as.matrix( histdata[, sapply(formula.list, function(f) all.vars(f)[1])] )
      X0 = model.matrix(formula.list[[1]], data = histdata)
      
      smpl = mvlm_sample_pp_rcpp(Y, X, Y0, X0, d0, nsmpl, a0, keep_sigma = keep.sigma)
    }
    else {
      stop("An error occurred")
    }
  }
  colnames(smpl$beta_sample) = colnames(Y)
  rownames(smpl$beta_sample) = colnames(X)
  smpl
}




compute_phat_sameX = function(
  Y, X, Y0 = NULL, X0 = NULL, a0 = 0, d0 = 2 * (ncol(Y) + 1), trt.name
  ) {
    J       = ncol(Y)
    XtX     = crossprod(X)
    XtX_inv = chol2inv(chol(XtX))
    Bhat    = XtX_inv %*% crossprod(X,Y)
    A       = crossprod(Y - X %*% Bhat)
    n       = nrow(X)
    p       = ncol(X)
    if ( a0 > 0 ) {
      X0tX0     = crossprod(X0)
      X0tX0_inv = chol2inv(chol(X0tX0))
      Bhat0     = X0tX0_inv %*% crossprod(X0,Y0)
      A0        = crossprod(Y0 - X0 %*% Bhat0)
      n0        = nrow(X0)
      
      Omega          = chol2inv( chol( XtX + a0 * X0tX0 ) )
      Lambda         = Omega %*% (a0 * X0tX0)
      B.post         = ( diag(p) - Lambda ) %*% Bhat + Lambda %*% Bhat0
      V              = A + a0 * A0 + t(Bhat - Bhat0) %*% crossprod(Lambda, XtX) %*% (Bhat - Bhat0)
      nu.Sigma       = n + a0 * n0 + d0 - p - J - 1
      
    }
    if ( a0 == 0 ) {
      V        = A
      nu.Sigma = n + d0 - p - J - 1
      B.post   = Bhat
      Omega    = XtX_inv
    }
    rownames(B.post) = rownames(Omega) = colnames(Omega) = colnames(X)
    nu.B          = nu.Sigma - J + 1
    b1.post       = B.post[trt.name, ]
    b1.scale      = Omega[trt.name, trt.name] / nu.B * V
    b1.cor        = cov2cor(b1.scale)
    b1.cor        = b1.cor[lower.tri(b1.cor)]
    names(b1.cor) = paste0('cor_', c('12', '13', '23'))
    
    
    ## compute multivariate t probabilities
    ## marginals are nonstandard univariate t with same df
    p.indiv = pt(q = -b1.post / sqrt(diag(b1.scale)), df = nu.B, lower.tail = F)
    
    ## two-way outcomes
    combos   = cbind(c(1,2), c(1,3), c(2,3))
    p.twoway = numeric(ncol(combos))
    for ( i in 1:length(p.twoway) ) {
      temp = combos[, i]
      p.twoway[i] = 1 - mvtnorm::pmvt(
        upper = rep(0, 2), delta = b1.post[temp], df = floor(nu.B), sigma = b1.scale[temp, temp]
        , algorithm = GenzBretz(maxpts = 100000, abseps = 1e-6)
        , type = 'shifted'
      )
    }
    p.threeway = 1 - mvtnorm::pmvt(
      upper = rep(0, 3), delta = b1.post, df = floor(nu.B), sigma = b1.scale
      , algorithm = GenzBretz(maxpts = 1e6, abseps = 1e-6)
      , type = 'shifted'
    )[1]
    
    names(p.indiv)    = paste0('bayes_', 1:length(p.indiv))
    names(p.twoway)   = paste0('bayes_', c('1.2', '1.3', '2.3'))
    names(p.threeway) = 'bayes_1.2.3'
    
    c(p.indiv, p.twoway, p.threeway, b1.cor)
  }


