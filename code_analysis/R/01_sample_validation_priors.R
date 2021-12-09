remove(list = ls())

source('~/stan_glm_powerprior.R')
library(surbayes)
library(sas7bdat)


## 
## Get data
##
set.seed(123)
dir <- "P:/User_Workarea/Ethan_Mackenzie_Alt/Paper1/COMPASS"
compass1 <- sas7bdat::read.sas7bdat(paste0(dir, '/Phase1/blind_phase1.sas7bdat'))
compass2 <- sas7bdat::read.sas7bdat(paste0(dir, '/Phase2/blind_phase2.sas7bdat'))
compass <- data.frame(rbind(compass1, compass2))


#######################################
## params for sim
######################################
b02list <- c(0.50, 0.75, 1)
b01list <- c(0, 0.50, 0.75, 1)
a0list <- c(0, 0.50, 0.75, 1)
nfutlist <- seq(3000, 4000, length.out = 10)
B <- 10000
thin <- 5
burn <- 2000

b0grid <- expand.grid(b02 = b02list, b01 = b01list)
b0grid$x_pp_id = 1:nrow(b0grid)

## 
## Keep only relevant vars
##
yvars <- c(
  'SCORE_SIS16', 'RateHltn', 'SCORE_PROMIS'
)
xvars <- c(
  'ECAREPLAN', 'HisStk', 'HisTia', 'AGE', 'RACECAT1', 'SSCAT3', 'SELFPAY', 'STKDIAG'
)
other.vars <- c('PHASEN', 'hospclus')
keep.indx <- which(names(compass) %in% c(yvars, xvars, other.vars))
compass <- compass[, keep.indx]



## Change STKDIAG to Stroke / TIA
compass$STKDIAG <- ifelse(compass$STKDIAG == "Transient Ischemic Attack", 0, 1)


## Drop empty strings
print(paste0("Total number of obs before dropping missing: ", nrow(compass) ) )
for ( j in 1:ncol(compass) ) {
  x <- compass[, j]
  if(class(x) %in% c('character', 'factor')) {
    compass <- compass[x != "", ]
    compass[, j] <- as.factor(as.character(compass[, j]))
  }
}
compass <- compass[complete.cases(compass), ]
print(paste0("Total number of obs after dropping missing: ", nrow(compass) ) )

## make sure yvars are numeric
compass[, yvars] <- lapply(compass[, yvars], as.numeric)

## dichotomize SSCAT3
compass$SSCAT3 <- as.character(compass$SSCAT3)
compass$SSCAT3 <- with(compass, ifelse(SSCAT3 %in% c("=0", "1-4"), "0-4", "5+"))
compass$SSCAT3 <- as.factor(compass$SSCAT3)


## binary variables must be numeric
bin.vars <- c('ECAREPLAN', 'RACECAT1', 'SELFPAY', 'HisStk', 'SSCAT3')
compass[, bin.vars] <- lapply(compass[, bin.vars], as.numeric)
compass[, bin.vars] <- lapply(compass[, bin.vars], function(x) { ifelse(x == max(x), 1, 0) } )

# Separate out wake forest and compass data
wf.indx <- compass$PHASEN == 1 & compass$hospclus == 'E'
pilot <- compass[wf.indx, ]
compass <- compass[!wf.indx, ]
pilot <- pilot[complete.cases(pilot), ]
compass <- compass[complete.cases(compass), ]


##
## VALIDATION PRIOR SAMPLING
##
X.formula.list <- list(
  RACECAT1 ~ 1,
  AGE ~ RACECAT1,
  SELFPAY ~ RACECAT1 + AGE,
  HisStk ~ RACECAT1 + AGE + SELFPAY,
  SSCAT3 ~ RACECAT1 + AGE + SELFPAY + HisStk,
  STKDIAG ~ RACECAT1 + AGE + SELFPAY + HisStk + SSCAT3
)
X.family.list <- list(
  binomial(), gaussian(), binomial(), binomial(), binomial(), binomial()
)
rhs <- c('ECAREPLAN', 'AGE_C', 'RACECAT1', 'HisStk', "SSCAT3", 'STKDIAG')
rhs <- paste0(rhs, collapse = ' + ')
y.formula.list <- as.list( paste(yvars, rhs, sep = " ~ ") )
y.formula.list <- lapply(y.formula.list, as.formula)


x.sample.list = vector('list', nrow(b0grid))
for ( i in 1:nrow(b0grid) ) {
  print(paste0("i = ", i, ' of ', nrow(b0grid)))
  b01 = b0grid[i, 'b02']
  b00 = b0grid[i, 'b01']
  histdata = pilot
  if ( b00 == 0)
    histdata = NULL
  temp = vector('list', length(X.formula.list))
  for ( j in seq_along(X.formula.list) ) {
    print(paste0('j = ', j, ' of ', length(X.formula.list)))
    temp[[j]] = 
      pp_sample(
        formula = X.formula.list[[j]],
        family  = X.family.list[[j]],
        data = compass,
        histdata = histdata,
        a01 = b01,
        a00 = b00,
        iter   = 10000,
        warmup = 1000,
        thin   = 4
      )
  }
  x.sample.list[[i]] = temp
}
x.sample.list$pp_grid = b0grid
x.sample.list$x.formula.list = X.formula.list
saveRDS(x.sample.list, file = '~/x_sample.rds')


y.sample.list = vector('list', length(a0list))
pilot$AGE_C = as.numeric(scale(pilot$AGE))
compass$AGE_C = as.numeric(scale(compass$AGE))
for ( i in seq_along(a0list ) ) {
  histdata = NULL
  if(a0list[i] > 0) {
    histdata = pilot
    y.sample.list[[i]] = sur_sample(y.formula.list, data = compass, 1000 + 5 * 10000, histdata, a0 = a0list[i], burnin = 1000, thin = 5)
  } else {
    y.sample.list[[i]] = sur_sample(y.formula.list, data = compass, 10000, histdata = NULL, a0 = NULL)
  }
}
y.sample.list = list()
y.sample.list$sample = sur_sample(y.formula.list, data = compass, 10000, histdata = NULL, a0 = NULL)
y.sample.list$pp_grid = a0list
y.sample.list$y.formula.list = y.formula.list
saveRDS(y.sample.list, file = '~/y_sample_2.rds')

setwd('~/Redo')
for ( i in 1:nrow(b0grid) ) {
  smpl = x.sample.list[[i]]
  smpl$b01 = b0grid[i, 'b01']
  smpl$b02 = b0grid[i, 'b01']
  smpl$x_pp_id = b0grid[i, 'x_pp_id']
  saveRDS(smpl, paste0('x_sample_', i, '.rds'))
}

