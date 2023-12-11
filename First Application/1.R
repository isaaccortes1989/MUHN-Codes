rm(list=ls())
library(extraDistr)
library(latex2exp)
library(nortest)
library(moments)
library(LambertW)

y <- c(0.0903296,  0.2036540, 0.2043140, 0.2808870, 0.1976530, 0.3286410,
       0.1486220,  0.1623940, 0.2627270, 0.1794550, 0.3266350, 0.2300810,
       0.1833120,  0.1509440, 0.2000710, 0.1918020, 0.1541920, 0.4641250,
       0.1170630,  0.1481410, 0.1448100, 0.1330830, 0.2760160, 0.4204770, 
       0.1224170,  0.2285950, 0.1138520, 0.2252140, 0.1769690, 0.2007440,
       0.1670450,  0.2316230, 0.2910290, 0.3412730, 0.4387120, 0.2626510, 
       0.1896510,  0.1725670, 0.2400770, 0.3116460, 0.1635860, 0.1824530,
       0.1641270,  0.1534810, 0.1618650, 0.2760160, 0.2538320, 0.2004470)

#=====================#
#=======TABLE 2=======#
#=====================#

length(y)
summary(y)
sd(y)
skewness(y)
kurtosis(y)

#===========================================================#
#===Run the code lines for each distribution independently==#
#===========================================================#


#================================================#
#===UL Estimates, SEs and Criteria===============#
#================================================#

mUL   <- function(y,sigma){
llike <- 2*sigma-3*log(1-y)-log(1+exp(sigma))-((exp(sigma)*y)/(1-y)) 
l     <- -sum(llike)
         return(l)
}

mUL2  <- function(y,sigma){
llike <- 2*log(sigma)-3*log(1-y)-log(1+sigma)-((sigma*y)/(1-y)) 
l     <- -sum(llike)
return(l)
}

UL1    <- optim(c(0), y = y, fn = mUL, method = "L-BFGS-B", lower = c(-Inf), upper = c(Inf), hessian = TRUE)
UL2    <- optim(exp(UL1$par[1]), y = y, fn = mUL2, method = "L-BFGS-B", lower = c(0), upper = c(Inf), hessian = TRUE)

AIC1       <- 2*UL2$value+2*1
BIC1       <- 2*UL2$value+1*log(length(y))
estimates1 <- UL2$par
se1        <- sqrt(diag(solve(UL2$hessian)))

#=============#
#===RESULTS===#
#=============#
estimates1
se1
AIC1
BIC1


#================================================#
#===UHN Estimates, SEs and Criteria==============#
#================================================#

mUHN <- function(y,sigma){
  llike <- log(2)-sigma-2*log(1-y)+dnorm(y/(exp(sigma)*(1-y)),log = TRUE) 
  l     <- -sum(llike)
  return(l)
}

mUHN2 <- function(y,sigma){
  llike <- log(2)-log(sigma)-2*log(1-y)+dnorm(y/(sigma*(1-y)),log = TRUE) 
  l     <- -sum(llike)
  return(l)
}

UH1    <- optim(c(0),y = y,fn = mUHN, method = "L-BFGS-B",lower = c(-Inf), upper = c(Inf), hessian = TRUE)
UH2    <- optim(exp(UH1$par[1]),y = y,fn = mUHN2, method = "L-BFGS-B",lower = c(0), upper = c(Inf), hessian = TRUE)

AIC2       <- 2*UH2$value+2*1
BIC2       <- 2*UH2$value+1*log(length(y))
estimates2 <- UH2$par
se2        <- sqrt(diag(solve(UH2$hessian)))

#=============#
#===RESULTS===#
#=============#
estimates2
se2
AIC2
BIC2

#================================================#
#===MUHN Estimates, SEs and Criteria=============#
#================================================#

mMUHN <- function(y,sigma){
  llike <- log(2)-sigma-2*log(y)+dnorm((1-y)/(exp(sigma)*y),log = TRUE) 
  l     <- -sum(llike)
  return(l)
}

mMUHN2 <- function(y,sigma){
  llike <- log(2)-log(sigma)-2*log(y)+dnorm((1-y)/(sigma*y),log = TRUE) 
  l     <- -sum(llike)
  return(l)
}

M1         <- optim(c(0),y = y,fn = mMUHN, method = "L-BFGS-B",lower = c(-Inf), upper = c(Inf), hessian = TRUE)
M2         <- optim(exp(M1$par[1]),y = y,fn = mMUHN2, method = "L-BFGS-B",lower = c(0), upper = c(Inf), hessian = TRUE)
AIC3       <- 2*M2$value+2*1
BIC3       <- 2*M2$value+1*log(length(y))
estimates3 <- M2$par
se3        <- sqrt(diag(solve(M2$hessian)))

#=============#
#===RESULTS===#
#=============#
estimates3
se3
AIC3
BIC3

#================================================#
#===Beta Estimates, SEs and Criteria=============#
#================================================#

mB      <- function(y,theta){
  alpha <- theta[1] ; beta <- theta[2]
  llike <- dbeta(y,shape1 = exp(alpha), shape2 = exp(beta), log = TRUE)
  l     <- -sum(llike)
  return(l)
}

mB2 <- function(y,theta){
  alpha <- theta[1] ; beta <- theta[2]
  llike    <- dbeta(y,shape1 = alpha, shape2 = beta, log = TRUE)
  l        <- -sum(llike)
  return(l)
}

B1         <- optim(c(0,0),y = y,fn = mB, method = "L-BFGS-B",lower = c(-Inf,-Inf), upper = c(Inf,Inf), hessian = TRUE)
B2         <- optim(c(exp(B1$par[1]),exp(B1$par[2])),y = y,fn = mB2, method = "L-BFGS-B",lower = c(0,0), upper = c(Inf,Inf), hessian = TRUE)
AIC4       <- 2*B2$value+2*2
BIC4       <- 2*B2$value+2*log(length(y))
estimates4 <- B2$par
se4        <- sqrt(diag(solve(B2$hessian)))

#=============#
#===RESULTS===#
#=============#
estimates4
se4
AIC4
BIC4

#================================================#
#===KM Estimates, SEs and Criteria===============#
#================================================#

mKUM    <- function(y,theta){
  alpha <- theta[1] ; beta <- theta[2]
  llike <- dkumar(y,a = exp(alpha), b = exp(beta), log = TRUE)
  l     <- -sum(llike)
  return(l)
}

mKUM2   <- function(y,theta){
  alpha <- theta[1] ; beta <- theta[2]
  llike <- dkumar(y,a = alpha, b = beta, log = TRUE)
  l     <- -sum(llike)
  return(l)
}

K1         <- optim(c(1,1),y = y,fn = mKUM, method = "L-BFGS-B",lower = c(-Inf,-Inf), upper = c(Inf,Inf), hessian = TRUE)
K2         <- optim(c(exp(K1$par[1]),exp(K1$par[2])),y = y,fn = mKUM2, method = "L-BFGS-B",lower = c(0,0), upper = c(Inf,Inf), hessian = TRUE)
AIC5       <- 2*K2$value+2*2
BIC5       <- 2*K2$value+2*log(length(y))
estimates5 <- K2$par
se5        <- sqrt(diag(solve(K2$hessian)))

#=============#
#===RESULTS===#
#=============#
estimates5
se5
AIC5
BIC5










