rm(list=ls())
library(extraDistr)
library(latex2exp)
library(nortest)
library(moments)
library(LambertW)

#============================#
#========NQRs================#
#============================#

#============================#
#========KM==================#
#============================#
r1 <- qnorm(pkumar(y,2.719,44.661))  
shapiro.test(r1) 
cvm.test(r1)      
ad.test(r1)       

#============================#
#========Beta================#
#============================#
r2 <- qnorm(pbeta(y,shape1 = 5.942,shape2 = 21.206)) 
shapiro.test(r2)  
cvm.test(r2)   
ad.test(r2)   

#============================#
#========MUHN================#
#============================#
sigma <- 4.564618
r3    <- qnorm(2*pnorm((y-1)/(sigma*y),0,1)) 
shapiro.test(r3)  
cvm.test(r3)   
ad.test(r3)   

#============================#
#========UHN=================#
#============================#
sigma2 <- 0.3374166
r4     <- qnorm(2*pnorm(y/(sigma2*(1-y)))-1) 
shapiro.test(r4)  
cvm.test(r4)   
ad.test(r4)   

#============================#
#========UL==================#
#============================#
sigma3 <- 4.049086
r5     <- qnorm(1-(1-((sigma3*y)/(1+sigma3)*(y-1)))*exp(-(sigma3*y)/(1-y))) 
shapiro.test(r5)  
cvm.test(r5)   
ad.test(r5)   

