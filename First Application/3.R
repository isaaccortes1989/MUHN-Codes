rm(list=ls())
library(extraDistr)
library(latex2exp)
library(nortest)
library(moments)
library(LambertW)

#============================#
#========NQRs================#
#============================#

y <- c(0.0903296,  0.2036540, 0.2043140, 0.2808870, 0.1976530, 0.3286410,
       0.1486220,  0.1623940, 0.2627270, 0.1794550, 0.3266350, 0.2300810,
       0.1833120,  0.1509440, 0.2000710, 0.1918020, 0.1541920, 0.4641250,
       0.1170630,  0.1481410, 0.1448100, 0.1330830, 0.2760160, 0.4204770, 
       0.1224170,  0.2285950, 0.1138520, 0.2252140, 0.1769690, 0.2007440,
       0.1670450,  0.2316230, 0.2910290, 0.3412730, 0.4387120, 0.2626510, 
       0.1896510,  0.1725670, 0.2400770, 0.3116460, 0.1635860, 0.1824530,
       0.1641270,  0.1534810, 0.1618650, 0.2760160, 0.2538320, 0.2004470)

#===========================================================#
#===Run the code lines for each distribution independently==#
#===========================================================#

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

