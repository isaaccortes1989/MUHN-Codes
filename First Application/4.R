rm(list=ls())
library(extraDistr)
library(latex2exp)
library(nortest)
library(moments)
library(LambertW)

#===========================================================#
#===Run the code lines for each distribution independently==#
#===========================================================#

#----------------------#
#--QQ-PLOT MUHN--------# 
#----------------------#

sim <- function(n,sigma){ #Random Numbers
  Y <- runif(n)
  Z <- (1-sigma*qnorm(Y/2))^(-1)
  return(Z)
}

set.seed(231)  

y <- c(0.0903296,  0.2036540, 0.2043140, 0.2808870, 0.1976530, 0.3286410,
       0.1486220,  0.1623940, 0.2627270, 0.1794550, 0.3266350, 0.2300810,
       0.1833120,  0.1509440, 0.2000710, 0.1918020, 0.1541920, 0.4641250,
       0.1170630,  0.1481410, 0.1448100, 0.1330830, 0.2760160, 0.4204770, 
       0.1224170,  0.2285950, 0.1138520, 0.2252140, 0.1769690, 0.2007440,
       0.1670450,  0.2316230, 0.2910290, 0.3412730, 0.4387120, 0.2626510, 
       0.1896510,  0.1725670, 0.2400770, 0.3116460, 0.1635860, 0.1824530,
       0.1641270,  0.1534810, 0.1618650, 0.2760160, 0.2538320, 0.2004470)

n      <- 48 
sigma  <- 4.565 
Y      <- sim(n,sigma)

R      <- qnorm(2*pnorm((Y-1)/(sigma*Y)))
r      <- qnorm(2*pnorm((y-1)/(sigma*y)))


par(mar = c(4,5,3,3),las = 0,mgp = c(3.5,0.5,0))
qqplot(R,r,main="MUHN", plot.it = TRUE, xlab = "", ylab = "",xlim = c(-4,4),ylim=c(-4,4),
       pch = 19,cex.axis = 1.2,cex.lab = 2,font.axis = 2,cex = 1.6,cex.main=2.5)
curve(1*x,add = TRUE,col="black",lwd = 2)
mtext(text = "Theoretical Quantiles", side = 1, line = 2, cex = 2)
mtext(text = "Sample Quantiles", side = 2, line = 2, cex = 2)
text(2.19,-3.7, expression("AD:0.902")  , cex = 1.6)
text(1.97,-3.2, expression("CVM:0.904"), cex = 1.6)
text(2.09,-2.7, expression("SW:0.853") , cex = 1.6)



#----------------------#
#--QQ-PLOT BETA--------# 
#----------------------#

set.seed(231)  

y <- c(0.0903296,  0.2036540, 0.2043140, 0.2808870, 0.1976530, 0.3286410,
       0.1486220,  0.1623940, 0.2627270, 0.1794550, 0.3266350, 0.2300810,
       0.1833120,  0.1509440, 0.2000710, 0.1918020, 0.1541920, 0.4641250,
       0.1170630,  0.1481410, 0.1448100, 0.1330830, 0.2760160, 0.4204770, 
       0.1224170,  0.2285950, 0.1138520, 0.2252140, 0.1769690, 0.2007440,
       0.1670450,  0.2316230, 0.2910290, 0.3412730, 0.4387120, 0.2626510, 
       0.1896510,  0.1725670, 0.2400770, 0.3116460, 0.1635860, 0.1824530,
       0.1641270,  0.1534810, 0.1618650, 0.2760160, 0.2538320, 0.2004470)

n      <- 48 
Y      <- rbeta(n,shape1=5.942,shape2=21.206)

R      <- qnorm(pbeta(Y,shape1=5.942,shape2=21.206))
r      <- qnorm(pbeta(y,shape1=5.942,shape2=21.206))

par(mar = c(4,5,3,3),las = 0,mgp = c(3.5,0.5,0))
qqplot(R,r,main="Beta", plot.it = TRUE, xlab = "", ylab = "",xlim = c(-4,4),ylim=c(-4,4),
       pch = 19,cex.axis = 1.2,cex.lab = 2,font.axis = 2,cex = 1.6,cex.main=2.5)
curve(1*x,add = TRUE,col="black",lwd = 2)
mtext(text = "Theoretical Quantiles", side = 1, line = 2, cex = 2)
mtext(text = "Sample Quantiles", side = 2, line = 2, cex = 2)
text(2.19,-3.7, expression("AD:0.043"),cex = 1.6)
text(1.97,-3.2, expression("CVM:0.047"),cex = 1.6)
text(2.09,-2.7, expression("SW:0.044"), cex = 1.6)

#----------------------#
#--QQ-PLOT KUMAR-------#
#----------------------#

set.seed(231)  

y <- c(0.0903296,  0.2036540, 0.2043140, 0.2808870, 0.1976530, 0.3286410,
       0.1486220,  0.1623940, 0.2627270, 0.1794550, 0.3266350, 0.2300810,
       0.1833120,  0.1509440, 0.2000710, 0.1918020, 0.1541920, 0.4641250,
       0.1170630,  0.1481410, 0.1448100, 0.1330830, 0.2760160, 0.4204770, 
       0.1224170,  0.2285950, 0.1138520, 0.2252140, 0.1769690, 0.2007440,
       0.1670450,  0.2316230, 0.2910290, 0.3412730, 0.4387120, 0.2626510, 
       0.1896510,  0.1725670, 0.2400770, 0.3116460, 0.1635860, 0.1824530,
       0.1641270,  0.1534810, 0.1618650, 0.2760160, 0.2538320, 0.2004470)

n      <- 48 
Y      <- rkumar(n,2.719,44.661)

R      <- qnorm(pkumar(Y,2.719,44.661))
r      <- qnorm(pkumar(y,2.719,44.661))

par(mar = c(4,5,3,3),las = 0,mgp = c(3.5,0.5,0))
qqplot(R,r,main = "KM", plot.it = TRUE, xlab = "", ylab = "",xlim = c(-4,4),ylim=c(-4,4),
       pch = 19,cex.axis = 1.2,cex.lab = 2,font.axis = 2,cex = 1.6,cex.main=2.5)
curve(1*x,add = TRUE,col="black",lwd = 2)
mtext(text = "Theoretical Quantiles", side = 1, line = 2, cex = 2)
mtext(text = "Sample Quantiles", side = 2, line = 2, cex = 2)
text(2.19,-3.7, expression("AD:0.003"),cex = 1.6)
text(1.97,-3.2, expression("CVM:0.004"),cex = 1.6)
text(2.09,-2.7, expression("SW:0.002"), cex = 1.6)

#----------------------#
#--QQ-PLOT UHN---------#
#----------------------#

sim <- function(n,sigma){ #Random Numbers
  Y <- runif(n)
  Z <- (sigma*qnorm((Y+1)/2))/(1+sigma*qnorm((Y+1)/2))  
  return(Z)
}

y <- c(0.0903296,0.2036540,0.2043140,0.2808870,0.1976530,0.3286410,0.1486220,0.1623940,
       0.2627270,0.1794550,0.3266350,0.2300810,0.1833120,0.1509440,0.2000710,0.1918020,
       0.1541920,0.4641250,0.1170630,0.1481410,0.1448100,0.1330830,0.2760160,0.4204770,
       0.1224170,0.2285950,0.1138520,0.2252140,0.1769690,0.2007440,0.1670450,0.2316230,
       0.2910290,0.3412730,0.4387120,0.2626510,0.1896510,0.1725670,0.2400770,0.3116460,
       0.1635860,0.1824530,0.1641270,0.1534810,0.1618650,0.2760160,0.2538320,0.2004470)

set.seed(231)  
n      <- 48 
sigma  <- 0.337 
Y      <- sim(n,sigma)

R      <- qnorm(2*pnorm(Y/(sigma*(1-Y)))-1)
r      <- qnorm(2*pnorm(y/(sigma*(1-y)))-1)


par(mar = c(4,5,3,3),las = 0,mgp = c(3.5,0.5,0))
qqplot(R,r,main="UHN", plot.it = TRUE, xlab = "", ylab = "",xlim = c(-4,4),ylim=c(-4,4),
       pch = 19,cex.axis = 1.2,cex.lab = 2,font.axis = 2,cex = 1.6,cex.main=2.5)
curve(1*x,add = TRUE,col="black",lwd = 2)
mtext(text = "Theoretical Quantiles", side = 1, line = 2, cex = 2)
mtext(text = "Sample Quantiles", side = 2, line = 2, cex = 2)
text(2.19,-3.7, expression("AD:<0.001")  , cex = 1.6)
text(1.97,-3.2, expression("CVM:<0.001"), cex = 1.6)
text(2.09,-2.7, expression("SW:<0.001") , cex = 1.6)

#----------------------#
#--QQ-PLOT UL----------#
#----------------------#

sim <- function(n,sigma){ #Random Numbers
  U <- runif(n)
  Z <- (1+sigma)*(U-1)*exp(-(1+sigma))
  Y <- (1+sigma+W(Z,branch = -1))/(1+W(Z, branch = -1))
  return(Y)
}

y <- c(0.0903296,0.2036540,0.2043140,0.2808870,0.1976530,0.3286410,0.1486220,0.1623940,
       0.2627270,0.1794550,0.3266350,0.2300810,0.1833120,0.1509440,0.2000710,0.1918020,
       0.1541920,0.4641250,0.1170630,0.1481410,0.1448100,0.1330830,0.2760160,0.4204770,
       0.1224170,0.2285950,0.1138520,0.2252140,0.1769690,0.2007440,0.1670450,0.2316230,
       0.2910290,0.3412730,0.4387120,0.2626510,0.1896510,0.1725670,0.2400770,0.3116460,
       0.1635860,0.1824530,0.1641270,0.1534810,0.1618650,0.2760160,0.2538320,0.2004470)


set.seed(231)  
n      <- 48 
sigma  <- 4.049 
Y      <- sim(n,sigma)

R      <- qnorm(1-(1-((sigma*Y)/(1+sigma)*(Y-1)))*exp(-(sigma*Y)/(1-Y))) 
r      <- qnorm(1-(1-((sigma*y)/(1+sigma)*(y-1)))*exp(-(sigma*y)/(1-y))) 

par(mar = c(4,5,3,3),las = 0,mgp = c(3.5,0.5,0))
qqplot(R,r,main="UL", plot.it = TRUE, xlab = "", ylab = "",xlim = c(-4,4),ylim=c(-4,4),
       pch = 19,cex.axis = 1.2,cex.lab = 2,font.axis = 2,cex = 1.6,cex.main=2.5)
curve(1*x,add = TRUE,col="black",lwd = 2)
mtext(text = "Theoretical Quantiles", side = 1, line = 2, cex = 2)
mtext(text = "Sample Quantiles", side = 2, line = 2, cex = 2)
text(2.19,-3.7, expression("AD:0.007")  , cex = 1.6)
text(1.97,-3.2, expression("CVM:0.010"), cex = 1.6)
text(2.09,-2.7, expression("SW:0.005") , cex = 1.6)

