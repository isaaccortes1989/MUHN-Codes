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

set.seed(36)  

y <- c(0.666666667, 0.25, 0.2, 0.285714286, 0.3,
       0.333333333, 0.117647059 , 0.260869565, 
       0.303030303, 0.23255814, 0.295081967, 
       0.186666667, 0.519230769, 0.223880597, 
       0.155462185, 0.304093567 , 0.211981567, 
       0.191806331, 0.150316456, 0.152815013, 
       0.190889371, 0.192644483, 0.125574273, 
       0.188819876, 0.156626506, 0.107526882, 
       0.126582278, 0.105551497, 0.096667766, 
       0.109576968, 0.089108911, 0.101898582, 
       0.069335719, 0.071443406, 0.058835027, 
       0.077533357, 0.071332887, 0.081372097)

n      <- 38 
sigma  <- 7.230709
Y      <- sim(n,sigma)

R      <- qnorm(2*pnorm((Y-1)/(sigma*Y)))
r      <- qnorm(2*pnorm((y-1)/(sigma*y)))


par(mar = c(4,5,3,3),las = 0,mgp = c(3.5,0.5,0))
qqplot(R,r,main="MUHN", plot.it = TRUE, xlab = "", ylab = "",xlim = c(-4,4),ylim=c(-4,4),
       pch = 19,cex.axis = 1.2,cex.lab = 2,font.axis = 2,cex = 1.6,cex.main=2.5)
curve(1*x,add = TRUE,col="black",lwd = 2)
mtext(text = "Theoretical Quantiles", side = 1, line = 2, cex = 2)
mtext(text = "Sample Quantiles", side = 2, line = 2, cex = 2)
text(2.19,-3.7, expression("AD:0.360")  , cex = 1.6)
text(1.97,-3.2, expression("CVM:0.295"), cex = 1.6)
text(2.09,-2.7, expression("SW:0.589") , cex = 1.6)

#----------------------#
#--QQ-PLOT BETA--------# 
#----------------------#

set.seed(36)  

y <- c(0.0903296,  0.2036540, 0.2043140, 0.2808870, 0.1976530, 0.3286410,
       0.1486220,  0.1623940, 0.2627270, 0.1794550, 0.3266350, 0.2300810,
       0.1833120,  0.1509440, 0.2000710, 0.1918020, 0.1541920, 0.4641250,
       0.1170630,  0.1481410, 0.1448100, 0.1330830, 0.2760160, 0.4204770, 
       0.1224170,  0.2285950, 0.1138520, 0.2252140, 0.1769690, 0.2007440,
       0.1670450,  0.2316230, 0.2910290, 0.3412730, 0.4387120, 0.2626510, 
       0.1896510,  0.1725670, 0.2400770, 0.3116460, 0.1635860, 0.1824530,
       0.1641270,  0.1534810, 0.1618650, 0.2760160, 0.2538320, 0.2004470)

n      <- 38 
Y      <- rbeta(n,shape1 = 2.325258,shape2 = 9.437853)

R      <- qnorm(pbeta(Y,shape1 = 2.325258,shape2 = 9.437853))
r      <- qnorm(pbeta(y,shape1 = 2.325258,shape2 = 9.437853))

par(mar = c(4,5,3,3),las = 0,mgp = c(3.5,0.5,0))
qqplot(R,r,main="Beta", plot.it = TRUE, xlab = "", ylab = "",xlim = c(-4,4),ylim=c(-4,4),
       pch = 19,cex.axis = 1.2,cex.lab = 2,font.axis = 2,cex = 1.6,cex.main=2.5)
curve(1*x,add = TRUE,col="black",lwd = 2)
mtext(text = "Theoretical Quantiles", side = 1, line = 2, cex = 2)
mtext(text = "Sample Quantiles", side = 2, line = 2, cex = 2)
text(2.19,-3.7, expression("AD:0.074"),cex = 1.6)
text(1.97,-3.2, expression("CVM:0.217"),cex = 1.6)
text(2.09,-2.7, expression("SW:0.006"), cex = 1.6)

#----------------------#
#--QQ-PLOT KUMAR-------#
#----------------------#

set.seed(36) 
 
y <- c(0.0903296,  0.2036540, 0.2043140, 0.2808870, 0.1976530, 0.3286410,
       0.1486220,  0.1623940, 0.2627270, 0.1794550, 0.3266350, 0.2300810,
       0.1833120,  0.1509440, 0.2000710, 0.1918020, 0.1541920, 0.4641250,
       0.1170630,  0.1481410, 0.1448100, 0.1330830, 0.2760160, 0.4204770, 
       0.1224170,  0.2285950, 0.1138520, 0.2252140, 0.1769690, 0.2007440,
       0.1670450,  0.2316230, 0.2910290, 0.3412730, 0.4387120, 0.2626510, 
       0.1896510,  0.1725670, 0.2400770, 0.3116460, 0.1635860, 0.1824530,
       0.1641270,  0.1534810, 0.1618650, 0.2760160, 0.2538320, 0.2004470)


n      <- 38 
Y      <- rkumar(n,1.61014,10.51161)

R      <- qnorm(pkumar(Y,1.61014,10.51161))
r      <- qnorm(pkumar(y,1.61014,10.51161))

par(mar = c(4,5,3,3),las = 0,mgp = c(3.5,0.5,0))
qqplot(R,r,main = "KM", plot.it = TRUE, xlab = "", ylab = "",xlim = c(-4,4),ylim=c(-4,4),
       pch = 19,cex.axis = 1.2,cex.lab = 2,font.axis = 2,cex = 1.6,cex.main=2.5)
curve(1*x,add = TRUE,col="black",lwd = 2)
mtext(text = "Theoretical Quantiles", side = 1, line = 2, cex = 2)
mtext(text = "Sample Quantiles", side = 2, line = 2, cex = 2)
text(2.19,-3.7, expression("AD:0.028"),cex = 1.6)
text(1.97,-3.2, expression("CVM:0.105"),cex = 1.6)
text(2.09,-2.7, expression("SW:0.002"), cex = 1.6)

#----------------------#
#--QQ-PLOT UHN---------#
#----------------------#

sim <- function(n,sigma){
  Y <- runif(n)
  Z <- (sigma*qnorm((Y+1)/2))/(1+sigma*qnorm((Y+1)/2))  
  return(Z)
}

y <- c(0.0903296,  0.2036540, 0.2043140, 0.2808870, 0.1976530, 0.3286410,
       0.1486220,  0.1623940, 0.2627270, 0.1794550, 0.3266350, 0.2300810,
       0.1833120,  0.1509440, 0.2000710, 0.1918020, 0.1541920, 0.4641250,
       0.1170630,  0.1481410, 0.1448100, 0.1330830, 0.2760160, 0.4204770, 
       0.1224170,  0.2285950, 0.1138520, 0.2252140, 0.1769690, 0.2007440,
       0.1670450,  0.2316230, 0.2910290, 0.3412730, 0.4387120, 0.2626510, 
       0.1896510,  0.1725670, 0.2400770, 0.3116460, 0.1635860, 0.1824530,
       0.1641270,  0.1534810, 0.1618650, 0.2760160, 0.2538320, 0.2004470)

set.seed(36)  
n      <- 38 
sigma  <- 0.4425385 
Y      <- sim(n,sigma)

R      <- qnorm(2*pnorm(Y/(sigma*(1-Y)))-1)
r      <- qnorm(2*pnorm(y/(sigma*(1-y)))-1)


par(mar = c(4,5,3,3),las = 0,mgp = c(3.5,0.5,0))
qqplot(R,r,main="UHN", plot.it = TRUE, xlab = "", ylab = "",xlim = c(-4.5,4.5),ylim=c(-4.5,4.5),
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
sim <- function(n,sigma){
  U <- runif(n)
  Z <- (1+sigma)*(U-1)*exp(-(1+sigma))
  Y <- (1+sigma+W(Z,branch = -1))/(1+W(Z, branch = -1))
  return(Y)
}

y <- c(0.0903296,  0.2036540, 0.2043140, 0.2808870, 0.1976530, 0.3286410,
       0.1486220,  0.1623940, 0.2627270, 0.1794550, 0.3266350, 0.2300810,
       0.1833120,  0.1509440, 0.2000710, 0.1918020, 0.1541920, 0.4641250,
       0.1170630,  0.1481410, 0.1448100, 0.1330830, 0.2760160, 0.4204770, 
       0.1224170,  0.2285950, 0.1138520, 0.2252140, 0.1769690, 0.2007440,
       0.1670450,  0.2316230, 0.2910290, 0.3412730, 0.4387120, 0.2626510, 
       0.1896510,  0.1725670, 0.2400770, 0.3116460, 0.1635860, 0.1824530,
       0.1641270,  0.1534810, 0.1618650, 0.2760160, 0.2538320, 0.2004470)

set.seed(36)  
n      <- 38 
sigma  <- 4.140816 
Y      <- sim(n,sigma)

R      <- qnorm(1-(1-((sigma*Y)/(1+sigma)*(Y-1)))*exp(-(sigma*Y)/(1-Y))) 
r      <- qnorm(1-(1-((sigma*y)/(1+sigma)*(y-1)))*exp(-(sigma*y)/(1-y))) 

par(mar = c(4,5,3,3),las = 0,mgp = c(3.5,0.5,0))
qqplot(R,r,main="UL", plot.it = TRUE, xlab = "", ylab = "",xlim = c(-4,4),ylim=c(-4,4),
       pch = 19,cex.axis = 1.2,cex.lab = 2,font.axis = 2,cex = 1.6,cex.main=2.5)
curve(1*x,add = TRUE,col="black",lwd = 2)
mtext(text = "Theoretical Quantiles", side = 1, line = 2, cex = 2)
mtext(text = "Sample Quantiles", side = 2, line = 2, cex = 2)
text(2.19,-3.7, expression("AD:0.003")  , cex = 1.6)
text(1.97,-3.2, expression("CVM:0.016"), cex = 1.6)
text(2.09,-2.7, expression("SW:<0.001") , cex = 1.6) 


