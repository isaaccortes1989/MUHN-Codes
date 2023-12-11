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

#=================================#
#=========HISTOGRAM===============#
#=================================#
x      <- seq(0,1,0.001)

#==============#
#====UL========#
#==============#

sigma1 <- 4.049
f1     <- (sigma1^2/((1+sigma1)*(1-x)^3))*exp(-((sigma1*x)/(1-x))) 

#==============#
#====UHN=======#
#==============#
sigma2 <- 0.337 
f2     <- (2/(sigma2*(1-x)^2))*dnorm(x/(sigma2*(1-x)))

#==============#
#====MUHN======#
#==============#
sigma3 <- 4.565
f3     <- (2/(sigma3*x^2))*dnorm((1-x)/(sigma3*x))

#==============#
#====Beta======#
#==============#
f4   <- dbeta(x,shape1 = 5.942,shape2 = 21.206)

#==============#
#====KM========#
#==============#
f5   <- dkumar(x,2.719,44.661)

par(mar = c(4,5,3,3),las = 0,mgp = c(3.5,0.5,0))
hist(y,breaks=5,main = "",col = "white",freq = F, xlim = c(0,0.6), ylim = c(0,5), xlab = "",ylab = "",cex.axis = 1.2,cex.lab = 2,font.axis = 2,cex = 1.6)
box()

lines(x,f1,type ="l",  lwd = 2, main="", col = "red")
lines(x,f2,type ="l",  lwd = 2, main="", col = "blue")
lines(x,f3,type ="l",  lwd = 2, main="", col = "green")
lines(x,f4,type ="l",  lwd = 2, main="", col = "black")
lines(x,f5,type ="l",  lwd = 2, main="", col = "orange")


legend(0.4,4,c(as.expression(substitute("UL")),
               as.expression(substitute("UHN")),
               as.expression(substitute("MUHN")),
               as.expression(substitute("Beta")),
               as.expression(substitute("KM"))),
       col = c("red","blue","green","black","orange"), lwd = 2)
mtext(text = "z", side = 1, line = 2, cex = 1.2)
mtext(text = "Density", side = 2, line = 2, cex = 1.8)
