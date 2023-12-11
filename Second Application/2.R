rm(list=ls())
library(extraDistr)
library(latex2exp)
library(nortest)
library(moments)
library(LambertW)


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

#=================================#
#=========HISTOGRAM===============#
#=================================#
x      <- seq(0,1,0.001)

#==============#
#====UL========#
#==============#

sigma1 <- 4.140816                 
f1     <- (sigma1^2/((1+sigma1)*(1-x)^3))*exp(-((sigma1*x)/(1-x)))

#==============#
#====UHN=======#
#==============#

sigma2 <- 0.4425385                  
f2     <- (2/(sigma2*(1-x)^2))*dnorm(x/(sigma2*(1-x)))

#==============#
#====MUHN======#
#==============#

sigma3 <- 7.230709                  
f3     <- (2/(sigma3*x^2))*dnorm((1-x)/(sigma3*x)) 

#==============#
#====MUHN======#
#==============#
f4   <- dbeta(x,shape1 = 2.325258,shape2 = 9.437853)  

#==============#
#====KM========#
#==============#
f5   <- dkumar(x,1.61014,10.51161)   




par(mar = c(4,5,3,3),las = 0,mgp = c(3.5,0.5,0))
hist(y,breaks=13,main = "",col = "white",freq = F, xlim = c(0,1), ylim = c(0,6), xlab = "",ylab = "",cex.axis = 1.2,cex.lab = 2,font.axis = 2,cex = 1.6)
box()
lines(x,f1,type ="l",  lwd = 2, main="", col = "red")
lines(x,f2,type ="l",  lwd = 2, main="", col = "blue")
lines(x,f3,type ="l",  lwd = 2, main="", col = "green")
lines(x,f4,type ="l",  lwd = 2, main="", col = "black")
lines(x,f5,type ="l",  lwd = 2, main="", col = "orange")


legend(0.6,5,c(as.expression(substitute("UL")),
               as.expression(substitute("UHN")),
               as.expression(substitute("MUHN")),
               as.expression(substitute("Beta")),
               as.expression(substitute("KM"))),
       col = c("red","blue","green","black","orange"), lwd = 2)
mtext(text = "z", side = 1, line = 2, cex = 1.2)
mtext(text = "Density", side = 2, line = 2, cex = 1.8)
