sim2 <- function(n,sigma){
  Y <- rnorm(n)
  X <- sigma*abs(Y)
  Z <- 1/(1+X)
  return(Z)
}

muestras = c(50,80,100)
sigma    = 2.5  #Here the parameter is modified

for(nss in 1:length(muestras)) {
#========================================
  set.seed(29)  #Here the seed is modified
  n <- muestras[nss]
  estimados <- c()
  see       <- c()
  replicas  <- 1000
#=========================================
  for(i in 1 :replicas){
  flag <- 0
  while(flag==0){
  Z   <- sim2(n,sigma)

  tryCatch({

ML  <- sqrt(mean((((Z-1)/Z)^2))) 
SE  <- sqrt(ML^2/(2*n)) 

     if((sum(SE=="NA")==0)&(sum(SE=="NaN")==0)&(sum(SE=="Error")==0))
        {
          estimados <- rbind(estimados,ML)			
          see       <- rbind(see,SE)
          flag      <- 1
          print(i)
        }else{flag <- 0}  

      },error = function(e) flag = 0) #trycatch
    } #while
  }#for
  write.table(estimados, paste("C:/Users/Isaac/Documents/Paulina/Sigma25/MESTIMADOS",n,".txt"), sep="\t",row.names=FALSE)
  write.table(see,paste("C:/Users/Isaac/Documents/Paulina/Sigma25/MSE",n,".txt"),sep="\t",row.names=FALSE)
}



