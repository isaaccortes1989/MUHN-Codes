sim <- function(n,sigma){
  Y <- runif(n)
  Z <- (1-sigma*qnorm(Y/2))^(-1)
  return(Z)
}

muestras = c(50,80,100)
sigma = 2.5  #Here the parameter is modified

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
  Z   <- sim(n,sigma)

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

      },error = function(e) flag = 0) 
    } 
  }
  write.table(estimados, paste("C:/Users/Isaac/Documents/Paulina/Sigma25/ESTIMADOS",n,".txt"), sep="\t",row.names=FALSE)
  write.table(see,paste("C:/Users/Isaac/Documents/Paulina/Sigma25/SE",n,".txt"),sep="\t",row.names=FALSE)
}



