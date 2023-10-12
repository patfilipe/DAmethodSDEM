##----------------------------------------------------------##----------------------------------------------------------##7
##  Inference on stochastic differential the Stochastic Gompertz model:
## 1) The fixed effects closed formula
## 2) The delta approximation method for random alpha case
## 3) The delta approximation method for random beta case 
## 4) The delta approximation method both alpha and beta random                                                                                       ##
#----------------------------------------------------------#----------------------------------------------------------#
#                                            Data                                          #
#----------------------------------------------------------#----------------------------------------------------------#


## Usage
# DA.fit(Mat_Time, Mat_Obs,  random, c(a,t,b,o,s))

## Mat_Time - matrix of m observation times. Each column represents the observation time of an individual.  
## Mat_Obs - matrix of m trajectories. Each column represents  the observations of an individual. 
## Notes: Mat_Time and Mat_Obs must  have the same dimensions. If different animals had different observation times, the remain values should have NA 

## random  - the random effects in the drift. If random=0, no random effects. If random=1,  a random effect on alpha
##                                         If random=11,  a random effect on beta. If random=2,  a random effect on both alpha and beta
## starting parameter values for the minimization
## c(a,t,b,o,s): initial values for a - alpha, t - theta (give the value 0 if no random effect is used)
##                                  b - beta, o - omega (give the value 0 if no random effect is used)
##                                  s - sigma (the diffusion coefficient assumed fixed
## If random=0 use c(a,b,s). if random=1 use c(a,t,b,s). if random=11 use c(a,b,o,s). if random=2 use c(a,t,b,o,s).    


DA.fit <- function(I,P,rnd, y) {
  library("numDeriv")
  
  d <- dim(I)  # Matrix dimension 'D';
  n <- d[2]  # Number of animals;
  m <- d[1]    # Number of weightings per animal;
  
  Nc<-vector(length=n)
  for (j in 1:n){
    for (i in 1:m){
      ifelse(I[i,j] == 'NA',Nc,Nc[j]<-Nc[j]+1)}}
  N<-Nc-1
  
  
  #----------------------------------------------------------#----------------------------------------------------------#
  #                                           Fixed effects case                                        #
  #----------------------------------------------------------#----------------------------------------------------------#
  if (rnd==0) {
    
    if (length(y)!=3) {
      print(paste("Error: wrong number of parameter"))
    }
    else {
      # - Likelihood function
      L<-function(I1,I2,P1,P2,N,x){
        -(-(N/2)*log(2*pi)-(N/2)*log(x[3]^2/(2*x[2]))-
            (1/2)*(sum(log(1-(exp(-x[2]*(I2-I1)))^2)))-
            (x[2]/(x[3]^2))*(sum((P2-x[1]-(P1-x[1])*(exp(-x[2]*(I2-I1))))^2/(1-(exp(-x[2]*(I2-I1)))^2))))
      }
      
      LTG<-function(x){
        LTG<-0
        for (i in 1:n) {
          LTG<-LTG+L(I[1:(Nc[i]-1),i],I[2:Nc[i],i],log(P[1:(Nc[i]-1),i]),log(P[2:Nc[i],i]),N[i],c(x[1],x[2],x[3]))
        }
        return(LTG)
      }
      
      #LTG is the  symmetric of the likelihood function and nlm (or nlminb) minimizes symmetric of the likelihood function. 
      mLTG<-nlm(LTG,c(y[1],y[2],y[3]),hessian=TRUE)
      uG <- mLTG$estimate[1]
      UG <- exp(uG)
      lG <- mLTG$estimate[2]
      sG <- abs(mLTG$estimate[3])
      
      # - Empirical Fisher information matrix; Inverse of the Hessian LT at point (a,b1,b2,s)
      vhG<-solve(mLTG$hessian)
      MEuG <- sqrt(vhG[1,1])
      MEUG <- sqrt(vhG[1,1])*exp(uG)
      MElG <- sqrt(vhG[2,2])
      MEsG <- sqrt(vhG[3,3])
      
      print(paste("SDE with fixed effects estimates"))
      print(paste("alpha", round(uG,5)))
      print(paste("beta", round(lG,5)))
      print(paste("sigma", round(sG,5)))
      print(paste("____"))
      print(paste("Standard deviation estimate from the empirical Fisher information matrix"))
      print(paste("alpha", round(MEuG,5)))
      print(paste("beta", round(MElG,5)))
      print(paste("sigma",round(MEsG,5)))

    }
    
    
  }
  
  #----------------------------------------------------------#----------------------------------------------------------#
  #                                            DA method for random alpha case                                          #
  #----------------------------------------------------------#----------------------------------------------------------#
  
  else if (rnd==1) {
    
    if (length(y)!=4) {
      print(paste("Error: wrong number of parameter"))
    }
    else {
      # - Likelihood function
      LG <- function(I1,I2,P1,P2,N,x){
        f <- (-1/2)*sum((P2 - x[1] - (P1 - x[1])*(exp(- x[3]*(I2 - I1))))^2/((x[4]^2/(2*x[3]))*(1 - (exp(- x[3]*(I2 - I1)))^2)))
        f1 <- (sum((1 - (exp(- x[3]*(I2 - I1))))*(P2 - x[1] - (P1 - x[1])*(exp(- x[3]*(I2 - I1))))/((x[4]^2/(2*x[3]))*(1 - (exp(- x[3]*(I2 - I1)))^2))))
        f2 <- (- sum((1 - (exp(- x[3]*(I2 - I1))))^2/((x[4]^2/(2*x[3]))*(1 - (exp(- x[3]*(I2 - I1)))^2))))
        MV <- -((-1/2)*sum(log(2*pi*((x[4]^2/(2*x[3]))*(1 - (exp(- x[3]*(I2 - I1)))^2)))) + f + log(1 + (1/2)*(x[2]^2*(f2 + (f1)^2))))
        MV}
      
      
      LTG <- function(x){
        LTG <- 0
        for (k in 1:n) {
          LTG<-LTG+LG(I[1:(Nc[k]-1),k],I[2:Nc[k],k],log(P[1:(Nc[k]-1),k]),log(P[2:Nc[k],k]),N[k],c(x[1],x[2],x[3],x[4]))
        }
        return(LTG)
      }
      #LTG is the  symmetric of the likelihood function and nlm (or nlminb) minimizes symmetric of the likelihood function.
      mLTG<-nlm(LTG,c(y[1],y[2],y[3], y[4]),hessian=TRUE)

      uG <- mLTG$estimate[1]
      UG <- exp(uG)
      tG <- mLTG$estimate[2]
      lG <- mLTG$estimate[3]
      sG <- abs(mLTG$estimate[4])

      # - Empirical Fisher information matrix; Inverse of the Hessian LT at point (a,b1,b2,s)
      vhG<-solve(mLTG$hessian)
      MEuG <- sqrt(vhG[1,1])
      MEUG <- sqrt(vhG[1,1])*exp(uG)
      MEtG <- sqrt(vhG[2,2])
      MElG <- sqrt(vhG[3,3])
      MEsG <- sqrt(vhG[4,4])

      print(paste("SDE with random effect on alpha estimates"))
      print(paste("mu", round(uG,5), "theta", round(tG,5)))
      print(paste("beta", round(lG,5)))
      print(paste("sigma", round(sG,5)))
      print(paste("___"))
      print(paste("Standard deviation estimate from the empirical Fisher information matrix"))
      print(paste("mu", round(MEuG,5), "theta", round(MEtG,5)))
      print(paste("beta", round(MElG,5)))
      print(paste("sigma",round(MEsG,5)))
    }
    
  }
  
  #----------------------------------------------------------#----------------------------------------------------------#
  #                                            DA method for random beta case                                           #
  #----------------------------------------------------------#----------------------------------------------------------#
  
  else if (rnd==11) {
    if (length(y)!=4) {
      print(paste("Error: wrong number of parameter"))
    }
    else {
      # - Likelihood function
      LG<-function(I1,I2,P1,P2,N,x){
        
        f <- (-1/2)*sum((P2 - x[1] - (P1 - x[1])*(exp(- x[2]*(I2 - I1))))^2/(((x[4]^2)/(2*x[2]))*(1 - (exp(- x[2]*(I2 - I1)))^2)) + log(2*pi*(((x[4]^2)/(2*x[2]))*(1 - (exp(- x[2]*(I2 - I1)))^2))))
        
        f1 <- (sum(((I1 - I2)*(P1 - x[1])*(exp(- x[2]*(I2 - I1)))*(P2 - x[1] - (P1 - x[1])*(exp(- x[2]*(I2 - I1)))))/(((x[4]^2)/(2*x[2]))*(1 - (exp(- x[2]*(I2 - I1)))^2))) -
                 sum((P2 - x[1] - (P1 - x[1])*(exp(- x[2]*(I2 - I1))))^2/(x[4]^2*(1 - (exp(- x[2]*(I2 - I1)))^2))) +
                 2*sum(((P2 - x[1] - (P1 - x[1])*(exp(- x[2]*(I2 - I1))))^2*(I2 - I1)*(exp(- x[2]*(I2 - I1)))^2)/(((x[4]^2)/(x[2]))*(1 - (exp(- x[2]*(I2 - I1)))^2)^2)) +
                 N/(2*x[2]) - sum(((I2 - I1)*(exp(- x[2]*(I2 - I1)))^2)/(1 - (exp(- x[2]*(I2 - I1)))^2)))
        
        f2 <- (2*sum(((I2 - I1)^2*(P1 - x[1])*(exp(- x[2]*(I2 - I1)))*(P2 - x[1] - (P1 - x[1])*(exp(- x[2]*(I2 - I1)))))/(((x[4]^2)/(x[2]))*(1 - (exp(- x[2]*(I2 - I1)))^2))) -
                 2*sum(((I2 - I1)^2*(P1 - x[1])^2*(exp(- x[2]*(I2 - I1)))^2)/(((x[4]^2)/(x[2]))*(1 - (exp(- x[2]*(I2 - I1)))^2))) -
                 4*sum(((I2 - I1)*(P1 - x[1])*(exp(- x[2]*(I2 - I1)))*(P2 - x[1] - (P1 - x[1])*(exp(- x[2]*(I2 - I1)))))/(x[4]^2*(1 - (exp(- x[2]*(I2 - I1)))^2))) +
                 8*sum(((I2 - I1)^2*(P1 - x[1])*(exp(- x[2]*(I2 - I1)))^3*(P2 - x[1] - (P1 - x[1])*(exp(- x[2]*(I2 - I1)))))/(((x[4]^2)/(x[2]))*(1 - (exp(- x[2]*(I2 - I1)))^2)^2)) +
                 4*sum(((I2 - I1)*(exp(- x[2]*(I2 - I1)))^2*(P2 - x[1] - (P1 - x[1])*(exp(- x[2]*(I2 - I1))))^2)/(x[4]^2*(1 - (exp(- x[2]*(I2 - I1)))^2)^2)) -
                 4*sum(((I2 - I1)^2*(exp(- x[2]*(I2 - I1)))^2*(P2 - x[1] - (P1 - x[1])*(exp(- x[2]*(I2 - I1))))^2)/((x[4]^2/x[2])*(1 - (exp(- x[2]*(I2 - I1)))^2)^2)) -
                 8*sum(((I2 - I1)^2*(exp(- x[2]*(I2 - I1)))^4*(P2 - x[1] - (P1 - x[1])*(exp(- x[2]*(I2 - I1))))^2)/((x[4]^2/x[2])*(1 - (exp(- x[2]*(I2 - I1)))^2)^3)) +
                 2*sum(((I2 - I1)^2*(exp(- x[2]*(I2 - I1)))^4)/(1 - (exp(- x[2]*(I2 - I1)))^2)^2) +
                 2*sum(((I2 - I1)^2*(exp(- x[2]*(I2 - I1)))^2)/(1 - (exp(- x[2]*(I2 - I1)))^2)) -
                 N/(2*x[2]^2))
        
        MV<- - ( f + log(1 + (1/2)*(x[3]^2)*(f2 + (f1)^2)))
        MV}
      
      LTG <- function(x){
        LTG <- 0
        for (k in 1:n) {
          LTG <- LTG + LG(I[1:(Nc[k]-1),k],I[2:Nc[k],k],log(P[1:(Nc[k]-1),k]),log(P[2:Nc[k],k]),N[k],c(x[1],x[2],x[3],x[4]))
        }
        return(LTG)
      }
      
      #LTG is the  symmetric of the likelihood function and nlm (or nlminb) minimizes symmetric of the likelihood function.
      mLTG <- nlm(LTG,c(y[1],y[2],y[3], y[4]),hessian=TRUE)

      uG <- mLTG$estimate[1]
      UG <- exp(uG)
      lG <- mLTG$estimate[2]
      oG <- mLTG$estimate[3]
      sG <- abs(mLTG$estimate[4])

      # - Empirical Fisher information matrix; Inverse of the Hessian LT at point (a,b1,b2,s)
      vhG<-solve(mLTG$hessian)

      MEuG <- sqrt(vhG[1,1])
      MEUG <- sqrt(vhG[1,1])*exp(uG)
      MElG <- sqrt(vhG[2,2])
      MEoG <- sqrt(vhG[3,3])
      MEsG <- sqrt(vhG[4,4])


      print(paste("SDE with random effect on beta estimates"))
      print(paste("alpha", round(uG,5) ))
      print(paste("lambda", round(lG,5), "omega", round(oG,5)))
      print(paste("sigma", round(sG,5)))
      print(paste("___"))
      print(paste("Standard deviation estimate from the empirical Fisher information matrix"))
      print(paste("alpha", round(MEuG,5)))
      print(paste("lambda", round(MElG,5), "omega", round(MEoG,5)))
      print(paste("sigma",round(MEsG,5)))

    }
    
    
  }
  
  
  #----------------------------------------------------------#----------------------------------------------------------#
  #                                               DA method for both alpha and beta random                              #                       
  #----------------------------------------------------------#----------------------------------------------------------#
  
  else if (rnd==2) {
    
    if (length(y)!=5) {
      print(paste("Error: wrong number of parameter"))
    }
    else {
      # - Likelihood function
      LG <- function(I1,I2,P1,P2,N,x){
        f <- ((-x[3]/x[5]^2)*sum((P2-x[1]-(P1-x[1])*(exp(-x[3]*(I2-I1))))^2/(1-(exp(-x[3]*(I2-I1)))^2))-(1/2)*sum(log(((pi*x[5]^2)/x[3])*(1-(exp(-x[3]*(I2-I1)))^2))))
        fa1 <- (((-2*x[3])/x[5]^2)*sum(((P2-x[1]-(P1-x[1])*(exp(-x[3]*(I2-I1))))*((exp(-x[3]*(I2-I1)))-1))/(1-(exp(-x[3]*(I2-I1)))^2)))
        fa2 <- (((-2*x[3])/x[5]^2)*sum(((exp(-x[3]*(I2-I1)))-1)^2/(1-(exp(-x[3]*(I2-I1)))^2)))
        fb1 <- (((-2*x[3])/x[5]^2)*sum(((I2-I1)*(P1-x[1])*(exp(-x[3]*(I2-I1)))*(P2-x[1]-(P1-x[1])*(exp(-x[3]*(I2-I1)))))/(1-(exp(-x[3]*(I2-I1)))^2)) -
                  (1/x[5]^2)*sum((P2-x[1]-(P1-x[1])*(exp(-x[3]*(I2-I1))))^2/(1-(exp(-x[3]*(I2-I1)))^2)) +
                  ((2*x[3])/x[5]^2)*sum(((I2-I1)*(exp(-x[3]*(I2-I1)))^2*(P2-x[1]-(P1-x[1])*(exp(-x[3]*(I2-I1))))^2)/(1-(exp(-x[3]*(I2-I1)))^2)^2) -
                  sum(((I2-I1)*(exp(-x[3]*(I2-I1)))^2)/(1-(exp(-x[3]*(I2-I1)))^2)) + N/(2*x[3]))
        fb2 <- (((2*x[3])/x[5]^2)*sum(((I2-I1)^2*(P1-x[1])*(exp(-x[3]*(I2-I1)))*(P2-x[1]-(P1-x[1])*(exp(-x[3]*(I2-I1)))))/(1-(exp(-x[3]*(I2-I1)))^2)) -
                  ((2*x[3])/x[5]^2)*sum(((I2-I1)^2*(P1-x[1])^2*(exp(-x[3]*(I2-I1)))^2)/(1-(exp(-x[3]*(I2-I1)))^2)) -
                  (4/x[5]^2)*sum(((I2-I1)*(P1-x[1])*(exp(-x[3]*(I2-I1)))*(P2-x[1]-(P1-x[1])*(exp(-x[3]*(I2-I1)))))/(1-(exp(-x[3]*(I2-I1)))^2)) +
                  ((8*x[3])/x[5]^2)*sum(((I2-I1)^2*(P1-x[1])*(exp(-x[3]*(I2-I1)))^3*(P2-x[1]-(P1-x[1])*(exp(-x[3]*(I2-I1)))))/(1-(exp(-x[3]*(I2-I1)))^2)^2) +
                  (4/x[5]^2)*sum(((I2-I1)*(exp(-x[3]*(I2-I1)))^2*(P2-x[1]-(P1-x[1])*(exp(-x[3]*(I2-I1))))^2)/(1-(exp(-x[3]*(I2-I1)))^2)^2) -
                  ((4*x[3])/x[5]^2)*sum(((I2-I1)^2*(exp(-x[3]*(I2-I1)))^2*(P2-x[1]-(P1-x[1])*(exp(-x[3]*(I2-I1))))^2)/(1-(exp(-x[3]*(I2-I1)))^2)^2) -
                  ((8*x[3])/x[5]^2)*sum(((I2-I1)^2*(exp(-x[3]*(I2-I1)))^4*(P2-x[1]-(P1-x[1])*(exp(-x[3]*(I2-I1))))^2)/(1-(exp(-x[3]*(I2-I1)))^2)^3) +
                  2*sum(((I2-I1)^2*(exp(-x[3]*(I2-I1)))^4)/(1-(exp(-x[3]*(I2-I1)))^2)^2) +
                  2*sum(((I2-I1)^2*(exp(-x[3]*(I2-I1)))^2)/(1-(exp(-x[3]*(I2-I1)))^2)) - N/(2*x[3]^2))
        MV <- -(-sum(P2)+f+log(1+(x[2]^2/2)*(fa2+(fa1)^2)+(x[4]^2/2)*(fb2+(fb1)^2)))
        MV}
      # - Sum of likelihood functions
      LTG<-function(x){
        LTG<-0
        for (k in 1:n) {
          LTG <- LTG + LG(I[1:(Nc[k]-1),k],I[2:Nc[k],k],log(P[1:(Nc[k]-1),k]),log(P[2:Nc[k],k]),N[k],c(x[1],x[2],x[3],x[4],x[5]))
        }
        return(LTG)
      }
      
      #LTG is the  symmetric of the likelihood function and nlm (or nlminb) minimizes symmetric of the likelihood function.
      mLTG <- nlm(LTG,c(y[1],y[2],y[3], y[4],y[5]),hessian=TRUE)
      uG <- mLTG$estimate[1]
      UG <- exp(uG)
      tG <- mLTG$estimate[2]
      lG <- mLTG$estimate[3]
      oG <- mLTG$estimate[4]
      sG <- abs(mLTG$estimate[5])

      # - Empirical Fisher information matrix; Inverse of the Hessian LT at point (a,b1,b2,s)
      vhG<-solve(mLTG$hessian)
      MEuG <- sqrt(vhG[1,1])
      MEUG <- sqrt(vhG[1,1])*exp(uG)
      MEtG <- sqrt(vhG[2,2])
      MElG <- sqrt(vhG[3,3])
      MEoG <- sqrt(vhG[4,4])
      MEsG <- sqrt(vhG[5,5])

      print(paste("SDE with random effect on alpha and beta estimates"))
      print(paste("mu", round(uG,5), "theta", round(tG,5)))
      print(paste("lambda", round(lG,5), "omega", round(oG,5)))
      print(paste("sigma", round(sG,5)))
      print(paste("___"))
      print(paste("Standard deviation estimate from the empirical Fisher information matrix"))
      print(paste("mu", round(MEuG,5), "theta", round(MEtG,5)))
      print(paste("lambda", round(MElG,5), "omega", round(MEoG,5)))
      print(paste("sigma",round(MEsG,5)))

    }
    
  }
  
  else {
    print("Choice of parameter rnd incorrect")
  }
  
}
