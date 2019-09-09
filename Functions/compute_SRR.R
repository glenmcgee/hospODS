####################################################################
#                                                                  #
#                          Compute SRRs                            #
#                                                                  #
####################################################################
#                             March 17                            #
####################################################################
## Compute SRRs for a given model fit
library(tidyverse)
library(glmmML)
library(gaussquad)

# expit
expit <- function(y){
  return(exp(y)/(1+exp(y)))
}

compute_SRR <- function(beta,sigma,bk,sample,frmla,testid,offset=0,Mki=NA){
  
  ## set cluster id at observation level 
  clstr <- sample[,testid]

  # must reassign name of frmla in order for it to be taken into other functions
  form <- formula(frmla)

  ## gauss hermite weights + abscissas
  nodes <- 25
  h <- as.matrix(hermite.h.quadrature.rules(nodes)[[nodes]])[,1]
  w <- as.matrix(hermite.h.quadrature.rules(nodes)[[nodes]])[,2]
  
  ## get design matrix
  X <- model.matrix(object=form,data=sample)
  
  ## define terms
  #beta <- betasig[1:length(betasig)-1]
  #sigma <- exp(betasig[length(betasig)])
  Xkibeta <- X%*%beta+offset
  
  if(length(Mki)==1){
    Mki <- rep(1,length(Xkibeta))
  }
  
  ## loop over OBSERVATIONS
  num_ki <- c(expit(Xkibeta+bk[clstr]))
  #denom_ki <- c(expit(Xkibeta+0)) old way with bk=0
  denom_ki <- rep(NA,length(num_ki))
  for (ki in (1:length(num_ki))){
    
    marg_mean <- 0
    for (q in 1:nodes){
      marg_mean <- marg_mean+(1/sqrt(pi))*w[q]*expit(Xkibeta[ki]+sqrt(2)*sigma*h[q])
    }
    
    denom_ki[ki] <- marg_mean
  }
  
  
  
  
  df <- as_tibble(cbind(denom_ki,num_ki,Mki,clstr))
  df <- df %>% group_by(clstr) %>% summarise(num=sum(num_ki*Mki),
                                             denom=sum(denom_ki*Mki),
                                             SRR=num/denom)
  
  return(c(df$SRR))
  
  
}


