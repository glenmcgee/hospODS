####################################################################
#                                                                  #
#                     WGLMM -   CSCC Weights                       #
#                                                                  #
####################################################################
#                             April 2017                            #
####################################################################
### Computes components of sandwich-form estimators


## Note that we need to somehow fix the fact that we get NA/inf
## Right now were just replacing log(0) with -750
library(tidyverse)
library(glmmML)
library(gaussquad)

# expit
expit <- function(y){
  return(exp(y)/(1+exp(y)))
}


####################################################################
#                   Weighted  Log-Likelihood                       #
####################################################################
## calculate negative log likelihood (to be minimized)
## note: we drop the 1/sqrt(pi) because it is simply a constant and we are minimizing
wglmm_loglik_CSCC <- function(betasig,clstr,X,Y,M_ki){  #M0,M1){

  ## gauss hermite weights + abscissas
  nodes <- 25
  h <- as.matrix(hermite.h.quadrature.rules(nodes)[[nodes]])[,1]
  w <- as.matrix(hermite.h.quadrature.rules(nodes)[[nodes]])[,2]

  ## define terms
  beta <- betasig[1:length(betasig)-1]
  sigma <- betasig[length(betasig)] ## glmmML version maximizes wrt sigma
  

  compute_gk <- function(q)  {

    ## Xk0%*%beta
    Xkibeta <- X%*%beta

    ## compute g
    gki <-  c(exp(Y*(Xkibeta+sqrt(2)*sigma*h[q])) /(1+exp(Xkibeta+sqrt(2)*sigma*h[q])) )
    gki[!is.finite(gki)] <- 0   ## error handling

    # ## get weights at the patient level
    # M0ki <- M0[clstr]
    # M1ki <- M1[clstr]

    ## take product of cluster memebers
    df <- as_tibble(cbind(clstr,Y,gki,M_ki)) #M0ki,M1ki))
    ## apply sampling weights
    df$wg <-  gki^M_ki  #(Y*M1ki+(1-Y)*M0ki)
    ## apply GH weight, and apply IPW weights for case or control
    df <- df %>% group_by(clstr) %>% summarise(g_k=w[q]*prod(wg))



    ## apply sampling weights and report cluster specific vals for the qth node
    return(c(df$g_k))


  }

  ## loop over nodes
  denom <- rep(0,max(clstr))
  for (q in 1:nodes){
    denom <- denom+compute_gk(q)
  }
  ## sum over clusters
  logdenom <- log(denom/sqrt(pi))
  logdenom[!is.finite(logdenom)] <- -744.4  ## error handling to set 0s to just zero

  return(-sum(logdenom))


}


####################################################################
#                      Weighted  Score Funcion                     #
####################################################################
## negative because we minimize  negative log likelihood 
wglmm_score_CSCC <- function(betasig,clstr,X,Y,M_ki){ #M0,M1){
  ## gauss hermite weights + abscissas
  nodes <- 25
  h <- as.matrix(hermite.h.quadrature.rules(nodes)[[nodes]])[,1]
  w <- as.matrix(hermite.h.quadrature.rules(nodes)[[nodes]])[,2]
  
  ## define terms
  beta <- betasig[1:length(betasig)-1]
  sigma <- betasig[length(betasig)] ## glmmML version maximizes wrt sigma
  
  ## Xk0%*%beta
  Xkibeta <- X%*%beta
  
  
  compute_gkbeta <- function(q)  {

    ## compute g with weights
    gki <-  c(exp(Y*(Xkibeta+sqrt(2)*sigma*h[q])) /(1+exp(Xkibeta+sqrt(2)*sigma*h[q])) )^(M_ki)  #(Y*M1ki+(1-Y)*M0ki)
    gki[!is.finite(gki)] <- 0   ## error handling

    # compute g_beta (without the g)
    Y_expit <- c(Y-expit(Xkibeta+sqrt(2)*sigma*h[q]))
    A_beta <- (X*Y_expit)*(M_ki)  #(Y*M1ki+(1-Y)*M0ki) # well multiply this by g later

    A_sigma <- c(sqrt(2)*h[q]*(Y-expit(Xkibeta+sqrt(2)*sigma*h[q]) )*(M_ki))  #(Y*M1ki+(1-Y)*M0ki)  # well multiply this by g later
    
    df <- aggregate(cbind(log(gki),A_sigma,A_beta),by=list(c(clstr)), sum)[,1+(1:(length(beta)+2))]
    
    g_k <- exp(df[,1])
    
    return(list(g_k=w[q]*g_k,
                g_k_beta= w[q]*g_k*df[2+(1:length(beta))],
                g_k_sigma=w[q]*g_k*df[,2]))
    

  }
  

  ## loop over nodes
  num <- matrix(0,nrow=max(clstr),ncol=length(beta))
  num_sig <- rep(0,max(clstr))
  denom <- rep(0,max(clstr))
  for (q in 1:nodes){
    cq <- compute_gkbeta(q)
    num <- num+cq$g_k_beta
    num_sig <- num_sig+cq$g_k_sigma
    denom <- denom+cq$g_k
  }
  
  
  ## error handling
  denom[abs(denom)<1e-323] <- 1e-323 
  num[abs(num)<1e-323] <- 1e-323
  num_sig[abs(num_sig)<1e-323] <- 1e-323
  
  ## compute score
  Ub <- apply(num/denom,2,sum)
  Us <- sum(num_sig/denom)
  return(-c(Ub,Us))     ## return negative so optim will maximize
  
  
}






####################################################################
#                      CHEESE- Middle Variance                     #
####################################################################
## negative because we minimize  negative log likelihood 
cheese_CSCC <- function(betasig,clstr,X,Y,M0,M1){
  ## gauss hermite weights + abscissas
  nodes <- 25
  h <- as.matrix(hermite.h.quadrature.rules(nodes)[[nodes]])[,1]
  w <- as.matrix(hermite.h.quadrature.rules(nodes)[[nodes]])[,2]
  
  ## define terms
  beta <- betasig[1:length(betasig)-1]
  sigma <- betasig[length(betasig)] ## glmmML version maximizes wrt sigma  
  
  ## collect data
  #X <- model.matrix(g.glmm.test,data=sample)
  #Y <- model.frame(g.glmm.test)[,1]
  
  ## loop over clusters
  cheese <- 0
  for (k in (unique(clstr))){
    ## specify cluster data
    #kcluster <- sample[clstr==k,,drop=F]
    Xk <- X[c(clstr==k),,drop=FALSE]
    Yk <- Y[clstr==k,drop=FALSE]
    
    ## separate cases and controls
    Xk0 <- Xk[Yk==0,,drop=F]
    Xk1 <- Xk[Yk==1,,drop=F]
    
    Yk0 <- Yk[Yk==0,drop=F]
    Yk1 <- Yk[Yk==1,drop=F]
    
    ## loop over nodes for Gauss-Hermite
    num <- 0
    denom <- 0
    num_sig <- 0
    for (q in 1:nodes){
      ## calculate g=joint bernoulli pdfs
      g0 <- (prod( exp(Yk0*(Xk0%*%beta+sqrt(2)*sigma*h[q])) /(1+exp(Xk0%*%beta+sqrt(2)*sigma*h[q]))     ))
      g1 <- (prod( exp(Yk1*(Xk1%*%beta+sqrt(2)*sigma*h[q])) /(1+exp(Xk1%*%beta+sqrt(2)*sigma*h[q]))     ))
      
      ## HANDLE NEG INFINITIES
      if (!is.finite(g0)){
        if ( sum(Yk0*(Xk0%*%beta+sqrt(2)*sigma*h[q])) -sum(log(1+exp(Xk0%*%beta+sqrt(2)*sigma*h[q])) ) < (-800) ){
          g0 <- 0
        }
      }
      ## HANDLE NEG INFINITIES
      if (!is.finite(g1)){
        if ( sum(Yk1*(Xk1%*%beta+sqrt(2)*sigma*h[q])) -sum(log(1+exp(Xk1%*%beta+sqrt(2)*sigma*h[q])) ) < (-800) ){
          g1 <- 0
        }
      }
      
      g <- (g0^M0[k])*(g1^M1[k])
      
      ##  calculate other functions were integrating (ie numerator chunks)
      g_beta0 <- M0[k]*g*( t(Xk0)%*%(Yk0-expit(Xk0%*%beta+sqrt(2)*sigma*h[q])     ))
      g_beta1 <- M1[k]*g*( t(Xk1)%*%(Yk1-expit(Xk1%*%beta+sqrt(2)*sigma*h[q])     ))
      g_beta <- g_beta0+g_beta1
      
      ## HERE WE HAVE MULTIPLIED BY SIGMA DUE TO THE CHAIN RULE NECESSARY (because we are using log to constrain the sign of sigma)
      g_sig0 <- M0[k]*g*( sum(sqrt(2)*h[q]*(Yk0-expit(Xk0%*%beta+sqrt(2)*sigma*h[q])     )))
      g_sig1 <- M1[k]*g*( sum(sqrt(2)*h[q]*(Yk1-expit(Xk1%*%beta+sqrt(2)*sigma*h[q])     )))
      g_sig <- g_sig0+g_sig1
      
      ## sum elements of approximated integrals
      num <- num+w[q]*g_beta
      num_sig <- num_sig+w[q]*g_sig
      denom <- denom+w[q]*g      
      
    }
    ### MAYBE?
    for (i in 1:length(num)){
      if (abs(num[i])<1e-323){num[i]=1e-323 }
    }
    if (abs(num_sig)<1e-323){num_sig=1e-323 }
    if (abs(denom)<1e-323){denom=1e-323 }
    ## HAVE NOT ACCOUNTED FOR THE TRICKY CASES OF NA/INFTY/0
    cheese <- cheese +as.matrix(c((num/denom),(num_sig/denom)))%*%as.matrix(t(c((num/denom),(num_sig/denom))))
  }
  ## return the negative so that optim will maximize
  return(cheese)
  
}



  

####################################################################
#                          FULL FUNCTION                           #
####################################################################

wglmm_robust <- function(full,sample,frmla,testid,param,SRS=FALSE){
  clstr <- sample[,testid]
  full_clstr <- full[,testid]
  
  # must reassign name of frmla in order for it to be taken into other functions
  form <- formula(frmla)
  fullY <- model.frame(formula=form,data=full)[,1]
  fullX <- model.matrix(object=form,data=full)  
  
  Y <- model.frame(formula=form,data=sample)[,1]
  X <- model.matrix(object=form,data=sample)
  
  # Weights
  if (SRS==TRUE){
    ## SRS Weights
    full_tibble <- as_tibble(cbind(fullY,full_clstr))
    sample_tibble <- as_tibble(cbind(Y,clstr))

    N_k <- full_tibble %>% group_by(full_clstr) %>% summarise(Nk=n())
    n_k <- sample_tibble %>% group_by(clstr) %>% summarise(nk=n())

    M0 <- M1 <- c(N_k$Nk/n_k$nk)

  } else {
    ##   CSCC Weights   ##
    full_tibble <- as_tibble(cbind(fullY,full_clstr))
    sample_tibble <- as_tibble(cbind(Y,clstr))  
    
    N_yk <- full_tibble %>% group_by(full_clstr) %>% summarise(N0k=sum(1-fullY),N1k=sum(fullY))
    n_yk <- sample_tibble %>% group_by(clstr) %>% summarise(n0k=sum(1-Y),n1k=sum(Y))
    
    M0 <- c(N_yk$N0k/n_yk$n0k)  
    M1 <- c(N_yk$N1k/n_yk$n1k)
  }
    ## weights
    M0_ki <- M0[clstr] ## assign to cases and controls separately
    M1_ki <- M1[clstr]
    M_ki <- M1_ki*Y+M0_ki*(1-Y) ## make observation-level weights
    
    
    
    chCC <- cheese_CSCC(param,clstr=clstr,X=X,Y=Y,M0=M0,M1=M1)
    invI <- try(solve(optimHess(par=param,wglmm_loglik_CSCC,gr=wglmm_score_CSCC,clstr=clstr,X=X,Y=Y,M_ki=M_ki)),silent=T)
    return(list(chCC,invI))

  
}




