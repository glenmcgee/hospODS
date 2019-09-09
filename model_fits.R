####################################################################
#                                                                  #
#                            Fit Models                            #
#                                                                  #
####################################################################
#                               April 3                            #
####################################################################
## Compute SRRs for a given model fit
library(tidyverse)
library(glmmML)


####################################################################
#                             Fit with ML                          #
#                       For: Full, RS, CS, SRS                     #
####################################################################

fit_ML <- function(samp,fixed_only=FALSE){      ## TESTING THE ML METHOD ON SRS
  # make data frame
  samp <- data.frame(samp)
  
  
  reg <- glmmML(event~xw.cont+xw.bin+xbetween,cluster=clinic,family=binomial,data=samp)
  est <- c(reg$coefficients,reg$sigma)
  SE <- c(reg$coef.sd,reg$sigma.sd) 
  ## error handling:
  if (length(SE)<length(est)){SE <- rep(NA,length(est))}
  pred <- reg$posterior.modes
  
  
  if (fixed_only==TRUE){    ## report only fixed effects
    return(c(nrow(samp),est,SE))
  } else{                   ## report fixed effects + MSEP
    
    SRR <- compute_SRR(reg$coefficients,reg$sigma,pred,samp,'event~xw.cont+xw.bin+xbetween','clinic')
          # old: compute_SRR(est,pred,samp,'event~xw.cont+xw.bin+xbetween','clinic')
    MSEP <- mean((pred-true_rand)^2)
    MSEP_SRR <- mean((SRR-true_SRR)^2)
    misclass <- sum((sign(SRR-1)!=sign(true_SRR-1)))
    
    return(c(nrow(samp),est,SE,MSEP,MSEP_SRR,misclass))
  }
  
}


####################################################################
#                    Fit CSCC with WL or Offset                    #
####################################################################
fit_CSCC <- function(samp,offset_method=FALSE,naive_method=FALSE,fixed_only=FALSE,scaled=FALSE,fullsizes=NA,SRS=FALSE){
  ## default is weighted (pseudolikelihood)
  ## offset=TRUE gives offset method
  
  full_tibble <- as_tibble(full)
  sample_tibble <- as_tibble(samp)  
  
  if(SRS==TRUE){
    ## SRS weights
    N_k <- full_tibble %>% group_by(clinic) %>% summarise(Nk=n())
    n_k <- sample_tibble %>% group_by(clinic) %>% summarise(nk=n())
    
    Mk <- c(N_k$Nk/n_k$nk)  
    
    M_ki <- Mk[samp$clinic]
    log_offset <- rep(0,length(M_ki))
  }else{
    ## CSCC 
    N_yk <- full_tibble %>% group_by(clinic) %>% summarise(N0k=sum(1-event),N1k=sum(event))
    n_yk <- sample_tibble %>% group_by(clinic) %>% summarise(n0k=sum(1-event),n1k=sum(event))
    
    M0 <- c(N_yk$N0k/n_yk$n0k)  
    M1 <- c(N_yk$N1k/n_yk$n1k)
    
    ###################################
    ## weights
    M0_ki <- M0[samp$clinic] ## assign to cases and controls separately
    M1_ki <- M1[samp$clinic]
    M_ki <- M1_ki*samp$event+M0_ki*(1-samp$event) ## make observation-level weights
    
    ## Method 2 scaling weights (nk):
    if(scaled==TRUE){
      a_ki <- rep(1,length(M_ki))
      for (ii in unique(samp$clinic)){
        a_ki[samp$clinic==ii] <- sum(samp$clinic==ii)/sum(M_ki[samp$clinic==ii])
      }
      M_ki <- a_ki*M_ki
    }
    # ## Nk Alternative Method 2 scaling weights: 
    # if(scaled==TRUE){
    #   a_ki <- rep(1,length(M_ki))
    #   for (ii in unique(samp$clinic)){
    #     a_ki[samp$clinic==ii] <- sum(samp$clinic==ii)/(M_ki[samp$clinic==ii])
    #   }
    #   M_ki <- a_ki*M_ki
    # }
    
    ###################################
    ## offset terns
    logoff <- log(M0/M1)
    log_offset <- logoff[samp$clinic]   ## make observation-level offset terms
    
  }

  

  

  

  # make data frame
  samp <- data.frame(samp)
  
  
  if (offset_method==TRUE){
    ## fit model
    reg <- glmmML(event~xw.cont+xw.bin+xbetween,cluster=clinic,offset=log_offset,family=binomial,data=samp)
    ## collect estimates
    est <- c(reg$coefficients,reg$sigma)
    SE <- c(reg$coef.sd,reg$sigma.sd) 
  }else if(naive_method==TRUE){
    reg <- glmmML(event~xw.cont+xw.bin+xbetween,cluster=clinic,family=binomial,data=samp)
    ## collect estimates
    est <- c(reg$coefficients,reg$sigma)
    SE <- c(reg$coef.sd,reg$sigma.sd) 
  }else if(SRS==TRUE){
    ## fit model
    reg <- glmmML(event~xw.cont+xw.bin+xbetween,cluster=clinic,weights=M_ki,family=binomial,data=samp)
    ## collect estimates
    est <- c(reg$coefficients,reg$sigma)
    #invI <- reg$variance
    robust <- wglmm_robust(full,samp,'(event~xw.cont+xw.bin+xbetween)',testid='clinic',est,SRS=TRUE)
    cheese <- robust[[1]]
    invI <- robust[[2]]
    SE <- sqrt(diag(invI%*%cheese%*%invI)) 
    #SE <- c(reg$coef.sd,reg$sigma.sd)  # model based (naive)
  }else {
    ## fit model
    reg <- glmmML(event~xw.cont+xw.bin+xbetween,cluster=clinic,weights=M_ki,family=binomial,data=samp)
    ## collect estimates
    est <- c(reg$coefficients,reg$sigma)
    #invI <- reg$variance
    robust <- wglmm_robust(full,samp,'(event~xw.cont+xw.bin+xbetween)',testid='clinic',est)
    cheese <- robust[[1]]
    invI <- robust[[2]]
    SE <- sqrt(diag(invI%*%cheese%*%invI)) 
    #SE <- c(reg$coef.sd,reg$sigma.sd)  # model based (naive)
  }
  

  ## error handling:
  if (length(SE)<length(est)){SE <- rep(NA,length(est))}
  pred <- reg$posterior.modes
  
  
  
  if (fixed_only==TRUE){    ## report only fixed effects
    return(c(nrow(samp),est,SE))
  } else{                   ## report fixed effects + MSEP
    
    ## compute SRRs                                                           ## offset=0 if weighted method, =log_offset if offset method
    SRR <- compute_SRR(reg$coefficients,reg$sigma,pred,samp,'event~xw.cont+xw.bin+xbetween','clinic',offset = offset_method*log_offset,Mki=M_ki)
    # old: compute_SRR(est,pred,samp,'event~xw.cont+xw.bin+xbetween','clinic',offset = offset_method*log_offset)
    MSEP <- mean((pred-true_rand)^2)
    MSEP_SRR <- mean((SRR-true_SRR)^2)
    misclass <- sum((sign(SRR-1)!=sign(true_SRR-1)))
    
    return(c(nrow(samp),est,SE,MSEP,MSEP_SRR,misclass))
  }
  
  

}






####################################################################
#                     Fit CC with WL or Offset                     #
####################################################################
fit_CC <- function(samp,offset_method=FALSE,fixed_only=TRUE){   
  ## fitting CC with offset method
  ## weighted method is commented out (not as good)
  
  N1 <- sum(full$event); N0 <- sum(1-full$event)
  n1 <- sum(samp$event); n0 <- sum(1-samp$event)
  
  ####################
  ## weights
  M1 <- N1/n1; M0 <- N0/n0
  ## obs level
  M_ki <- M1*samp$event+M0*(1-samp$event)
  
  ####################
  ## offset
  logoff <- log((n1/N1)/(n0/N0))
  log_offset <- rep(logoff,nrow(samp))  ## make observation-level offset terms
  
  ## make data frame
  samp <- data.frame(samp)
  
  if (offset_method==TRUE){
    ## fit model
    reg <- glmmML(event~xw.cont+xw.bin+xbetween,cluster=clinic,offset=log_offset,family=binomial,data=samp)
  }else {
    ## fit model
    reg <- glmmML(event~xw.cont+xw.bin+xbetween,cluster=clinic,weights=M_ki,family=binomial,data=samp)
  }
  
  est <- c(reg$coefficients,reg$sigma)
  SE <- c(reg$coef.sd,reg$sigma.sd) 
  ## error handling:
  if (length(SE)<length(est)){SE <- rep(NA,length(est))}

  
  if (fixed_only==TRUE){    ## report only fixed effects
    return(c(nrow(samp),est,SE))
  } else{                   ## report fixed effects + MSEP
    print("Doesn't make sense to make predictions for a sample that doesn't draw from every cluster")
  }
  
  
  
}