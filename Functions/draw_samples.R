####################################################################
#                                                                  #
#                         Sampling Functions                       #
#                                                                  #
####################################################################
#                              Feb 2017                            #
####################################################################
## code to draw samples

library(tidyverse)

####################################################################
#                   Function to Draw CSCC Sample                   #
####################################################################

draw_CSCC_each <- function(dat,nk){
  ncontrols <- ncases <- nk/2
  dat <- as_tibble(dat)
  
  sum_dat <- dat %>% 
    group_by(clinic) %>% 
    summarise(#Nk=n(),
              Nk1=sum(event),
              Nk0=sum(1-event),
              #nk=min(Nk,ncases+ncontrols),
              
              nk1=Nk1*(Nk1<ncases) + ## if there arent enough cases, take them all
                ncases*(Nk1>=ncases & Nk0>=ncontrols) + ## if there are enough of both, only take whats needed
                min(Nk1,ncases+ncontrols-Nk0)*(Nk1>=ncases & Nk0<ncontrols), ## if we have enough cases but not enough controls, take extra
              
              nk0=Nk0*(Nk0<ncontrols) + ## if there arent enough controls, take them all
                ncontrols*(Nk1>=ncases & Nk0>=ncontrols) + ## if there are enough of both, only take whats needed
                min(Nk0,ncontrols+ncases-Nk1)*(Nk0>=ncontrols & Nk1<ncases)#, ## if we have enough controls but not enough cases, take extra
              
              #checknk=((nk1+nk0)==nk), ## check that were sampling the right amount
              #checknky=((nk1<=Nk1) & (nk0<=Nk0)) ## check that were not asking for too many
    )
  
  #sum_dat
  
  ## assign everyone a within strata sampling id
  dat <- dat %>% group_by(clinic, event) %>% mutate(sample_id=sample(n()))
  ## assign the controls a sampling limit
  dat$nk <- sum_dat$nk0[dat$clinic] 
  ## assign cases a sampling limit
  dat$nk[dat$event==1] <- (sum_dat$nk1[dat$clinic])[dat$event==1]
  
  
  dat <- dat %>% filter(sample_id<=nk)  ## if sampling id is less than or equal to sampling limit for
                                        ## event/clinic strata, then sample it
  
  return(dat)
}



####################################################################
#                   Function to Draw CSCC Sample                   #
#                  Proportional to Cluster Sizes                   #
####################################################################

draw_CSCCpropto_each <- function(dat,frac){

  dat <- as_tibble(dat)
  
  sum_dat <- dat %>% 
    group_by(clinic) %>% 
    summarise(Nk=n(),
      nyk=round(Nk*frac/2,0),
      Nk1=sum(event),
      Nk0=sum(1-event),
      
      nk1=Nk1*(Nk1<nyk) + ## if there arent enough cases, take them all
        nyk*(Nk1>=nyk & Nk0>=nyk) + ## if there are enough of both, only take whats needed
        min(Nk1,nyk+nyk-Nk0)*(Nk1>=nyk & Nk0<nyk), ## if we have enough cases but not enough controls, take extra
      
      nk0=Nk0*(Nk0<nyk) + ## if there arent enough controls, take them all
        nyk*(Nk1>=nyk & Nk0>=nyk) + ## if there are enough of both, only take whats needed
        min(Nk0,nyk+nyk-Nk1)*(Nk0>=nyk & Nk1<nyk)#, ## if we have enough controls but not enough cases, take extra
     )
  

  
  ## assign everyone a within strata sampling id
  dat <- dat %>% group_by(clinic, event) %>% mutate(sample_id=sample(n()))
  ## assign the controls a sampling limit
  dat$nk <- sum_dat$nk0[dat$clinic] 
  ## assign cases a sampling limit
  dat$nk[dat$event==1] <- (sum_dat$nk1[dat$clinic])[dat$event==1]
  
  
  dat <- dat %>% filter(sample_id<=nk)  ## if sampling id is less than or equal to sampling limit for
  ## event/clinic strata, then sample it
  
  return(dat)
}


####################################################################
#                   Function to Draw CSCC Sample                   #
#                 Drawing More Heavily at the Margin               #
####################################################################
## split clusters into 3 tertiles based on raw risk
## lowest and highest tertiles get rates of 0.75*frac
## middle tier gets 1.5*frac

## same input as CSCCbalanced
draw_CSCCmargin_each <- function(dat,nk){
  nk <- nk/2 ## number of cases and controls
  dat <- as_tibble(dat)
  

  sum_risk <- dat %>% 
    group_by(clinic) %>% 
    summarise(risk_k=mean(event))
  
  risk_cat <- cut(sum_risk$risk_k, breaks = quantile(sum_risk$risk_k, probs = seq(0, 1, 1/3)),include.lowest=TRUE,labels=c(1,2,3))    ## split into tertiles
  dat$risk_cat <- risk_cat[dat$clinic]
  
  sum_dat <- dat %>% 
    group_by(clinic) %>% 
    summarise(Nk1=sum(event),
              Nk0=sum(1-event),
              nyk=round(nk*(0.75*(first(risk_cat)!=2)+1.5*(first(risk_cat)==2)   )),
              nk1=Nk1*(Nk1<nyk) + ## if there arent enough cases, take them all
              nyk*(Nk1>=nyk & Nk0>=nyk) + ## if there are enough of both, only take whats needed
              min(Nk1,nyk+nyk-Nk0)*(Nk1>=nyk & Nk0<nyk), ## if we have enough cases but not enough controls, take extra
              
              nk0=Nk0*(Nk0<nyk) + ## if there arent enough controls, take them all
                nyk*(Nk1>=nyk & Nk0>=nyk) + ## if there are enough of both, only take whats needed
                min(Nk0,nyk+nyk-Nk1)*(Nk0>=nyk & Nk1<nyk)#, ## if we have enough controls but not enough cases, take extra
    )
        
  
  
  ## assign everyone a within strata sampling id
  dat <- dat %>% group_by(clinic, event) %>% mutate(sample_id=sample(n()))
  ## assign the controls a sampling limit
  dat$nk <- sum_dat$nk0[dat$clinic] 
  ## assign cases a sampling limit
  dat$nk[dat$event==1] <- (sum_dat$nk1[dat$clinic])[dat$event==1]
  
  
  dat <- dat %>% filter(sample_id<=nk)  ## if sampling id is less than or equal to sampling limit for
  ## event/clinic strata, then sample it
  
  return(dat)
}

####################################################################
#                   Function to Draw CSCC Sample                   #
#                 Drawing More Heavily at the Margin               #
#                           Based on SRR                           #
####################################################################
## split clusters into 3 tertiles based on raw risk
## lowest and highest tertiles get rates of 0.75*frac
## middle tier gets 1.5*frac

## same input as CSCCbalanced
draw_CSCCmarginSRR_each <- function(dat,SRRs,nk){
  nk <- nk/2 ## number of cases and controls
  dat <- as_tibble(dat)
  
  
  SRR_cat <- cut(SRRs, breaks = quantile(SRRs, probs = seq(0, 1, 1/3)),include.lowest=TRUE,labels=c(1,2,3))    ## split into tertiles
  dat$SRR_cat <- SRR_cat[dat$clinic]
  
  sum_dat <- dat %>% 
    group_by(clinic) %>% 
    summarise(Nk1=sum(event),
              Nk0=sum(1-event),
              nyk=round(nk*(0.75*(first(SRR_cat)!=2)+1.5*(first(SRR_cat)==2)   )),
              nk1=Nk1*(Nk1<nyk) + ## if there arent enough cases, take them all
                nyk*(Nk1>=nyk & Nk0>=nyk) + ## if there are enough of both, only take whats needed
                min(Nk1,nyk+nyk-Nk0)*(Nk1>=nyk & Nk0<nyk), ## if we have enough cases but not enough controls, take extra
              
              nk0=Nk0*(Nk0<nyk) + ## if there arent enough controls, take them all
                nyk*(Nk1>=nyk & Nk0>=nyk) + ## if there are enough of both, only take whats needed
                min(Nk0,nyk+nyk-Nk1)*(Nk0>=nyk & Nk1<nyk)#, ## if we have enough controls but not enough cases, take extra
    )
  
  
  
  ## assign everyone a within strata sampling id
  dat <- dat %>% group_by(clinic, event) %>% mutate(sample_id=sample(n()))
  ## assign the controls a sampling limit
  dat$nk <- sum_dat$nk0[dat$clinic] 
  ## assign cases a sampling limit
  dat$nk[dat$event==1] <- (sum_dat$nk1[dat$clinic])[dat$event==1]
  
  
  dat <- dat %>% filter(sample_id<=nk)  ## if sampling id is less than or equal to sampling limit for
  ## event/clinic strata, then sample it
  
  return(dat)
}


####################################################################
#                   Function to Draw SRS Sample                    #
####################################################################
draw_SRS_each <- function(dat,n_k){
  dat <- as_tibble(dat)
  
  sum_dat <- dat %>% 
    group_by(clinic) %>% 
    summarise(Nk=n(),
      nk=min(Nk,n_k) )
  
  ## assign everyone a within strata sampling id
  dat <- dat %>% group_by(clinic) %>% mutate(sample_id=sample(n()))
  
  ## assign sampling limit
  dat$nk <- sum_dat$nk[dat$clinic] 

  dat <- dat %>% filter(sample_id<=nk)  ## if sampling id is less than or equal to sampling limit for
  ## event/clinic strata, then sample it
  
  return(dat)
}



####################################################################
#                   Function to Draw SRS Sample                    #
#                  Proportional to Cluster Sizes                   #
####################################################################
draw_SRSpropto_each <- function(dat,frac){
  
  dat <- as_tibble(dat)
  
  sum_dat <- dat %>% 
    group_by(clinic) %>% 
    summarise(Nk=n(),
              nk=round(frac*Nk,0) )
  
  ## assign everyone a within strata sampling id
  dat <- dat %>% group_by(clinic) %>% mutate(sample_id=sample(n()))
  
  ## assign sampling limit
  dat$nk <- sum_dat$nk[dat$clinic] 
  
  dat <- dat %>% filter(sample_id<=nk)  ## if sampling id is less than or equal to sampling limit for
  ## event/clinic strata, then sample it
  
  return(dat)
}

####################################################################
#                   Function to Draw SRS Sample                    #
#                  Drawing More Heavily at the Margin              #
####################################################################
draw_SRSmargin_each <- function(dat,nnk){
  
  dat <- as_tibble(dat)
  
  sum_risk <- dat %>% 
    group_by(clinic) %>% 
    summarise(risk_k=mean(event))
  
  risk_cat <- cut(sum_risk$risk_k, breaks = quantile(sum_risk$risk_k, probs = seq(0, 1, 1/3)),include.lowest=TRUE,labels=c(1,2,3))    ## split into tertiles
  dat$risk_cat <- risk_cat[dat$clinic]
  
  
  
  sum_dat <- dat %>% 
    group_by(clinic) %>% 
    summarise(nk=round(nnk*(0.75*(first(risk_cat)!=2)+1.5*(first(risk_cat)==2)   )) )
  
  ## assign everyone a within strata sampling id
  dat <- dat %>% group_by(clinic) %>% mutate(sample_id=sample(n()))
  
  ## assign sampling limit
  dat$nk <- sum_dat$nk[dat$clinic] 
  
  dat <- dat %>% filter(sample_id<=nk)  ## if sampling id is less than or equal to sampling limit for
  ## event/clinic strata, then sample it
  
  return(dat)
}

####################################################################
#                   Function to Draw SRS Sample                    #
#                  Drawing More Heavily at the Margin              #
#                           Based on SRR                           #
####################################################################
draw_SRSmarginSRR_each <- function(dat,SRRs,nnk){
  
  dat <- as_tibble(dat)
  
  
  SRR_cat <- cut(SRRs, breaks = quantile(SRRs, probs = seq(0, 1, 1/3)),include.lowest=TRUE,labels=c(1,2,3))    ## split into tertiles
  dat$SRR_cat <- SRR_cat[dat$clinic]
  
  
  
  sum_dat <- dat %>% 
    group_by(clinic) %>% 
    summarise(nk=round(nnk*(0.75*(first(SRR_cat)!=2)+1.5*(first(SRR_cat)==2)   )) )
  
  ## assign everyone a within strata sampling id
  dat <- dat %>% group_by(clinic) %>% mutate(sample_id=sample(n()))
  
  ## assign sampling limit
  dat$nk <- sum_dat$nk[dat$clinic] 
  
  dat <- dat %>% filter(sample_id<=nk)  ## if sampling id is less than or equal to sampling limit for
  ## event/clinic strata, then sample it
  
  return(dat)
}

####################################################################
#             Function to Conduct Cluster Sampling                 #
####################################################################
draw_CS_each <- function(dat,K){
  ## draw K clusters
  clust_id <- sample(unique(dat$clinic),K,replace=FALSE)
  
  ## make tibble
  dat <- as_tibble(dat)
  
  ## select only those with sampled cluster
  dat <- dat %>% filter(clinic %in% clust_id)
  
  return(dat)
}

####################################################################
#                   Function to Draw Random Sample                 #
####################################################################
draw_RS_each <- function(dat,nn){
  
  ## make tibble
  dat <- as_tibble(dat)
  
  ## assign everyone a sampling id
  dat <- dat %>% mutate(sample_id=sample(n()))
  

  ## select only those sampled
  dat <- dat %>% filter(sample_id<=nn) %>% arrange(clinic)
  #dat$clinic <- as.numeric(factor(clean_full$clinic))
  return(dat)
}

####################################################################
#                     Function to Draw CC Sample                  #
####################################################################
draw_CC_each <- function(dat,nn){
  
  ncontrols <- ncases <- round(nn/2)
  dat <- as_tibble(dat)
  
  sum_dat <- dat %>% 
    summarise(#Nk=n(),
      N1=sum(event),
      N0=sum(1-event),

      
      n1=N1*(N1<ncases) + ## if there arent enough cases, take them all
        ncases*(N1>=ncases & N0>=ncontrols) + ## if there are enough of both, only take whats needed
        min(N1,ncases+ncontrols-N0)*(N1>=ncases & N0<ncontrols), ## if we have enough cases but not enough controls, take extra
      
      n0=N0*(N0<ncontrols) + ## if there arent enough controls, take them all
        ncontrols*(N1>=ncases & N0>=ncontrols) + ## if there are enough of both, only take whats needed
        min(N0,ncontrols+ncases-N1)*(N0>=ncontrols & N1<ncases)#, ## if we have enough controls but not enough cases, take extra
      
      #checknk=((nk1+nk0)==nk), ## check that were sampling the right amount
      #checknky=((nk1<=Nk1) & (nk0<=Nk0)) ## check that were not asking for too many
    )
  
  
  ## assign everyone a within strata sampling id
  dat <- dat %>% group_by(event) %>% mutate(sample_id=sample(n()))
  ## assign the controls a sampling limit
  dat$nny <- sum_dat$n0
  ## assign cases a sampling limit
  dat$nny[dat$event==1] <- (sum_dat$n1)
  
  ## if sampling id is less than or equal to sampling limit for
  dat <- dat %>% filter(sample_id<=nny)  
  
  
  ## old method that doesnt sample extra controls when lacking cases 
  # nny <- round(nn/2)
  # ## make tibble
  # dat <- as_tibble(dat)
  # 
  # ## assign everyone a sampling id
  # dat <- dat %>% group_by(event) %>% mutate(sample_id=sample(n()))
  # 
  # ## select only those sampled
  # dat <- dat %>% filter(sample_id<=nny) %>% arrange(clinic)
  # #dat$clinic <- as.numeric(factor(clean_full$clinic))
  
  return(dat)
}
