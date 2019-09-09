
####################################################################
#                                                                  #
#                  Simulation Study 3 for BASE Model               #
#                                                                  #
####################################################################
#                         Simulated Data                           #
####################################################################
#                             March 2017                           #
####################################################################

## necessary packages
library(glmmML)
source('draw_samples.R')
source('compute_SRR.R')
source('wglmm_test.R')
source('model_fits.R')
samp_sizes <- read.table("samp_sizes.txt",header=TRUE)[,1] ## cluster


## Set Parameters
Kclust <- 150   ## K: number of clusters
set.seed(1234); samp_sizes <- sample(samp_sizes,Kclust)  ## draw relevant sample sizes
N <- sum(samp_sizes)
ev.rate <- 0.2     ## Event rate
sigma <- 0.2    ## sigma: SD of random effects
bin.risk <- 0.2 ## % dual eligible

## sampling fractions
samp_frac <- c(0.10,0.20,0.30,0.50)
nk_avg <- mean(samp_sizes); nk_avg
samp_nk <- round(samp_frac*nk_avg); samp_nk
samp_n <- samp_frac*N
samp_K <- round(samp_frac*Kclust); samp_K

## set seed for reproducibility
set.seed(1000)
num_iter <- 2000

## expit
expit <- function(y){
  return(exp(y)/(1+exp(y)))
}

####################################################################
#                       Generating Covariates                      #
####################################################################

## generating covariate data:
clinic <- rep(1:Kclust,times=samp_sizes)   ## cluster identifier
event <- rep(0,length(clinic)) ## placeholder until simulation begins
xw.cont <- rnorm(N,0,2) ## within cluster variable
xw.bin <- rbinom(N,1,0.2) ## within cluster binary
xbetween <- rep(runif(Kclust,0,1),times=samp_sizes)  ## between cluster covariate: uniform on 0 to 1
full <- data.frame(cbind(clinic,event,xw.cont,xw.bin,xbetween))


## regression parameters:
beta.xw.cont <- 0.2  # 0.1
beta.xw.bin <- 0.2   # 0.2
beta.xbetween <- 0.5  # 0.5

## get intercept ( should vary with event rate )
beta.int <- -1.7; beta_covs <- c(beta.xw.cont,beta.xw.bin,beta.xbetween)
step <- 0.00005
while(abs(ev.rate-mean(expit(beta.int+cbind(xw.cont,xw.bin,xbetween)%*%beta_covs)))>0.0001){
  if ((ev.rate-mean(expit(beta.int+cbind(xw.cont,xw.bin,xbetween)%*%beta_covs)))<0){
    beta.int <- beta.int-step  
  } else {
    beta.int <- beta.int+step    
  }
}

true_coefs <- c(beta.int,beta_covs)
true_sigma <- sigma

## check coefs
c(true_coefs,true_sigma)
## check risk
ev.rate


## compute fitted values for true SRRs
linpred <- cbind(1,xw.cont,xw.bin,xbetween)%*%true_coefs    ##Xbeta





####################################################################
#                        Simulation - NonParallel                     #
####################################################################

df <- c()
for(jj in (1:num_iter)){
  try({
    
  set.seed(1234+99*jj)
  
  ####################################################################
  #                    Data Generating Mechanism                     #
  ####################################################################
  ## draw random intercepts
  true_rand <- rnorm(Kclust,0,true_sigma)
  full$true_rand_int <- rep(true_rand,times=samp_sizes)
  
  ## draw outcomes
  ## note the 1 on the end is since we multiply by random intercept
  ## -c(1,2) are to exclude clinic id and outcome
  full$event <- rbinom(nrow(full),1,expit(true_coefs[1] + (data.matrix(full)[,-c(1,2)])%*%c(true_coefs[-1],1) ))  
  
  ## true SRRs
  true_SRR <- compute_SRR(true_coefs,true_sigma,true_rand,full,'event~xw.cont+xw.bin+xbetween','clinic')
  ## old:
  # denom_ki <- c(expit(linpred+0))
  # num_ki <- c(expit(linpred+full$true_rand_int))
  # 
  # df <- as_tibble(cbind(denom_ki,num_ki,clinic))
  # df <- df %>% group_by(clinic) %>% summarise(num=sum(num_ki),
  #                                             denom=sum(denom_ki),
  #                                             SRR=num/denom)
  # true_SRR <- c(df$SRR)
  
  
  ####################################################################
  #                         Full Data Naive                          #
  #################################################################### 
  naivereg <- glmmML(event~xw.cont+xbetween,cluster=clinic,family=binomial,data=full) ## without +xw.bin
  naiveest <- c(naivereg$coefficients,naivereg$sigma)
  naiveSE <- c(naivereg$coef.sd,naivereg$sigma.sd) 
  #prednaive <- EBpred_naive(naiveest,naive,'event~xw.cont+xw.bin+xbetween','clinic')
  prednaive <- naivereg$posterior.modes
  SRRnaive <- compute_SRR(naivereg$coefficients,naivereg$sigma,prednaive,full,'event~xw.cont+xbetween','clinic')
  # old compute_SRR(naiveest,prednaive,full,'event~xw.cont+xbetween','clinic') ## +xw.bin
  MSPEnaive <- mean((prednaive-true_rand)^2)
  SRR_MSEPnaive <- mean((SRRnaive-true_SRR)^2)
  misclassnaive <- sum((sign(SRRnaive-1)!=sign(true_SRR-1)))
  
  results_naive <- c(nrow(full),naiveest,naiveSE,MSPEnaive,SRR_MSEPnaive,misclassnaive)
  
  ####################################################################
  #                            Draw Samples                          #
  ####################################################################
  
  ## Random Sample
  set.seed(1234+99*jj) ## set same seed each time so that samples are cumulative
  sample_RS_10 <- draw_RS_each(full,samp_n[1])  ## 10%
  set.seed(1234+99*jj)
  sample_RS_20 <- draw_RS_each(full,samp_n[2])  ## 20%
  set.seed(1234+99*jj)
  sample_RS_30 <- draw_RS_each(full,samp_n[3])  ## 30%
  set.seed(1234+99*jj)
  sample_RS_50 <- draw_RS_each(full,samp_n[4])  ## 50%
  
  ## CC Sample
  set.seed(1234+99*jj) ## set same seed each time so that samples are cumulative
  sample_CC_10 <- draw_CC_each(full,samp_n[1])  ## 10%
  set.seed(1234+99*jj)
  sample_CC_20 <- draw_CC_each(full,samp_n[2])  ## 20%
  set.seed(1234+99*jj)
  sample_CC_30 <- draw_CC_each(full,samp_n[3])  ## 30%
  set.seed(1234+99*jj)
  sample_CC_50 <- draw_CC_each(full,samp_n[4])  ## 50%
  
  ## Cluster Sampling
  set.seed(1234+99*jj) ## set same seed each time so that samples are cumulative
  sample_CS_10 <- draw_CS_each(full,samp_K[1])  ## 10%
  set.seed(1234+99*jj)
  sample_CS_20 <- draw_CS_each(full,samp_K[2])  ## 20%
  set.seed(1234+99*jj)
  sample_CS_30 <- draw_CS_each(full,samp_K[3])  ## 30%
  set.seed(1234+99*jj)
  sample_CS_50 <- draw_CS_each(full,samp_K[4])  ## 50%
  
  ## CSCC balanced
  set.seed(1234+99*jj) ## set same seed each time so that samples are cumulative
  sample_CSCC_10 <- draw_CSCC_each(full,samp_nk[1])  ## n1k=10
  set.seed(1234+99*jj)
  sample_CSCC_20 <- draw_CSCC_each(full,samp_nk[2])  ## n1k=20
  set.seed(1234+99*jj)
  sample_CSCC_30 <- draw_CSCC_each(full,samp_nk[3])  ## n1k=30
  set.seed(1234+99*jj)
  sample_CSCC_50 <- draw_CSCC_each(full,samp_nk[4])  ## n1k=50
  
  ## CSCC proportional to Nk
  set.seed(1234+99*jj) ## set same seed each time so that samples are cumulative
  sample_CSCCpropto_10 <- draw_CSCCpropto_each(full,samp_frac[1])  ## 10%
  set.seed(1234+99*jj)
  sample_CSCCpropto_20 <- draw_CSCCpropto_each(full,samp_frac[2])  ## 20%
  set.seed(1234+99*jj)
  sample_CSCCpropto_30 <- draw_CSCCpropto_each(full,samp_frac[3])  ## 30%
  set.seed(1234+99*jj)
  sample_CSCCpropto_50 <- draw_CSCCpropto_each(full,samp_frac[4])  ## 50%
  
  ## CSCC at the margin 
  set.seed(1234+99*jj) ## set same seed each time so that samples are cumulative
  sample_CSCCmargin_10 <- draw_CSCCmargin_each(full,samp_nk[1])  ## n1k=10
  set.seed(1234+99*jj)
  sample_CSCCmargin_20 <- draw_CSCCmargin_each(full,samp_nk[2])  ## n1k=20
  set.seed(1234+99*jj)
  sample_CSCCmargin_30 <- draw_CSCCmargin_each(full,samp_nk[3])  ## n1k=30
  set.seed(1234+99*jj)
  sample_CSCCmargin_50 <- draw_CSCCmargin_each(full,samp_nk[4])  ## n1k=50

  ## CSCC at the margin (SRR)
  set.seed(1234+99*jj) ## set same seed each time so that samples are cumulative
  sample_CSCCmarginSRR_10 <- draw_CSCCmarginSRR_each(full,SRRnaive,samp_nk[1])  ## n1k=10
  set.seed(1234+99*jj)
  sample_CSCCmarginSRR_20 <- draw_CSCCmarginSRR_each(full,SRRnaive,samp_nk[2])  ## n1k=20
  set.seed(1234+99*jj)
  sample_CSCCmarginSRR_30 <- draw_CSCCmarginSRR_each(full,SRRnaive,samp_nk[3])  ## n1k=30
  set.seed(1234+99*jj)
  sample_CSCCmarginSRR_50 <- draw_CSCCmarginSRR_each(full,SRRnaive,samp_nk[4])  ## n1k=50
  
  ## SRS balanced
  set.seed(1234+99*jj) ## set same seed each time so that samples are cumulative
  sample_SRS_10 <- draw_SRS_each(full,samp_nk[1])  ## 10%
  set.seed(1234+99*jj)
  sample_SRS_20 <- draw_SRS_each(full,samp_nk[2])  ## 20%
  set.seed(1234+99*jj)
  sample_SRS_30 <- draw_SRS_each(full,samp_nk[3])  ## 30%
  set.seed(1234+99*jj)
  sample_SRS_50 <- draw_SRS_each(full,samp_nk[4])  ## 50%
  
  ## SRS proportional to Nk
  set.seed(1234+99*jj) ## set same seed each time so that samples are cumulative
  sample_SRSpropto_10 <- draw_SRSpropto_each(full,samp_frac[1])  ## 10%
  set.seed(1234+99*jj)
  sample_SRSpropto_20 <- draw_SRSpropto_each(full,samp_frac[2])  ## 20%
  set.seed(1234+99*jj)
  sample_SRSpropto_30 <- draw_SRSpropto_each(full,samp_frac[3])  ## 30%
  set.seed(1234+99*jj)
  sample_SRSpropto_50 <- draw_SRSpropto_each(full,samp_frac[4])  ## 50%

  ## SRSmargin at the margin
  set.seed(1234+99*jj) ## set same seed each time so that samples are cumulative
  sample_SRSmargin_10 <- draw_SRSmargin_each(full,samp_nk[1])  ## 10%
  set.seed(1234+99*jj)
  sample_SRSmargin_20 <- draw_SRSmargin_each(full,samp_nk[2])  ## 20%
  set.seed(1234+99*jj)
  sample_SRSmargin_30 <- draw_SRSmargin_each(full,samp_nk[3])  ## 30%
  set.seed(1234+99*jj)
  sample_SRSmargin_50 <- draw_SRSmargin_each(full,samp_nk[4])  ## 50%
  
  ## SRSmargin at the margin (SRR)
  set.seed(1234+99*jj) ## set same seed each time so that samples are cumulative
  sample_SRSmarginSRR_10 <- draw_SRSmarginSRR_each(full,SRRnaive,samp_nk[1])  ## 10%
  set.seed(1234+99*jj)
  sample_SRSmarginSRR_20 <- draw_SRSmarginSRR_each(full,SRRnaive,samp_nk[2])  ## 20%
  set.seed(1234+99*jj)
  sample_SRSmarginSRR_30 <- draw_SRSmarginSRR_each(full,SRRnaive,samp_nk[3])  ## 30%
  set.seed(1234+99*jj)
  sample_SRSmarginSRR_50 <- draw_SRSmarginSRR_each(full,SRRnaive,samp_nk[4])  ## 50%
  
  ####################################################################
  #                         Full Data Analysis                       #
  #################################################################### 
  fullreg <- glmmML(event~xw.cont+xw.bin+xbetween,cluster=clinic,family=binomial,data=full)
  fullest <- c(fullreg$coefficients,fullreg$sigma)
  fullSE <- c(fullreg$coef.sd,fullreg$sigma.sd) 
  
  results_full <- c(nrow(full),fullest,fullSE)
  
  ####################################################################
  #                           Random Sample                          #
  ####################################################################    
  
  results_RS <- c(fit_ML(sample_RS_10,fixed_only=TRUE),
                  fit_ML(sample_RS_20,fixed_only=TRUE),
                  fit_ML(sample_RS_30,fixed_only=TRUE),
                  fit_ML(sample_RS_50,fixed_only=TRUE) )  
  
  ####################################################################
  #                         Cluster Sampling                         #
  ####################################################################    
  
  results_CS <- c(fit_ML(sample_CS_10,fixed_only=TRUE),
                  fit_ML(sample_CS_20,fixed_only=TRUE),
                  fit_ML(sample_CS_30,fixed_only=TRUE),
                  fit_ML(sample_CS_50,fixed_only=TRUE) ) 
 
  ####################################################################
  #                   Case Control Sample - WEIGHTED                 #
  ####################################################################    

  
  results_CC <- c(fit_CC(sample_CC_10,fixed_only=TRUE),
                  fit_CC(sample_CC_20,fixed_only=TRUE),
                  fit_CC(sample_CC_30,fixed_only=TRUE),
                  fit_CC(sample_CC_50,fixed_only=TRUE) ) 
  
  
   
  ####################################################################
  #                   Case Control Sample - OFFSET                   #
  ####################################################################    
  
  results_CC_offset <- c(fit_CC(sample_CC_10,offset_method=TRUE,fixed_only=TRUE),
                         fit_CC(sample_CC_20,offset_method=TRUE,fixed_only=TRUE),
                         fit_CC(sample_CC_30,offset_method=TRUE,fixed_only=TRUE),
                         fit_CC(sample_CC_50,offset_method=TRUE,fixed_only=TRUE) )  
  
  
 

  
  ####################################################################
  #                           SRS ML                                 #
  ####################################################################   
  
  results_SRS <- c(fit_ML(sample_SRS_10,fixed_only=TRUE),
                   fit_ML(sample_SRS_20,fixed_only=TRUE),
                   fit_ML(sample_SRS_30,fixed_only=TRUE),
                   fit_ML(sample_SRS_50,fixed_only=TRUE) )
  
  
  ####################################################################
  #                   SRS - propto nk ML                             #
  ####################################################################   
  
  results_SRSpropto <- c(fit_ML(sample_SRSpropto_10,fixed_only=TRUE),
                         fit_ML(sample_SRSpropto_20,fixed_only=TRUE),
                         fit_ML(sample_SRSpropto_30,fixed_only=TRUE),
                         fit_ML(sample_SRSpropto_50,fixed_only=TRUE) )
  
  ####################################################################
  #                   SRS - at the margin ML                         #
  ####################################################################   
  
  results_SRSmargin <- c(fit_ML(sample_SRSmargin_10,fixed_only=TRUE),
                         fit_ML(sample_SRSmargin_20,fixed_only=TRUE),
                         fit_ML(sample_SRSmargin_30,fixed_only=TRUE),
                         fit_ML(sample_SRSmargin_50,fixed_only=TRUE) )
  
  ####################################################################
  #                   SRS - at the margin (SRR) ML                   #
  ####################################################################   
  
  results_SRSmarginSRR <- c(fit_ML(sample_SRSmarginSRR_10,fixed_only=TRUE),
                            fit_ML(sample_SRSmarginSRR_20,fixed_only=TRUE),
                            fit_ML(sample_SRSmarginSRR_30,fixed_only=TRUE),
                            fit_ML(sample_SRSmarginSRR_50,fixed_only=TRUE) )
  
  ####################################################################
  #                           SRS WL                          #
  ####################################################################   
  
  results_SRS_WL <- c(fit_CSCC(sample_SRS_10,SRS=TRUE,fixed_only=TRUE),
                           fit_CSCC(sample_SRS_20,SRS=TRUE,fixed_only=TRUE),
                           fit_CSCC(sample_SRS_30,SRS=TRUE,fixed_only=TRUE),
                           fit_CSCC(sample_SRS_50,SRS=TRUE,fixed_only=TRUE) )
  
  
  ####################################################################
  #                   SRS- propto nk - WL                     #
  ####################################################################   
  
  
  results_SRSpropto_WL <- c(fit_CSCC(sample_SRSpropto_10,SRS=TRUE,fixed_only=TRUE),
                                 fit_CSCC(sample_SRSpropto_20,SRS=TRUE,fixed_only=TRUE),
                                 fit_CSCC(sample_SRSpropto_30,SRS=TRUE,fixed_only=TRUE),
                                 fit_CSCC(sample_SRSpropto_50,SRS=TRUE,fixed_only=TRUE) )
  
  
  ####################################################################
  #                   SRS- at the margin - WL                 #
  ####################################################################   
  
  
  results_SRSmargin_WL <- c(fit_CSCC(sample_SRSmargin_10,SRS=TRUE,fixed_only=TRUE),
                                 fit_CSCC(sample_SRSmargin_20,SRS=TRUE,fixed_only=TRUE),
                                 fit_CSCC(sample_SRSmargin_30,SRS=TRUE,fixed_only=TRUE),
                                 fit_CSCC(sample_SRSmargin_50,SRS=TRUE,fixed_only=TRUE) )
  
  ####################################################################
  #                   SRS- at the margin (SRR) - WL           #
  ####################################################################   
  
  
  results_SRSmarginSRR_WL <- c(fit_CSCC(sample_SRSmarginSRR_10,SRS=TRUE,fixed_only=TRUE),
                                    fit_CSCC(sample_SRSmarginSRR_20,SRS=TRUE,fixed_only=TRUE),
                                    fit_CSCC(sample_SRSmarginSRR_30,SRS=TRUE,fixed_only=TRUE),
                                    fit_CSCC(sample_SRSmarginSRR_50,SRS=TRUE,fixed_only=TRUE) )
  
  
  ####################################################################
  #                           CSCC WEIGHTED                          #
  ####################################################################   
  
  results_CSCC <- c(fit_CSCC(sample_CSCC_10,fixed_only=TRUE),
                    fit_CSCC(sample_CSCC_20,fixed_only=TRUE),
                    fit_CSCC(sample_CSCC_30,fixed_only=TRUE),
                    fit_CSCC(sample_CSCC_50,fixed_only=TRUE) )
  
  
  ####################################################################
  #                   CSCC- propto nk - WEIGHTED                     #
  ####################################################################   
  
  
  results_CSCCpropto <- c(fit_CSCC(sample_CSCCpropto_10,fixed_only=TRUE),
                          fit_CSCC(sample_CSCCpropto_20,fixed_only=TRUE),
                          fit_CSCC(sample_CSCCpropto_30,fixed_only=TRUE),
                          fit_CSCC(sample_CSCCpropto_50,fixed_only=TRUE) )
  
  ####################################################################
  #                   CSCC- at the margin - WEIGHTED                 #
  ####################################################################   
  
  
  results_CSCCmargin <- c(fit_CSCC(sample_CSCCmargin_10,fixed_only=TRUE),
                          fit_CSCC(sample_CSCCmargin_20,fixed_only=TRUE),
                          fit_CSCC(sample_CSCCmargin_30,fixed_only=TRUE),
                          fit_CSCC(sample_CSCCmargin_50,fixed_only=TRUE) )
  
  ####################################################################
  #                   CSCC- at the margin (SRR) - WEIGHTED           #
  ####################################################################   
  
  
  results_CSCCmarginSRR <- c(fit_CSCC(sample_CSCCmarginSRR_10,fixed_only=TRUE),
                             fit_CSCC(sample_CSCCmarginSRR_20,fixed_only=TRUE),
                             fit_CSCC(sample_CSCCmarginSRR_30,fixed_only=TRUE),
                             fit_CSCC(sample_CSCCmarginSRR_50,fixed_only=TRUE) )
  
  
  
  
  ####################################################################
  #                           CSCC OFFSET                          #
  ####################################################################   
  
  results_CSCC_offset <- c(fit_CSCC(sample_CSCC_10,offset_method=TRUE,fixed_only=TRUE),
                           fit_CSCC(sample_CSCC_20,offset_method=TRUE,fixed_only=TRUE),
                           fit_CSCC(sample_CSCC_30,offset_method=TRUE,fixed_only=TRUE),
                           fit_CSCC(sample_CSCC_50,offset_method=TRUE,fixed_only=TRUE) )
  
  
  ####################################################################
  #                   CSCC- propto nk - OFFSET                     #
  ####################################################################   
  
  
  results_CSCCpropto_offset <- c(fit_CSCC(sample_CSCCpropto_10,offset_method=TRUE,fixed_only=TRUE),
                                 fit_CSCC(sample_CSCCpropto_20,offset_method=TRUE,fixed_only=TRUE),
                                 fit_CSCC(sample_CSCCpropto_30,offset_method=TRUE,fixed_only=TRUE),
                                 fit_CSCC(sample_CSCCpropto_50,offset_method=TRUE,fixed_only=TRUE) )
  
  
  ####################################################################
  #                   CSCC- at the margin - OFFSET                 #
  ####################################################################   
  
  
  results_CSCCmargin_offset <- c(fit_CSCC(sample_CSCCmargin_10,offset_method=TRUE,fixed_only=TRUE),
                                 fit_CSCC(sample_CSCCmargin_20,offset_method=TRUE,fixed_only=TRUE),
                                 fit_CSCC(sample_CSCCmargin_30,offset_method=TRUE,fixed_only=TRUE),
                                 fit_CSCC(sample_CSCCmargin_50,offset_method=TRUE,fixed_only=TRUE) )
  
  ####################################################################
  #                   CSCC- at the margin (SRR) - OFFSET           #
  ####################################################################   
  
  
  results_CSCCmarginSRR_offset <- c(fit_CSCC(sample_CSCCmarginSRR_10,offset_method=TRUE,fixed_only=TRUE),
                                    fit_CSCC(sample_CSCCmarginSRR_20,offset_method=TRUE,fixed_only=TRUE),
                                    fit_CSCC(sample_CSCCmarginSRR_30,offset_method=TRUE,fixed_only=TRUE),
                                    fit_CSCC(sample_CSCCmarginSRR_50,offset_method=TRUE,fixed_only=TRUE) )
  
  ####################################################################
  #                           Output Results                         #
  ####################################################################  
  ## report jj then for each formulation, report sample size, estimates, SEs, MSEP
  
  
  ####################################################################
  #                           Output Results                         #
  ####################################################################  
  ## report jj then for each formulation, report sample size, estimates, SEs, MSEP
  
  df <- rbind(df,c(jj,results_full,#results_naive,
                   results_RS,results_CS,
                   results_CC,results_CC_offset,
                   results_SRS,results_SRSpropto,results_SRSmargin,results_SRSmarginSRR,
                   results_SRS_WL,results_SRSpropto_WL,results_SRSmargin_WL,results_SRSmarginSRR_WL,
                   results_CSCC,results_CSCCpropto,results_CSCCmargin,results_CSCCmarginSRR,
                   results_CSCC_offset,results_CSCCpropto_offset,results_CSCCmargin_offset,results_CSCCmarginSRR_offset))
  
  },silent=TRUE)
  
}  

## report which iterations didn't converge 
print((1:num_iter)[!(1:num_iter) %in% df[,1]])

## exclude iteration number
df <- df[,-1]

## returns results
write.table(x=df, file="sim3_results.txt", row.names=F, col.names=F, quote=F)

