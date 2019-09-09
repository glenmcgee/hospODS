
####################################################################
#                                                                  #
#                  Simulation Study for BASE Model                 #
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
set.seed(1234)
samp_sizes <- sample(samp_sizes,Kclust)  ## draw relevant sample sizes
N <- sum(samp_sizes)
ev.rate <- 0.2     ## Event rate
sigma <- 0.2    ## sigma: SD of random effects
bin.risk <- 0.2 ## % dual eligible


## set seed for reproducibility
set.seed(1000)
num_iter <- 5000

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

  
  
  ####################################################################
  #                            Draw Samples                          #
  ####################################################################
  set.seed(1234+99*jj) ## set same seed each time so that samples are cumulative
  sample_CSCC_10 <- draw_CSCC_each(full,20)  ## n1k=10
  set.seed(1234+99*jj)
  sample_CSCC_20 <- draw_CSCC_each(full,40)  ## n1k=20
  set.seed(1234+99*jj)
  sample_CSCC_30 <- draw_CSCC_each(full,60)  ## n1k=30
  set.seed(1234+99*jj)
  sample_CSCC_50 <- draw_CSCC_each(full,100)  ## n1k=50
  
  
  ####################################################################
  #                         Full Data Analysis                       #
  #################################################################### 
  fullreg <- glmmML(event~xw.cont+xw.bin+xbetween,cluster=clinic,family=binomial,data=full)
  fullest <- c(fullreg$coefficients,fullreg$sigma)
  fullSE <- c(fullreg$coef.sd,fullreg$sigma.sd) 
  predfull <- fullreg$posterior.modes
  SRRfull <- compute_SRR(fullreg$coefficients,fullreg$sigma,predfull,full,'event~xw.cont+xw.bin+xbetween','clinic')
  MSPEfull <- mean((predfull-true_rand)^2)
  SRR_MSEPfull <- mean((SRRfull-true_SRR)^2)
  misclassfull <- sum((sign(SRRfull-1)!=sign(true_SRR-1)))
  
  results_full <- c(nrow(full),fullest,fullSE,MSPEfull,SRR_MSEPfull,misclassfull)

  ####################################################################
  #                           CSCC NAIVE                             #
  ####################################################################  
  results_naive <- c(fit_CSCC(sample_CSCC_10,naive_method=TRUE),#SCALED:,scaled=TRUE,fullsizes=full %>% group_by(clinic) %>% summarize(count=n())),
                    fit_CSCC(sample_CSCC_20,naive_method=TRUE),
                    fit_CSCC(sample_CSCC_30,naive_method=TRUE),
                    fit_CSCC(sample_CSCC_50,naive_method=TRUE))
  
  
  ####################################################################
  #                           CSCC WEIGHTED                          #
  ####################################################################   

  results_CSCC <- c(fit_CSCC(sample_CSCC_10),
                    fit_CSCC(sample_CSCC_20),
                    fit_CSCC(sample_CSCC_30),
                    fit_CSCC(sample_CSCC_50))
  
  ####################################################################
  #                            CSCC OFFSET                           #
  ####################################################################   
  
  results_CSCC_offset <- c(fit_CSCC(sample_CSCC_10,offset_method=TRUE),
                    fit_CSCC(sample_CSCC_20,offset_method=TRUE),
                    fit_CSCC(sample_CSCC_30,offset_method=TRUE),
                    fit_CSCC(sample_CSCC_50,offset_method=TRUE) )
  
  
  ####################################################################
  #                           Output Results                         #
  ####################################################################  

  df <- rbind(df,c(jj,results_full,results_CSCC,results_CSCC_offset,results_naive))
  
  },silent=TRUE)
 
}   

## report which iterations didn't converge 
print((1:num_iter)[!(1:num_iter) %in% df[,1]])

## exclude iteration number
df <- df[,-1]

## returns results
write.table(x=df, file="sim1_results_K150.txt", row.names=F, col.names=F, quote=F)

