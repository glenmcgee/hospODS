####################################################################
#                                                                  #
#                  Simulation Study 1 for BASE Model               #
#                                                                  #
####################################################################
#                         Results Tables                            #
####################################################################
#                            March 2017                            #
####################################################################
library(xtable)
library(stargazer)

####################################################################
#                     Function to Make Tables                      #
####################################################################
sim1.table <- function(true_param,sim_results,mod_names,
                             var_names=c("Int","Within/cont","Within/bin","Between/cont","Sigma"),
                             latex=TRUE,mod_name="Sim1 - Base"){
  n_it <- nrow(sim_results)
  lenpar <- length(true_param)
  nlen <- 2*lenpar+1+1+1+1 # *2 because SE and +1 for MSPE +1 for SRR MSEP +1 for misclass +1 for sample sizes 
  nsizes <- ncol(sim_results)/nlen # number of sample size simulations (5= 10,20,30,50,Full)
  true_mat <- matrix(true_param,nrow=n_it,ncol=lenpar,byrow=TRUE) ## matrix of true values for coverage
  
  mod_names <- c("Full",mod_names) ## add full data
  
  ## get results
  get_results <- function(sim_results,i){
    temp_results <- sim_results[,nlen*i+ (1:nlen)]
    return(list(n_k=temp_results[,1],
                est=temp_results[,1+(1:lenpar)],
                SE=temp_results[,1+((lenpar+1):(2*lenpar))],
                MSEP=temp_results[,(2*lenpar+2):nlen]))
  }
  

  mean_tab <- c()
  bias_tab <- c()
  pctbias_tab <- c()
  meanSE_tab <- c()
  SD_tab <- c()
  cvg_tab <- c()  
  SESD_tab <- c()
  MSEP_tab <- c()
  
  
  for(ii in (0:(nsizes-1))){
    
    ## collect results
    results <- get_results(sim_results,ii)
    n_k <- results$n_k
    est <- results$est
    SE <- results$SE
    MSEP <- results$MSEP
    
    mean_tab <- rbind(mean_tab,c(mean(n_k),apply(est,2,mean,na.rm=TRUE)))
    bias_tab <- rbind(bias_tab,c(mean(n_k),apply(est-true_mat,2,mean,na.rm=TRUE)))
    pctbias_tab <- rbind(pctbias_tab,c(mean(n_k),100*apply((est-true_mat)/true_mat,2,mean,na.rm=TRUE)))
    meanSE_tab <- rbind(meanSE_tab,c(mean(n_k),apply(SE,2,mean,na.rm=TRUE)))
    SD_tab <- rbind(SD_tab,c(mean(n_k),apply(est,2,sd,na.rm=TRUE)))
    cvg_tab <- rbind(cvg_tab,c(mean(n_k),100*apply((est-1.96*SE < true_mat) & (true_mat < est+1.96*SE),2,mean,na.rm=TRUE)))
    SESD_tab <- rbind(SESD_tab,c(mean(n_k),apply(SE,2,mean,na.rm=TRUE)/apply(est,2,sd)))
    MSEP_tab <- rbind(MSEP_tab,c(mean(n_k),apply(MSEP,2,mean,na.rm=TRUE)))
    
  }
  
  ## relative efficiency table:
  eff_tab <- SD_tab
  opt_SE <- matrix(rep(eff_tab[1,2:ncol(eff_tab)],nrow(eff_tab)),ncol=lenpar,byrow=TRUE)
  eff_tab[,2:ncol(eff_tab)] <- opt_SE/eff_tab[,2:ncol(eff_tab)]
  
  
  rownames(mean_tab) <- mod_names; colnames(mean_tab) <- c('n',var_names)
  rownames(bias_tab) <- mod_names; colnames(bias_tab) <- c('n',var_names)
  rownames(pctbias_tab) <- mod_names; colnames(pctbias_tab) <- c('n',var_names)
  rownames(meanSE_tab) <- mod_names; colnames(meanSE_tab) <- c('n',var_names)
  rownames(cvg_tab) <- mod_names; colnames(cvg_tab) <- c('n',var_names)
  rownames(SESD_tab) <- mod_names; colnames(SESD_tab) <- c('n',var_names)
  rownames(MSEP_tab) <- mod_names;  colnames(MSEP_tab) <- c('n','MSEP','SRR','misclass') 
  rownames(eff_tab) <- mod_names;  colnames(eff_tab) <- c('n',var_names)
  
  
  if (latex==TRUE){
    return(list(xtable(mean_tab,caption=paste("Mean Estimates"," - ",mod_name),digits=2),
                xtable(bias_tab,caption=paste("Mean Bias"," - ",mod_name),digits=3),
                xtable(pctbias_tab,caption=paste("Mean Percent Bias"," - ",mod_name),digits=0),
                xtable(meanSE_tab,caption=paste("Mean SE"," - ",mod_name),digits=2),
                xtable(cvg_tab,caption=paste("Percent Coverage"," - ",mod_name),digits=0),
                xtable(SESD_tab,caption=paste("SE/SD Ratio"," - ",mod_name),digits=2),
                xtable(MSEP_tab,caption=paste("MSEP"," - ",mod_name),digits=4)))
    
  } else {
    return(list(mean_tab,bias_tab,pctbias_tab,meanSE_tab,cvg_tab,SESD_tab,MSEP_tab))
  }
  
  
}


####################################################################
#                          Tables for Sim1                         #
####################################################################
# ev.rate <- 0.2     ## Event rate
# sigma <- 0.2       ## sigma: SD of random effects
# bin.risk <- 0.2    ## % dual eligible
# R=5000, 2% didn't converge...
model_names <- paste(rep(c("CSCC (Weighted) ","CSCC (Offset) ","CSCC (Naive)"),each=4),rep(c(10,20,30,50),2))

## 150 Clusters (N=40751, n=(3000,5883,8544,13316)  )
sim1_results_K150 <- read.table("sim1_results_K150.txt",header=FALSE)
true_results_K150 <- c(-1.7407,0.2,0.2,0.5,0.2)
sizes_K150 <- c(3000,5883,8544,13316,40751)
sim1table_K150 <- sim1.table(true_results_K150,sim1_results_K150,
                             mod_names=model_names,
                        var_names=c("Int","Within/cont","Within/bin","Between/cont","Sigma"),
                        latex=FALSE,mod_name="Sim1 - Base")

## 200 Clusters (N=53913, n=4000,7881,11480,17877)  )
sim1_results_K200 <- read.table("sim1_results_K200.txt",header=FALSE)
true_results_K200 <- c(-1.7397,0.2,0.2,0.5,0.2)
sizes_K200 <- c(4000,7881,11480,17877,53913)
sim1table_K200 <- sim1.table(true_results_K200,sim1_results_K200,
                             mod_names=model_names,
                             var_names=c("Int","Within/cont","Within/bin","Between/cont","Sigma"),
                             latex=FALSE,mod_name="Sim1 - Base")



## 300 Clusters (N=83523, n=(6000,11844,17298,26949)  )
sim1_results_K300 <- read.table("sim1_results_K300.txt",header=FALSE)
true_results_K300 <- c(-1.7368,0.2,0.2,0.5,0.2)
sizes_K300 <- c(6000,11844,17298,26949,83523)
sim1table_K300 <-  sim1.table(true_results_K300,sim1_results_K300,
                              mod_names=model_names,
                              var_names=c("Int","Within/cont","Within/bin","Between/cont","Sigma"),
                              latex=FALSE,mod_name="Sim1 - Base")


# ## Combined tables
# ## weighted
# for (tt in (1:7)){
#   tab <- rbind((sim1table_K150[[tt]])[c(2:5,1),], ## reorder rows of each one so full comes after samples
#                (sim1table_K200[[tt]])[c(2:5,1),],
#                (sim1table_K300[[tt]])[c(2:5,1),])
#   stargazer(tab,digits=c(3,3,0,2,0,2,4)[tt])
#   
# }
# 
# ## offset
# for (tt in (1:7)){
#   tab <- rbind((sim1table_K150[[tt]])[c(6:9,1),], ## reorder rows of each one so full comes after samples
#                (sim1table_K200[[tt]])[c(6:9,1),],
#                (sim1table_K300[[tt]])[c(6:9,1),])
#   stargazer(tab,digits=c(3,3,0,2,0,2,4)[tt])
#   
# }
# 
# ## naive
# for (tt in (1:7)){
#   tab <- rbind((sim1table_K150[[tt]])[c(10:13,1),], ## reorder rows of each one so full comes after samples
#                (sim1table_K200[[tt]])[c(10:13,1),],
#                (sim1table_K300[[tt]])[c(10:13,1),])
#   stargazer(tab,digits=c(3,3,0,2,0,2,4)[tt])
#   
# }


#####################################
## PRINT TABLE 2:
tab_est <- round(rbind((sim1table_K150[[1]])[c(2:5,1),], ## reorder rows of each one so full comes after samples
                       (sim1table_K200[[1]])[c(2:5,1),],
                       (sim1table_K300[[1]])[c(2:5,1),]),2)

tab_cvg <- round(rbind((sim1table_K150[[5]])[c(2:5,1),], ## reorder rows of each one so full comes after samples
                       (sim1table_K200[[5]])[c(2:5,1),],
                       (sim1table_K300[[5]])[c(2:5,1),]),0)
tab1 <- cbind(tab_est[,1:2],tab_cvg[,2],tab_est[,3],tab_cvg[,3],tab_est[,4],tab_cvg[,4],tab_est[,5],tab_cvg[,5],tab_est[,6],tab_cvg[,6])
stargazer(tab1,digits=2)

## extra TABLE OFFSET:
tab_est_off <- round(rbind((sim1table_K150[[1]])[c(6:9,1),], ## reorder rows of each one so full comes after samples
                       (sim1table_K200[[1]])[c(6:9,1),],
                       (sim1table_K300[[1]])[c(6:9,1),]),2)

tab_cvg_off <- round(rbind((sim1table_K150[[5]])[c(6:9,1),], ## reorder rows of each one so full comes after samples
                       (sim1table_K200[[5]])[c(6:9,1),],
                       (sim1table_K300[[5]])[c(6:9,1),]),0)
tab1_off <- cbind(tab_est_off[,1:2],tab_cvg_off[,2],tab_est_off[,3],tab_cvg_off[,3],tab_est_off[,4],tab_cvg_off[,4],tab_est_off[,5],tab_cvg_off[,5],tab_est_off[,6],tab_cvg_off[,6])
stargazer(tab1_off,digits=2)



