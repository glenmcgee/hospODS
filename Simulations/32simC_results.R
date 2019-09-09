####################################################################
#                                                                  #
#                  Simulation Study 3 for BASE Model               #
#                                                                  #
####################################################################
#                         Results Tables                            #
####################################################################
#                          March 2017                           #
####################################################################
library(xtable)
library(stargazer)
library(wesanderson)
library(ggplot2)
library(tidyr)

##### Colors for plots
# blue red green yellow
pal <- c(wes_palette(n=5, name="Zissou")[c(1,5)],wes_palette(n=5, name="Darjeeling")[c(2,3)])

wes_blue <- wes_palette(n=5, name="Zissou")[1]
wes_red <- wes_palette(n=5, name="Zissou")[5]
wes_green <- wes_palette(n=5, name="Darjeeling")[2]
wes_gold <- wes_palette(n=5, name="Darjeeling")[3]
wes_yellow <- wes_palette(n=5, name="Zissou")[3]
wes_orange <- wes_palette(n=5, name="Darjeeling")[4]
wes_navy <- wes_palette(n=5, name="Rushmore")[4]
wes_pink <- wes_palette(n=5, name="Royal2")[3]#wes_pink <- wes_palette(n=4, name="GrandBudapest2")[1]#wes_pink <- wes_palette(n=4, name="GrandBudapest")[2]
wes_darkgreen <- wes_palette(n=4, name="Chevalier")[1]

####################################################################
#                     Function to Make Tables                      #
####################################################################
sim3.table <- function(true_param,sim_results,mod_names,
                       var_names=c("Int","Within/cont","Within/bin","Between/cont","Sigma"),
                       latex=TRUE,mod_name="sim3 - Base"){
  n_it <- nrow(sim_results)
  lenpar <- length(true_param)
  nlen <- 2*lenpar+1 # *2 because SE and + 1 for sample sizes 
  nsizes <- ncol(sim_results)/nlen # number of sample size simulations (5= 10,20,30,50,Full)
  true_mat <- matrix(true_param,nrow=n_it,ncol=lenpar,byrow=TRUE) ## matrix of true values for coverage
  
  mod_names <- c("Full",mod_names) ## add full data
  
  ## get results
  get_results <- function(sim_results,i){
    temp_results <- sim_results[,nlen*i+ (1:nlen)]
    return(list(n_k=temp_results[,1],
                est=temp_results[,1+(1:lenpar)],
                SE=temp_results[,1+((lenpar+1):(2*lenpar))]))
  }
  
  
  mean_tab <- c()
  bias_tab <- c()
  pctbias_tab <- c()
  meanSE_tab <- c()
  SD_tab <- c()
  cvg_tab <- c()  
  SESD_tab <- c()
  MSE_tab <- c()

  
  
  for(ii in (0:(nsizes-1))){
    
    ## collect results
    results <- get_results(sim_results,ii)
    n_k <- results$n_k
    est <- results$est
    SE <- results$SE
 
    
    mean_tab <- rbind(mean_tab,c(mean(n_k),apply(est,2,mean,na.rm=TRUE)))
    bias_tab <- rbind(bias_tab,c(mean(n_k),apply(est-true_mat,2,mean,na.rm=TRUE)))
    pctbias_tab <- rbind(pctbias_tab,c(mean(n_k),100*apply((est-true_mat)/true_mat,2,mean,na.rm=TRUE)))
    meanSE_tab <- rbind(meanSE_tab,c(mean(n_k),apply(SE,2,mean,na.rm=TRUE)))
    SD_tab <- rbind(SD_tab,c(mean(n_k),apply(est,2,sd,na.rm=TRUE)))
    cvg_tab <- rbind(cvg_tab,c(mean(n_k),100*apply((est-1.96*SE < true_mat) & (true_mat < est+1.96*SE),2,mean,na.rm=TRUE)))
    SESD_tab <- rbind(SESD_tab,c(mean(n_k),apply(SE,2,mean,na.rm=TRUE)/apply(est,2,sd)))
    MSE_tab <- rbind(MSE_tab,c(mean(n_k),apply((est-true_mat)^2,2,mean,na.rm=TRUE)))
    
    
    
  }
  
  
  ## relative efficiency table:
  eff_tab <- SD_tab
  opt_SE <- matrix(rep(eff_tab[1,2:ncol(eff_tab)],nrow(eff_tab)),ncol=lenpar,byrow=TRUE)
  eff_tab[,2:ncol(eff_tab)] <- opt_SE/eff_tab[,2:ncol(eff_tab)]
  
  ## report true sampling fractions
  true_fracs <- mean_tab[-1,1]/mean_tab[1,1]
  
  rownames(mean_tab) <- mod_names; colnames(mean_tab) <- c('n',var_names)
  rownames(bias_tab) <- mod_names; colnames(bias_tab) <- c('n',var_names)
  rownames(pctbias_tab) <- mod_names; colnames(pctbias_tab) <- c('n',var_names)
  rownames(meanSE_tab) <- mod_names; colnames(meanSE_tab) <- c('n',var_names)
  rownames(cvg_tab) <- mod_names; colnames(cvg_tab) <- c('n',var_names)
  rownames(SESD_tab) <- mod_names; colnames(SESD_tab) <- c('n',var_names)
  rownames(eff_tab) <- mod_names;  colnames(eff_tab) <- c('n',var_names)
  rownames(MSE_tab) <- mod_names;  colnames(MSE_tab) <- c('n',var_names)
  
  if (latex==TRUE){
    return(list(xtable(mean_tab,caption=paste("Mean Estimates"," - ",mod_name),digits=2),
                xtable(bias_tab,caption=paste("Mean Bias"," - ",mod_name),digits=3),
                xtable(pctbias_tab,caption=paste("Mean Percent Bias"," - ",mod_name),digits=0),
                xtable(meanSE_tab,caption=paste("Mean SE"," - ",mod_name),digits=2),
                xtable(cvg_tab,caption=paste("Percent Coverage"," - ",mod_name),digits=0),
                xtable(SESD_tab,caption=paste("SE/SD Ratio"," - ",mod_name),digits=2),
                xtable(eff_tab,caption=paste("Relative Efficiency"," - ",mod_name),digits=2)))
    
  } else {
    return(list(mean=mean_tab,
                bias=bias_tab,
                pctbias=pctbias_tab,
                meanSE=meanSE_tab,
                cvg=cvg_tab,
                SESD=SESD_tab,
                eff=eff_tab,
                MSE=MSE_tab,
                true_fracs=true_fracs))
  }
  
  
}


####################################################################
#                       Function to Make Plots                     #
####################################################################
## Uses tables created above as input:


sim3.plot <- function(sim3table,title="",report_true_frac=FALSE,which_mod=seq(1,21)){
  ## collect data
  eff <- sim3table$eff
  MSE <- sim3table$MSE
  true_fracs <- sim3table$true_fracs
  

  
  ## repeat full data to be a flat line
  eff <- data.frame(rbind(eff[1,],eff[1,],eff[1,],eff))
  MSE <- data.frame(rbind(MSE[1,],MSE[1,],MSE[1,],MSE))
  
  ## add x variable (indexing sampling fraction)
  ## add grouping variable (for method)
  method <- rep(c("black",wes_green,wes_orange,wes_gold,wes_yellow,wes_red,wes_red,wes_red,wes_red,wes_pink,wes_pink,wes_pink,wes_pink,wes_blue,wes_blue,wes_blue,wes_blue,"grey","grey","grey","grey"),each=4)  
  ## full - black, RS - green, CS - orange, CC weighted gold, CC offset yellow, SRS red, SRS WL pink, CSCCWEIGHTED - blue, CSCC offset "grey"
  propto <- rep(c("solid","solid","solid","solid","solid","solid","dashed","dotted","dotdash","solid","dashed","dotted","dotdash","solid","dashed","dotted","dotdash","solid","dashed","dotted","dotdash"),each=4)  ##  propto <- rep(c(1,2,3,2,3),each=4) ## rep(c("Balanced","Balanced","Propto","Balanced","Propto"),each=4)
  id <- paste(method,propto)
  
  ## toggle true fractions or approximate fractions
  if (report_true_frac==TRUE){
    fracs <- c(c(0.10,0.20,0.30,0.50),true_fracs) 
  } else {
    fracs <- rep(c(0.10,0.20,0.30,0.50),21)
  }

  eff <- data.frame(cbind(id,method,propto,fracs,eff))
  MSE <- data.frame(cbind(id,method,propto,fracs,MSE))
  
  ## only show requested models in plots
  all_mod <- seq(1,nrow(eff)/4)
  which_row <- rep(all_mod %in% which_mod,each=4)
  eff <- eff[which_row,]
  MSE <- MSE[which_row,]
  
  ## axis limits for MSE plots:
  MSElim <- apply(MSE[6:10],2,max)   #old: [2:ncol(MSE)]
  
  var_names <- colnames(eff)[6:ncol(eff)]## actual titles
  #label_names <- c("(Intercept)","(a)","(b)","(c)","(d)")
  label_names <- c(expression(beta[0]),expression(beta[1]),expression(beta[2]),expression(beta[3]),expression(sigma))
  eff_plots <- list()
  MSE_plots <- list()
  for (p in (1:length(var_names))){
    temp_plot <-     ggplot(data=eff,aes_string(x="fracs",group="id",y=var_names[p]))+
      geom_line(aes(colour=as.factor(method),linetype=as.factor(propto)),size=1) +
      scale_colour_manual(values=setNames(method, method),guide=FALSE)+
      scale_linetype_manual(values=setNames(propto, propto),guide=FALSE)+
      theme_classic()+
      scale_y_continuous("Relative Efficiency",limits=c(0,1)) + 
      scale_x_continuous("Sampling Fraction")+
      ggtitle((label_names[p])) + 
      theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
      #ggtitle(paste(label_names[p],":",var_names[p])) #+ this would be to identify each one
      #ggtitle(paste(title,var_names[p])) #+ 
      #theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
    eff_plots[[p]] <- (temp_plot)
    
    temp_plot_MSE <-     ggplot(data=MSE,aes_string(x="fracs",group="id",y=var_names[p]))+
      geom_line(aes(colour=as.factor(method),linetype=as.factor(propto)),size=1) +
      scale_colour_manual(values=setNames(method, method),guide=FALSE)+
      scale_linetype_manual(values=setNames(propto, propto),guide=FALSE)+
      theme_classic()+
      scale_y_continuous("MSE",limits=c(0,MSElim[p])) + 
      scale_x_continuous("Sampling Fraction")+
      ggtitle((label_names[p])) + 
      theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
      #ggtitle(paste(label_names[p],":",var_names[p])) #+ this would be to identify each one
      #ggtitle(paste(title,var_names[p])) #+ 
      #ggtitle(paste(title,var_names[p])) #+ 
    #theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
    MSE_plots[[p]] <- (temp_plot_MSE)
    
    
  }
  
  
  return(list(eff_plots=eff_plots,MSE_plots=MSE_plots))
  
}







####################################################################
#                          Tables for sim3                         #
####################################################################
# ev.rate <- 0.2     ## Event rate
# sigma <- 0.2       ## sigma: SD of random effects
# bin.risk <- 0.2    ## % dual eligible
# R=1000, 
model_names <- paste(rep(c("RS","CS","CC (Weighted)","CC (Offset)",rep(c("SRS ","SRS (WL)", "CSCC (Weighted) ","CSCC (Offset) "),each=4)),each=4),
                     c(rep("",16),rep(rep(c("Bal ","Prop ","Rate ","SRR "),each=4),4)),
                     rep(c(10,20,30,50),16))


## 150 Clusters 
sim3_results_K150 <- read.table("sim3_results.txt",header=FALSE)
true_results_K150 <- c(-1.7407,0.2,0.2,0.5,0.2)

sim3table_K150 <- sim3.table(true_results_K150,sim3_results_K150,
                             mod_names=model_names,
                             var_names=c("Int","Within/cont","Within/bin","Between/cont","Sigma"),
                             latex=FALSE,mod_name="Sim3 - Base")

#######################################################
### Make relative uncertainty table (from efficiency table)
longtab_relunc <- data.frame(sim3table_K150$eff)
longtab_relunc <- longtab_relunc[-1,-1] ## exclude full 
methods <- c(rownames(longtab_relunc)) ## id for methods
methods <- c(substr(methods,1,nchar(methods)-2))
longtab_relunc$methods <- methods
longtab_relunc$size <- rep(c(10,20,30,50),nrow(longtab_relunc)/4) ## sampling rate (%)

## long to wide
tab0 <- longtab_relunc %>% select(methods,Int,size) %>% spread(key=methods,value=Int)
tab1 <- longtab_relunc %>% select(methods,Within.cont,size) %>% spread(key=methods,value=Within.cont)
tab2 <- longtab_relunc %>% select(methods,Within.bin,size) %>% spread(key=methods,value=Within.bin)
tab3 <- longtab_relunc %>% select(methods,Between.cont,size) %>% spread(key=methods,value=Between.cont)
tab4 <- longtab_relunc %>% select(methods,Sigma,size) %>% spread(key=methods,value=Sigma)
tab_relunc <- rbind(tab0,tab1,tab2,tab3,tab4)

tab_relunc <- cbind(rep(c("b0","b1","b2","b3","sig"),each=4),tab_relunc$size,100*tab_relunc[,c(13,2,4,14,15,9,10,5,6)]) ## select only relevant methods
colnames(tab_relunc)[1:2] <- c("Param","(%)")
#xtable(tab_relunc,digits=0)

#######################################################
### Table 3
stargazer(tab_relunc,digits=0)



####################################################################
#                          Plots for sim3                         #
####################################################################
# ev.rate <- 0.2     ## Event rate
# sigma <- 0.2       ## sigma: SD of random effects
# bin.risk <- 0.2    ## % dual eligible
# R=1000, 


## 150 Clusters 
## CSCC is Blue, SRS is red
## balanced is dashed, propto is dotted
sim3plots_K150 <- sim3.plot(sim3table_K150,which_mod=c(1,2,3,5,6:9,14:17))# excluding offset method and CC weighted and CSCC offset and SRS weighted

## Efficiency
sim3plots_K150$eff_plots[[1]]
ggsave("Sim3_effplot_1.pdf",width=4,height=4)
sim3plots_K150$eff_plots[[2]]
ggsave("Sim3_effplot_2.pdf",width=4,height=4)
sim3plots_K150$eff_plots[[3]]
ggsave("Sim3_effplot_3.pdf",width=4,height=4)
sim3plots_K150$eff_plots[[4]]
ggsave("Sim3_effplot_4.pdf",width=4,height=4)
sim3plots_K150$eff_plots[[5]]
ggsave("Sim3_effplot_5.pdf",width=4,height=4)


## MSE
sim3plots_K150$MSE_plots[[1]]
ggsave("Sim3_MSEplot_1.pdf",width=4,height=4)
sim3plots_K150$MSE_plots[[2]]
ggsave("Sim3_MSEplot_2.pdf",width=4,height=4)
sim3plots_K150$MSE_plots[[3]]
ggsave("Sim3_MSEplot_3.pdf",width=4,height=4)
sim3plots_K150$MSE_plots[[4]]
ggsave("Sim3_MSEplot_4.pdf",width=4,height=4)
sim3plots_K150$MSE_plots[[5]]
ggsave("Sim3_MSEplot_5.pdf",width=4,height=4)







