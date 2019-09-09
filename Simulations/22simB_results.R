####################################################################
#                                                                  #
#                  Simulation Study 2 for BASE Model               #
#                                                                  #
####################################################################
#                         Results Tables                           #
####################################################################
#                             March 2017                           #
####################################################################
library(xtable)
library(stargazer)
library(wesanderson)
library(ggplot2)

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
wes_pink <- wes_palette(n=5, name="Royal2")[3]
wes_darkgreen <- wes_palette(n=4, name="Chevalier")[1]


####################################################################
#                     Function to Make Tables                      #
####################################################################
sim2.table <- function(true_param,sim_results,mod_names,
                       var_names=c("Int","Within/cont","Within/bin","Between/cont","Sigma"),
                       latex=TRUE,mod_name="Sim2 - Base"){
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

  ## report true sampling fractions
  true_fracs <- mean_tab[-1,1]/mean_tab[1,1]
  
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
                xtable(MSEP_tab,caption=paste("MSEP"," - ",mod_name),digits=4),
                xtable(eff_tab,caption=paste("Relative Efficiency"," - ",mod_name),digits=2)))
    
  } else {
    return(list(mean=mean_tab,
                bias=bias_tab,
                pctbias=pctbias_tab,
                meanSE=meanSE_tab,
                cvg=cvg_tab,
                SESD=SESD_tab,
                MSEP=MSEP_tab,
                eff=eff_tab,
                true_fracs=true_fracs))
  }
  
  
}


####################################################################
#                       Function to Make Plots                     #
####################################################################
## Uses tables created above as input:


sim2.plot <- function(sim2table,title="",report_true_frac=FALSE,which_mod=seq(1,17)){
  ## collect data
  eff <- sim2table$eff
  MSEP <- sim2table$MSEP
  true_fracs <- sim2table$true_fracs
  
  ## repeat full data to be a flat line
  eff <- data.frame(rbind(eff[1,],eff[1,],eff[1,],eff))
  MSEP <- data.frame(rbind(MSEP[1,],MSEP[1,],MSEP[1,],MSEP))
  
  ## add x variable (indexing sampling fraction)
  ## add grouping variable (for method)
  method <- rep(c("black",wes_red,wes_red,wes_red,wes_red,wes_blue,wes_blue,wes_blue,wes_blue,wes_pink,wes_pink,wes_pink,wes_pink,"grey","grey","grey","grey"),each=4)  ## rep(c("Full","SRS","SRS","CSCC","CSCC"),each=4)
  propto <- rep(c("solid","solid","dashed","dotted","dotdash","solid","dashed","dotted","dotdash","solid","dashed","dotted","dotdash","solid","dashed","dotted","dotdash"),each=4)  ##  propto <- rep(c(1,2,3,2,3),each=4) ## rep(c("Balanced","Balanced","Propto","Balanced","Propto"),each=4)
  id <- paste(method,propto)
  
  ## toggle true fractions or approximate fractions
  if (report_true_frac==TRUE){
    fracs <- c(c(0.10,0.20,0.30,0.50),true_fracs) 
  } else {
    fracs <- rep(c(0.10,0.20,0.30,0.50),17)
  }
  
  eff <- data.frame(cbind(id,method,propto,fracs,eff))
  MSEP <- data.frame(cbind(id,method,propto,fracs,MSEP))  
  
  ## only show requested models in plots
  all_mod <- seq(1,nrow(eff)/4)
  which_row <- rep(all_mod %in% which_mod,each=4)
  eff <- eff[which_row,]
  MSEP <- MSEP[which_row,]
  
  var_names <- colnames(eff)[6:ncol(eff)] ## actual titles
  label_names <- c("(Intercept)","(a)","(b)","(c)","(d)")
  eff_plots <- list()
  for (p in (1:length(var_names))){
    temp_plot <-     ggplot(data=eff,aes_string(x="fracs",group="id",y=var_names[p]))+
                              geom_line(aes(colour=as.factor(method),linetype=as.factor(propto)),size=1) +
                              scale_colour_manual(values=setNames(method, method),guide=FALSE)+
                              scale_linetype_manual(values=setNames(propto, propto),guide=FALSE)+
                              theme_classic()+
                              scale_y_continuous("Relative Efficiency",limits=c(0,1)) + 
                              scale_x_continuous("Sampling Fraction")+
                              ggtitle(paste(label_names[p])) + 
                              theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
                              #ggtitle(paste(label_names[p],":",var_names[p])) #+ this would be to identify each one
                              #ggtitle(paste(title,var_names[p])) #+ 
                              #theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
    eff_plots[[p]] <- (temp_plot)
  }
    
  
  
  MSEP_plot <- ggplot(data=MSEP,aes(x=fracs,y=(MSEP),group=id,colour=as.factor(method),linetype=as.factor(propto)))+
    geom_line(size=1) +
    scale_colour_manual(values=setNames(method, method),guide=FALSE)+
    scale_linetype_manual(values=setNames(propto, propto),guide=FALSE)+
    theme_classic()+
    scale_y_continuous("MSEP",limits=c(0,max(MSEP$SRR,MSEP$MSEP))) + 
    scale_x_continuous("Sampling Fraction")+
    #ggtitle("(a)") + 
    ggtitle("Random Effects") + 
    theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
  
  SRR_plot <- ggplot(data=MSEP,aes(x=fracs,y=(SRR),group=id,colour=as.factor(method),linetype=as.factor(propto)))+
    geom_line(size=1) +
    scale_colour_manual(values=setNames(method, method),guide=FALSE)+
    scale_linetype_manual(values=setNames(propto, propto),guide=FALSE)+
    theme_classic()+
    scale_y_continuous("MSEP of RASRRs",limits=c(0,max(MSEP$SRR,MSEP$MSEP))) + 
    scale_x_continuous("Sampling Fraction")+
    # ggtitle("(b)") + 
    ggtitle("RASRRs") + 
    theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
  
  misclass_plot <- ggplot(data=MSEP,aes(x=fracs,y=(misclass),group=id,colour=as.factor(method),linetype=as.factor(propto)))+
    geom_line(size=1) +
    scale_colour_manual(values=setNames(method, method),guide=FALSE)+
    scale_linetype_manual(values=setNames(propto, propto),guide=FALSE)+
    theme_classic()+
    scale_y_continuous("No. Misclassified",limits=c(0,max(MSEP$misclass)))+#
    scale_x_continuous("Sampling Fraction")+
    # ggtitle(title) + 
    ggtitle("Misclassified Hospitals") + 
    theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
  
  
  return(list(eff_plots=eff_plots,MSEP_plot=MSEP_plot,SRR_plot=SRR_plot,misclass_plot=misclass_plot))
  
}



 



####################################################################
#                          Tables for Sim2                         #
####################################################################
# ev.rate <- 0.2     ## Event rate
# sigma <- 0.2       ## sigma: SD of random effects
# bin.risk <- 0.2    ## % dual eligible
# R=1000, 
## SRS, CSCC (Weighted and Offset)
## designs: balanced, proportional to cluster size, oversampling rates on the margin, oversampling SRRs on the margin
model_names <- paste(rep(c("SRS ","CSCC (Weighted) ","SRS (Weighted) ","CSCC (Offset) "),each=16),
                     rep(rep(c("Bal ","Prop ","Rate ","SRR "),each=4),3),
                     rep(c(10,20,30,50),12))


## 150 Clusters 
## this is the final one with proper SRRs. old ones are K150,K200,K300
sim2_results_K150 <- read.table("sim2_results.txt",header=FALSE) 
true_results_K150 <- c(-1.7407,0.2,0.2,0.5,0.2)



sim2table_K150 <- sim2.table(true_results_K150,sim2_results_K150,
                             mod_names =model_names,
                             var_names=c("Int","Within/cont","Within/bin","Between/cont","Sigma"),
                             latex=FALSE,mod_name="Sim2 - Base")





####################################################################
#                          Plots for Sim2                         #
####################################################################
# ev.rate <- 0.2     ## Event rate
# sigma <- 0.2       ## sigma: SD of random effects
# bin.risk <- 0.2    ## % dual eligible
# R=1000, 


## 150 Clusters 
## CSCC is Blue, SRS is red
## balanced is dashed, propto is dotted
sim2plots_K150 <- sim2.plot(sim2table_K150,which_mod=seq(1,9)) # for ICHPS ie no offset 1:13 to get SRS WL

sim2plots_K150$MSEP_plot
ggsave("sim2_MSEP_plot.pdf",width=4,height=4)
sim2plots_K150$SRR_plot
ggsave("sim2_MSEP_SRR_plot.pdf",width=4,height=4)
sim2plots_K150$misclass_plot
ggsave("sim2_misclass_plot.pdf",width=4,height=4)




## 150 Clusters 
## CSCC is Blue, SRS is red
## balanced is dashed, propto is dotted
## with offset 
sim2plots_K150_offset <- sim2.plot(sim2table_K150,which_mod=c(1:9,14:17)) # with offset

sim2plots_K150_offset$MSEP_plot
ggsave("sim2_MSEP_plot_offset.pdf",width=4,height=4)
sim2plots_K150_offset$SRR_plot
ggsave("sim2_MSEP_SRR_plot_offset.pdf",width=4,height=4)
sim2plots_K150_offset$misclass_plot
ggsave("sim2_misclass_plot_offset.pdf",width=4,height=4)




## 150 Clusters 
## CSCC is Blue, SRS is red
## balanced is dashed, propto is dotted
## with offset but no rank-based variants
sim2plots_K150_offset <- sim2.plot(sim2table_K150,which_mod=c(1:3,6:7,14:15)) 

sim2plots_K150_offset$MSEP_plot
ggsave("sim2_MSEP_plot_offset_norank.pdf",width=4,height=4)
sim2plots_K150_offset$SRR_plot
ggsave("sim2_MSEP_SRR_plot_offset_norank.pdf",width=4,height=4)
sim2plots_K150_offset$misclass_plot
ggsave("sim2_misclass_plot_offset_norank.pdf",width=4,height=4)
