####################################################################
#                                                                  #
#                 Output for Data Description Section              #
#                     Produce Regression Table                     #
#                           Produce Plots                          #
#                                                                  #
####################################################################
#                              Dec 2017                            #
####################################################################
library(stargazer)
library(xtable)
library(wesanderson)
library(RColorBrewer)
library(gridExtra)


## Set TRUE for plots of analyses with smaller sample size (n1k=20)
small20 <- FALSE

## results for both SES and COMORB versions are reported 
## for model checking plots: only COMORB is used
## SES profile plots are named SES_....
## COMORB versions of these are named PROF_...

####################################################################
#                            Load Data                             #
####################################################################
#source("process_data.R")


suffix <- ""
if (small20==TRUE){
  suffix <- "_20"
}

prefix <- ""

weights <- read.table("weights_results.txt",header=TRUE)
fixed_eff <- read.table(paste("fixed_eff_results",".txt",sep=""),header=TRUE)
pred <- read.table(paste("pred_results",".txt",sep=""),header=TRUE)
SRR <- read.table(paste("SRR_results",".txt",sep=""),header=TRUE)









####################################################################
#                         Colour Palettes                          #
####################################################################
# blue red green yellow
pal <- c(wes_palette(n=5, name="Zissou")[c(1,5)],wes_palette(n=5, name="Darjeeling")[c(2,3)])
# red to blue gradient
# very dark to very dark so differences are very apparent
colfunc <- colorRampPalette(c(brewer.pal(9,"Reds")[seq(9,3)],"#DBDBDB",brewer.pal(9,"Blues")[seq(3,9)]))
# light to light for SES plots
colfunc2 <- colorRampPalette(c(pal[2],pal[1]))

####################################################################
#                Create Dataframe of bks for Plotting              #
####################################################################
pred_gold <- pred$pred_gold; pred_WO <- pred$pred_WO;  pred_WOCOMORB <- pred$pred_WOCOMORB; 
if(small20==TRUE){
  pred_wglmm <- pred$pred_wglmm_20; pred_SRS <- pred$pred_SRS_20
  pred_off <- pred$pred_off_20
  pred_wglmm_prop <- pred$pred_wglmm_prop_20; pred_SRS_prop <- pred$pred_SRS_prop_20;
  pred_off_prop <- pred$pred_off_prop_20
}else{
  pred_wglmm <- pred$pred_wglmm; pred_SRS <- pred$pred_SRS; 
  pred_off <- pred$pred_off
  pred_wglmm_prop <- pred$pred_wglmm_prop; pred_SRS_prop <- pred$pred_SRS_prop; 
  pred_off_prop <- pred$pred_off_prop
}
pred <- data.frame(cbind(pred_gold,pred_WO,pred_WOCOMORB,pred_wglmm,pred_SRS,pred_off,pred_wglmm_prop,pred_SRS_prop,pred_off_prop))

## define ranks
pred$rank_gold <- rank(pred_gold)
pred$rank_WO<- rank(pred_WO)
pred$rank_WOCOMORB<- rank(pred_WOCOMORB)
pred$rank_wglmm <- rank(pred_wglmm)
pred$rank_SRS <- rank(pred_SRS)
pred$rank_wglmm_prop <- rank(pred_wglmm_prop)
pred$rank_SRS_prop <- rank(pred_SRS_prop)

## define changes in rank
pred$rankdif_WO_gold <- pred$rank_WO-pred$rank_gold
pred$rankdif_WOCOMORB_gold <- pred$rank_WOCOMORB-pred$rank_gold
pred$rankdif_wglmm_gold <- pred$rank_wglmm-pred$rank_gold
pred$rankdif_SRS_gold <- pred$rank_SRS-pred$rank_gold
pred$rankdif_wglmm_prop_gold <- pred$rank_wglmm_prop-pred$rank_gold
pred$rankdif_SRS_prop_gold <- pred$rank_SRS_prop-pred$rank_gold

## add hospital level covariates
pred <- data.frame(cbind(pred,hosp))

## define colours based on category (>/<= 1)
pred$colCat_gold <- rep("lightgray",length(pred_gold)); pred$colCat_gold[pred_gold<=1] <- pal[1];  pred$colCat_gold[pred_gold>1] <- pal[2]
pred$colCat_WO <- rep("lightgray",length(pred_WO)); pred$colCat_WO[pred_WO<=1] <- pal[1];  pred$colCat_WO[pred_WO>1] <- pal[2]
pred$colCat_WOCOMORB <- rep("lightgray",length(pred_WOCOMORB)); pred$colCat_WOCOMORB[pred_WOCOMORB<=1] <- pal[1];  pred$colCat_WOCOMORB[pred_WOCOMORB>1] <- pal[2]
pred$colCat_wglmm <- rep("lightgray",length(pred_wglmm)); pred$colCat_wglmm[pred_wglmm<=1] <- pal[1];  pred$colCat_wglmm[pred_wglmm>1] <- pal[2]
pred$colCat_SRS <- rep("lightgray",length(pred_SRS)); pred$colCat_SRS[pred_SRS<=1] <- pal[1];  pred$colCat_SRS[pred_SRS>1] <- pal[2]
pred$colCat_wglmm_prop <- rep("lightgray",length(pred_wglmm_prop)); pred$colCat_wglmm_prop[pred_wglmm_prop<=1] <- pal[1];  pred$colCat_wglmm_prop[pred_wglmm_prop>1] <- pal[2]
pred$colCat_SRS_prop <- rep("lightgray",length(pred_SRS_prop)); pred$colCat_SRS_prop[pred_SRS_prop<=1] <- pal[1];  pred$colCat_SRS_prop[pred_SRS_prop>1] <- pal[2]


## define colours based on category change with respect to gold
pred$colChange_WO <- rep("lightgray",length(pred_WO)); pred$colChange_WO[pred_WO>1 & pred_gold<=1] <- pal[1];  pred$colChange_WO[pred_WO<=1 & pred_gold>1] <- pal[2]
pred$colChange_WOCOMORB <- rep("lightgray",length(pred_WOCOMORB)); pred$colChange_WOCOMORB[pred_WOCOMORB>1 & pred_gold<=1] <- pal[1];  pred$colChange_WOCOMORB[pred_WOCOMORB<=1 & pred_gold>1] <- pal[2]
pred$colChange_wglmm <- rep("lightgray",length(pred_wglmm)); pred$colChange_wglmm[pred_wglmm>1 & pred_gold<=1] <- pal[1];  pred$colChange_wglmm[pred_wglmm<=1 & pred_gold>1] <- pal[2]
pred$colChange_SRS <- rep("lightgray",length(pred_SRS)); pred$colChange_SRS[pred_SRS>1 & pred_gold<=1] <- pal[1];  pred$colChange_SRS[pred_SRS<=1 & pred_gold>1] <- pal[2]
pred$colChange_wglmm_prop <- rep("lightgray",length(pred_wglmm_prop)); pred$colChange_wglmm_prop[pred_wglmm_prop>1 & pred_gold<=1] <- pal[1];  pred$colChange_wglmm_prop[pred_wglmm_prop<=1 & pred_gold>1] <- pal[2]
pred$colChange_SRS_prop <- rep("lightgray",length(pred_SRS_prop)); pred$colChange_SRS_prop[pred_SRS_prop>1 & pred_gold<=1] <- pal[1];  pred$colChange_SRS_prop[pred_SRS_prop<=1 & pred_gold>1] <- pal[2]


## define colours based on rank change with respect to gold
pred$colRankChange_WO <- colfunc2(length(pred_WO))[rank(pred$rankdif_WO_gold)]
pred$colRankChange_WOCOMORB <- colfunc2(length(pred_WOCOMORB))[rank(pred$rankdif_WOCOMORB_gold)]
pred$colRankChange_wglmm <- colfunc2(length(pred_wglmm))[rank(pred$rankdif_wglmm_gold)]
pred$colRankChange_SRS <- colfunc2(length(pred_SRS))[rank(pred$rankdif_SRS_gold)]
pred$colRankChange_wglmm_prop <- colfunc2(length(pred_wglmm_prop))[rank(pred$rankdif_wglmm_prop_gold)]
pred$colRankChange_SRS_prop <- colfunc2(length(pred_SRS_prop))[rank(pred$rankdif_SRS_prop_gold)]


## define colours based on change with respect to gold
pred$colAbsChange_WO <- colfunc(length(pred_WO))[rank(pred$pred_WO-pred$pred_gold)]
pred$colAbsChange_WOCOMORB <- colfunc(length(pred_WOCOMORB))[rank(pred$pred_WOCOMORB-pred$pred_gold)]
pred$colAbsChange_wglmm <- colfunc(length(pred_wglmm))[rank(pred$pred_wglmm-pred$pred_gold)]
pred$colAbsChange_SRS <- colfunc(length(pred_SRS))[rank(pred$pred_SRS-pred$pred_gold)]
pred$colAbsChange_wglmm_prop <- colfunc(length(pred_wglmm_prop))[rank(pred$pred_wglmm_prop-pred$pred_gold)]
pred$colAbsChange_SRS_prop <- colfunc(length(pred_SRS_prop))[rank(pred$pred_SRS_prop-pred$pred_gold)]


## define axis limits
minlim_pred <- -max(abs(pred[,1:5])) # minlim <- min(pred[,1:3])
maxlim_pred <- max(abs(pred[,1:5])) # maxlim <- max(pred[,1:3])

####################################################################
#                          bk-Based Plots                          #
####################################################################


###############################################
## Comparing WO to Gold 
ggplot(data=pred) + 
  geom_hline(aes(yintercept=0),alpha=0.5) + 
  geom_vline(aes(xintercept=0),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(x=pred_WO,y=pred_gold),colour=pred$colAbsChange_WO,alpha=0.6,size=1) + 
  theme_classic()+
  scale_y_continuous("Full With SES",limits=c(minlim_pred,maxlim_pred)) + 
  scale_x_continuous("Full Without SES",limits=c(minlim_pred,maxlim_pred)) + 
  ggtitle("(a)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"pred_gold_WO",suffix,".pdf",sep=""),width=4,height=4)

###############################################
## Comparing WOCOMORB to Gold 
ggplot(data=pred) + 
  geom_hline(aes(yintercept=0),alpha=0.5) + 
  geom_vline(aes(xintercept=0),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(x=pred_WOCOMORB,y=pred_gold),colour=pred$colAbsChange_WOCOMORB,alpha=0.6,size=1) + 
  theme_classic()+
  scale_y_continuous("Full With Comorb",limits=c(minlim_pred,maxlim_pred)) + 
  scale_x_continuous("Full Without Comorb",limits=c(minlim_pred,maxlim_pred)) + 
  ggtitle("(a)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"pred_gold_WOCOMORB",suffix,".pdf",sep=""),width=4,height=4)


###############################################
## Comparing wglmm to Gold 
ggplot(pred) + 
  geom_hline(aes(yintercept=0),alpha=0.5) + 
  geom_vline(aes(xintercept=0),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(x=pred_wglmm,y=pred_gold),colour=pred$colAbsChange_wglmm,alpha=0.6,size=1) +
  theme_classic()+
  scale_y_continuous("Full With SES",limits=c(minlim_pred,maxlim_pred)) + 
  scale_x_continuous("CSCC",limits=c(minlim_pred,maxlim_pred)) +
  ggtitle("(b)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"pred_gold_wglmm",suffix,".pdf",sep=""),width=4,height=4)

###############################################
## Comparing SRS to Gold 
ggplot(pred) + 
  geom_hline(aes(yintercept=0),alpha=0.5) + 
  geom_vline(aes(xintercept=0),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(x=pred_SRS,y=pred_gold),colour=pred$colAbsChange_SRS,alpha=0.6,size=1) +
  theme_classic()+
  scale_y_continuous("Full With SES",limits=c(minlim_pred,maxlim_pred)) + 
  scale_x_continuous("SRS",limits=c(minlim_pred,maxlim_pred)) +
  ggtitle("(c)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"pred_gold_SRS",suffix,".pdf",sep=""),width=4,height=4)


###############################################
## Comparing wglmm to Gold   (prop)
ggplot(pred) + 
  geom_hline(aes(yintercept=0),alpha=0.5) + 
  geom_vline(aes(xintercept=0),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(x=pred_wglmm_prop,y=pred_gold),colour=pred$colAbsChange_wglmm_prop,alpha=0.6,size=1) +
  theme_classic()+
  scale_y_continuous("Full With SES",limits=c(minlim_pred,maxlim_pred)) + 
  scale_x_continuous("CSCC",limits=c(minlim_pred,maxlim_pred)) +
  ggtitle("(b)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"pred_gold_wglmm_prop",suffix,".pdf",sep=""),width=4,height=4)

###############################################
## Comparing SRS to Gold  (prop)
ggplot(pred) + 
  geom_hline(aes(yintercept=0),alpha=0.5) + 
  geom_vline(aes(xintercept=0),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(x=pred_SRS_prop,y=pred_gold),colour=pred$colAbsChange_SRS_prop,alpha=0.6,size=1) +
  theme_classic()+
  scale_y_continuous("Full With SES",limits=c(minlim_pred,maxlim_pred)) + 
  scale_x_continuous("SRS",limits=c(minlim_pred,maxlim_pred)) +
  ggtitle("(c)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"pred_gold_SRS_prop",suffix,".pdf",sep=""),width=4,height=4)



####################################
## ordered plots
pred.ordered <- pred[(order(pred_gold)),]
pred_order <- seq(1,length(pred_gold)) 
pred.ordered <- data.frame(cbind(pred.ordered,pred_order))

## WO vs GOLD
ggplot(pred.ordered) + 
  geom_point( aes(x=pred_order,y=pred_WO),alpha=0.8,colour=pal[3],size=2) + 
  geom_point( aes(x=pred_order,y=pred_gold),alpha=0.8, colour=pal[4],size=2) +   
  #  ggtitle("Standardized Readmission Rates (WO vs Gold)") +
  theme_classic()+
  scale_x_continuous("Rank") + scale_y_continuous("Random Effect",limits=c(minlim_pred,maxlim_pred)) +
  ggtitle("(a)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"ordered_bk_gold_WO",suffix,".pdf",sep=""),width=4,height=4)

## WOCOMORB vs GOLD
ggplot(pred.ordered) + 
  geom_point( aes(x=pred_order,y=pred_WOCOMORB),alpha=0.8,colour=pal[3],size=2) + 
  geom_point( aes(x=pred_order,y=pred_gold),alpha=0.8, colour=pal[4],size=2) +   
  #  ggtitle("Standardized Readmission Rates (WO vs Gold)") +
  theme_classic()+
  scale_x_continuous("Rank") + scale_y_continuous("Random Effect",limits=c(minlim_pred,maxlim_pred)) +
  ggtitle("(a)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"ordered_bk_gold_WOCOMORB",suffix,".pdf",sep=""),width=4,height=4)


## wglmm vs GOLD
ggplot(pred.ordered) + 
  geom_point( aes(x=pred_order,y=pred_wglmm),alpha=0.8,colour=pal[3],size=2) + 
  geom_point( aes(x=pred_order,y=pred_gold),alpha=0.8, colour=pal[4],size=2) +   
  #  ggtitle("Standardized Readmission Rates (wglmm vs Gold)") +
  theme_classic()+
  scale_x_continuous("Rank") + scale_y_continuous("Random Effect",limits=c(minlim_pred,maxlim_pred)) +
  ggtitle("(b)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"ordered_bk_gold_wglmm",suffix,".pdf",sep=""),width=4,height=4)

## SRS vs GOLD
ggplot(pred.ordered) + 
  geom_point( aes(x=pred_order,y=pred_SRS),alpha=0.8,colour=pal[3],size=2) + 
  geom_point( aes(x=pred_order,y=pred_gold),alpha=0.8, colour=pal[4],size=2) +   
  #  ggtitle("Standardized Readmission Rates (SRS vs Gold)") +
  theme_classic()+
  scale_x_continuous("Rank") + scale_y_continuous("Random Effect",limits=c(minlim_pred,maxlim_pred)) +
  ggtitle("(c)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"ordered_bk_gold_SRS",suffix,".pdf",sep=""),width=4,height=4)

## wglmm prop vs GOLD
ggplot(pred.ordered) + 
  geom_point( aes(x=pred_order,y=pred_wglmm_prop),alpha=0.8,colour=pal[3],size=2) + 
  geom_point( aes(x=pred_order,y=pred_gold),alpha=0.8, colour=pal[4],size=2) +   
  #  ggtitle("Standardized Readmission Rates (wglmm vs Gold)") +
  theme_classic()+
  scale_x_continuous("Rank") + scale_y_continuous("Random Effect",limits=c(minlim_pred,maxlim_pred)) +
  ggtitle("(b)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"ordered_bk_gold_wglmm_prop",suffix,".pdf",sep=""),width=4,height=4)

## SRS prop vs GOLD
ggplot(pred.ordered) + 
  geom_point( aes(x=pred_order,y=pred_SRS_prop),alpha=0.8,colour=pal[3],size=2) + 
  geom_point( aes(x=pred_order,y=pred_gold),alpha=0.8, colour=pal[4],size=2) +   
  #  ggtitle("Standardized Readmission Rates (SRS vs Gold)") +
  theme_classic()+
  scale_x_continuous("Rank") + scale_y_continuous("Random Effect",limits=c(minlim_pred,maxlim_pred)) +
  ggtitle("(c)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"ordered_bk_gold_SRS_prop",suffix,".pdf",sep=""),width=4,height=4)

####################################################################
#                Create Dataframe of SRRs for Plotting             #
####################################################################
SRR_gold <- SRR$SRR_gold; SRR_WO <- SRR$SRR_WO; SRR_WOCOMORB <- SRR$SRR_WOCOMORB; 

if(small20==TRUE){
  SRR_wglmm <- SRR$SRR_wglmm_20; SRR_SRS<- SRR$SRR_SRS_20;
  SRR_off <- SRR$SRR_off_20
  SRR_wglmm_prop <- SRR$SRR_wglmm_prop_20; SRR_SRS_prop<- SRR$SRR_SRS_prop_20;
  SRR_off_prop <- SRR$SRR_off_prop_20
  
}else{
  SRR_wglmm <- SRR$SRR_wglmm; SRR_SRS<- SRR$SRR_SRS;
  SRR_off <- SRR$SRR_off
  SRR_wglmm_prop <- SRR$SRR_wglmm_prop; SRR_SRS_prop<- SRR$SRR_SRS_prop;
  SRR_off_prop <- SRR$SRR_off_prop
  
}
SRR <- data.frame(cbind(SRR_gold,SRR_WO,SRR_WOCOMORB,SRR_wglmm,SRR_SRS,SRR_off,SRR_wglmm_prop,SRR_SRS_prop,SRR_off_prop))

## define ranks
SRR$rank_gold <- rank(SRR_gold)
SRR$rank_WO<- rank(SRR_WO)
SRR$rank_WOCOMORB<- rank(SRR_WOCOMORB)
SRR$rank_wglmm <- rank(SRR_wglmm)
SRR$rank_SRS <- rank(SRR_SRS)
SRR$rank_off <- rank(SRR_off)
SRR$rank_wglmm_prop <- rank(SRR_wglmm_prop)
SRR$rank_SRS_prop <- rank(SRR_SRS_prop)
SRR$rank_off_prop <- rank(SRR_off_prop)


## define changes in rank
SRR$rankdif_WO_gold <- SRR$rank_WO-SRR$rank_gold
SRR$rankdif_WOCOMORB_gold <- SRR$rank_WOCOMORB-SRR$rank_gold
SRR$rankdif_wglmm_gold <- SRR$rank_wglmm-SRR$rank_gold
SRR$rankdif_SRS_gold <- SRR$rank_SRS-SRR$rank_gold
SRR$rankdif_off_gold <- SRR$rank_off-SRR$rank_gold
SRR$rankdif_wglmm_prop_gold <- SRR$rank_wglmm_prop-SRR$rank_gold
SRR$rankdif_SRS_prop_gold <- SRR$rank_SRS_prop-SRR$rank_gold
SRR$rankdif_off_prop_gold <- SRR$rank_off_prop-SRR$rank_gold


## add hospital level covariates
SRR <- data.frame(cbind(SRR,hosp))

## define colours based on category (>/<= 1)
SRR$colCat_gold <- rep("lightgray",length(SRR_gold)); SRR$colCat_gold[SRR_gold<=1] <- pal[1];  SRR$colCat_gold[SRR_gold>1] <- pal[2]
SRR$colCat_WO <- rep("lightgray",length(SRR_WO)); SRR$colCat_WO[SRR_WO<=1] <- pal[1];  SRR$colCat_WO[SRR_WO>1] <- pal[2]
SRR$colCat_WOCOMORB <- rep("lightgray",length(SRR_WOCOMORB)); SRR$colCat_WOCOMORB[SRR_WOCOMORB<=1] <- pal[1];  SRR$colCat_WOCOMORB[SRR_WOCOMORB>1] <- pal[2]
SRR$colCat_wglmm <- rep("lightgray",length(SRR_wglmm)); SRR$colCat_wglmm[SRR_wglmm<=1] <- pal[1];  SRR$colCat_wglmm[SRR_wglmm>1] <- pal[2]
SRR$colCat_SRS <- rep("lightgray",length(SRR_SRS)); SRR$colCat_SRS[SRR_SRS<=1] <- pal[1];  SRR$colCat_SRS[SRR_SRS>1] <- pal[2]
SRR$colCat_off <- rep("lightgray",length(SRR_off)); SRR$colCat_off[SRR_off<=1] <- pal[1];  SRR$colCat_off[SRR_off>1] <- pal[2]
SRR$colCat_wglmm_prop <- rep("lightgray",length(SRR_wglmm_prop)); SRR$colCat_wglmm_prop[SRR_wglmm_prop<=1] <- pal[1];  SRR$colCat_wglmm_prop[SRR_wglmm_prop>1] <- pal[2]
SRR$colCat_SRS_prop <- rep("lightgray",length(SRR_SRS_prop)); SRR$colCat_SRS_prop[SRR_SRS_prop<=1] <- pal[1];  SRR$colCat_SRS_prop[SRR_SRS_prop>1] <- pal[2]
SRR$colCat_off_prop <- rep("lightgray",length(SRR_off_prop)); SRR$colCat_off_prop[SRR_off_prop<=1] <- pal[1];  SRR$colCat_off_prop[SRR_off_prop>1] <- pal[2]


## define colours based on category change with respect to gold
SRR$colChange_WO <- rep("lightgray",length(SRR_WO)); SRR$colChange_WO[SRR_WO>1 & SRR_gold<=1] <- pal[2];  SRR$colChange_WO[SRR_WO<=1 & SRR_gold>1] <- pal[1]
SRR$colChange_WOCOMORB <- rep("lightgray",length(SRR_WOCOMORB)); SRR$colChange_WOCOMORB[SRR_WOCOMORB>1 & SRR_gold<=1] <- pal[2];  SRR$colChange_WOCOMORB[SRR_WOCOMORB<=1 & SRR_gold>1] <- pal[1]
SRR$colChange_wglmm <- rep("lightgray",length(SRR_wglmm)); SRR$colChange_wglmm[SRR_wglmm>1 & SRR_gold<=1] <- pal[2];  SRR$colChange_wglmm[SRR_wglmm<=1 & SRR_gold>1] <- pal[1]
SRR$colChange_SRS <- rep("lightgray",length(SRR_SRS)); SRR$colChange_SRS[SRR_SRS>1 & SRR_gold<=1] <- pal[2];  SRR$colChange_SRS[SRR_SRS<=1 & SRR_gold>1] <- pal[1]
SRR$colChange_off <- rep("lightgray",length(SRR_off)); SRR$colChange_off[SRR_off>1 & SRR_gold<=1] <- pal[2];  SRR$colChange_off[SRR_off<=1 & SRR_gold>1] <- pal[1]
SRR$colChange_wglmm_prop <- rep("lightgray",length(SRR_wglmm_prop)); SRR$colChange_wglmm_prop[SRR_wglmm_prop>1 & SRR_gold<=1] <- pal[2];  SRR$colChange_wglmm_prop[SRR_wglmm_prop<=1 & SRR_gold>1] <- pal[1]
SRR$colChange_SRS_prop <- rep("lightgray",length(SRR_SRS_prop)); SRR$colChange_SRS_prop[SRR_SRS_prop>1 & SRR_gold<=1] <- pal[2];  SRR$colChange_SRS_prop[SRR_SRS_prop<=1 & SRR_gold>1] <- pal[1]
SRR$colChange_off_prop <- rep("lightgray",length(SRR_off_prop)); SRR$colChange_off_prop[SRR_off_prop>1 & SRR_gold<=1] <- pal[2];  SRR$colChange_off_prop[SRR_off_prop<=1 & SRR_gold>1] <- pal[1]


## define colours based on rank change with respect to gold
SRR$colRankChange_WO <- colfunc2(length(SRR_WO))[rank(SRR$rankdif_WO_gold)]
SRR$colRankChange_WOCOMORB <- colfunc2(length(SRR_WOCOMORB))[rank(SRR$rankdif_WOCOMORB_gold)]
SRR$colRankChange_wglmm <- colfunc2(length(SRR_wglmm))[rank(SRR$rankdif_wglmm_gold)]
SRR$colRankChange_SRS <- colfunc2(length(SRR_SRS))[rank(SRR$rankdif_SRS_gold)]
SRR$colRankChange_off <- colfunc2(length(SRR_off))[rank(SRR$rankdif_off_gold)]
SRR$colRankChange_wglmm_prop <- colfunc2(length(SRR_wglmm_prop))[rank(SRR$rankdif_wglmm_prop_gold)]
SRR$colRankChange_SRS_prop <- colfunc2(length(SRR_SRS_prop))[rank(SRR$rankdif_SRS_prop_gold)]
SRR$colRankChange_off_prop <- colfunc2(length(SRR_off_prop))[rank(SRR$rankdif_off_prop_gold)]  


## define colours based on change with respect to gold
SRR$colAbsChange_WO <- colfunc(length(SRR_WO))[rank(SRR$SRR_WO-SRR$SRR_gold)]
SRR$colAbsChange_WOCOMORB <- colfunc(length(SRR_WOCOMORB))[rank(SRR$SRR_WOCOMORB-SRR$SRR_gold)]
SRR$colAbsChange_wglmm <- colfunc(length(SRR_wglmm))[rank(SRR$SRR_wglmm-SRR$SRR_gold)]
SRR$colAbsChange_SRS <- colfunc(length(SRR_SRS))[rank(SRR$SRR_SRS-SRR$SRR_gold)]
SRR$colAbsChange_off <- colfunc(length(SRR_off))[rank(SRR$SRR_off-SRR$SRR_gold)]
SRR$colAbsChange_wglmm_prop <- colfunc(length(SRR_wglmm_prop))[rank(SRR$SRR_wglmm_prop-SRR$SRR_gold)]
SRR$colAbsChange_SRS_prop <- colfunc(length(SRR_SRS_prop))[rank(SRR$SRR_SRS_prop-SRR$SRR_gold)]
SRR$colAbsChange_off_prop <- colfunc(length(SRR_off_prop))[rank(SRR$SRR_off_prop-SRR$SRR_gold)]


## define colours based on rank change, with a fraction decision rule (ie top 1/3, top 1/2, top 1/10)
## frac ie: 1/3
frac <- 3
locut <- length(SRR_gold)/frac
hicut <- length(SRR_gold)*(frac-1)/frac
SRR$colRankFract_WO <- rep("lightgray",length(SRR_WO))
SRR$colRankFract_WO[(SRR$rank_WO>locut & SRR$rank_gold<=locut) | (SRR$rank_WO>hicut & SRR$rank_gold<=hicut) ] <- pal[2]
SRR$colRankFract_WO[(SRR$rank_WO<=locut & SRR$rank_gold>locut) | (SRR$rank_WO<=hicut & SRR$rank_gold>hicut) ] <- pal[1]
SRR$colRankFract_WOCOMORB <- rep("lightgray",length(SRR_WOCOMORB))
SRR$colRankFract_WOCOMORB[(SRR$rank_WOCOMORB>locut & SRR$rank_gold<=locut) | (SRR$rank_WOCOMORB>hicut & SRR$rank_gold<=hicut) ] <- pal[2]
SRR$colRankFract_WOCOMORB[(SRR$rank_WOCOMORB<=locut & SRR$rank_gold>locut) | (SRR$rank_WOCOMORB<=hicut & SRR$rank_gold>hicut) ] <- pal[1]
SRR$colRankFract_wglmm <- rep("lightgray",length(SRR_wglmm))
SRR$colRankFract_wglmm[(SRR$rank_wglmm>locut & SRR$rank_gold<=locut) | (SRR$rank_wglmm>hicut & SRR$rank_gold<=hicut) ] <- pal[2]
SRR$colRankFract_wglmm[(SRR$rank_wglmm<=locut & SRR$rank_gold>locut) | (SRR$rank_wglmm<=hicut & SRR$rank_gold>hicut) ] <- pal[1]
SRR$colRankFract_SRS <- rep("lightgray",length(SRR_SRS))
SRR$colRankFract_SRS[(SRR$rank_SRS>locut & SRR$rank_gold<=locut) | (SRR$rank_SRS>hicut & SRR$rank_gold<=hicut) ] <- pal[2]
SRR$colRankFract_SRS[(SRR$rank_SRS<=locut & SRR$rank_gold>locut) | (SRR$rank_SRS<=hicut & SRR$rank_gold>hicut) ] <- pal[1]
SRR$colRankFract_off <- rep("lightgray",length(SRR_off))
SRR$colRankFract_off[(SRR$rank_off>locut & SRR$rank_gold<=locut) | (SRR$rank_off>hicut & SRR$rank_gold<=hicut) ] <- pal[2]
SRR$colRankFract_off[(SRR$rank_off<=locut & SRR$rank_gold>locut) | (SRR$rank_off<=hicut & SRR$rank_gold>hicut) ] <- pal[1]
SRR$colRankFract_wglmm_prop <- rep("lightgray",length(SRR_wglmm_prop))
SRR$colRankFract_wglmm_prop[(SRR$rank_wglmm_prop>locut & SRR$rank_gold<=locut) | (SRR$rank_wglmm_prop>hicut & SRR$rank_gold<=hicut) ] <- pal[2]
SRR$colRankFract_wglmm_prop[(SRR$rank_wglmm_prop<=locut & SRR$rank_gold>locut) | (SRR$rank_wglmm_prop<=hicut & SRR$rank_gold>hicut) ] <- pal[1]
SRR$colRankFract_SRS_prop <- rep("lightgray",length(SRR_SRS_prop))
SRR$colRankFract_SRS_prop[(SRR$rank_SRS_prop>locut & SRR$rank_gold<=locut) | (SRR$rank_SRS_prop>hicut & SRR$rank_gold<=hicut) ] <- pal[2]
SRR$colRankFract_SRS_prop[(SRR$rank_SRS_prop<=locut & SRR$rank_gold>locut) | (SRR$rank_SRS_prop<=hicut & SRR$rank_gold>hicut) ] <- pal[1]
SRR$colRankFract_off_prop <- rep("lightgray",length(SRR_off_prop))
SRR$colRankFract_off_prop[(SRR$rank_off_prop>locut & SRR$rank_gold<=locut) | (SRR$rank_off_prop>hicut & SRR$rank_gold<=hicut) ] <- pal[2]
SRR$colRankFract_off_prop[(SRR$rank_off_prop<=locut & SRR$rank_gold>locut) | (SRR$rank_off_prop<=hicut & SRR$rank_gold>hicut) ] <- pal[1]



## define axis limits
minlim <- 1-max(abs(SRR[,1:5]-1)) # minlim <- min(SRR[,1:3])
maxlim <- 1+max(abs(SRR[,1:5]-1)) # maxlim <- max(SRR[,1:3])
#minlim <- 0.78
#maxlim <- 1.26



####################################################################
#                         SRR-Based Plots                          #
####################################################################

###############################################
## Comparing WO to Gold 
ggplot(data=SRR) + 
  geom_hline(aes(yintercept=1),alpha=0.5) + 
  geom_vline(aes(xintercept=1),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(y=SRR_WO,x=SRR_gold),colour=SRR$colChange_WO,alpha=0.8,size=0.5) + 
  ## next line is for ZOOMBOX
  #geom_rect(aes(ymin=0.92, ymax=1.08, xmin=0.92, xmax=1.08),fill=alpha("grey",0),size=0.3,colour="darkgrey")+
  theme_classic()+
  scale_x_continuous("Complete Data",limits=c(minlim,maxlim), expand = c(0.0, 0.0)) + 
  scale_y_continuous("Full Without SES",limits=c(minlim,maxlim), expand = c(0.0, 0.0)) + 
  ggtitle("(e)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"SRR_gold_WO",suffix,".pdf",sep=""),width=4,height=4)

###############################################
## Comparing WOCOMORB to Gold 
ggplot(data=SRR) + 
  geom_hline(aes(yintercept=1),alpha=0.5) + 
  geom_vline(aes(xintercept=1),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(y=SRR_WOCOMORB,x=SRR_gold),colour=SRR$colChange_WOCOMORB,alpha=0.8,size=0.5) + 
  ## next line is for ZOOMBOX
  #geom_rect(aes(ymin=0.92, ymax=1.08, xmin=0.92, xmax=1.08),fill=alpha("grey",0),size=0.3,colour="darkgrey")+
  theme_classic()+
  scale_x_continuous("Complete Data",limits=c(minlim,maxlim), expand = c(0.0, 0.0)) + 
  scale_y_continuous("Full Without Comorb",limits=c(minlim,maxlim), expand = c(0.0, 0.0)) + 
  ggtitle("(e)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"SRR_gold_WOCOMORB",suffix,".pdf",sep=""),width=4,height=4)

###############################################
## Comparing wglmm to Gold 
ggplot(SRR) + 
  geom_hline(aes(yintercept=1),alpha=0.5) + 
  geom_vline(aes(xintercept=1),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(y=SRR_wglmm,x=SRR_gold),colour=SRR$colChange_wglmm,alpha=0.8,size=0.5) +
  ## next line is for ZOOMBOX
  #geom_rect(aes(ymin=0.92, ymax=1.08, xmin=0.92, xmax=1.08),fill=alpha("grey",0),size=0.3,colour="darkgrey")+
  theme_classic()+
  scale_x_continuous("Complete Data",limits=c(minlim,maxlim), expand = c(0.0, 0.0)) + 
  scale_y_continuous("CSCC",limits=c(minlim,maxlim), expand = c(0.0, 0.0)) +
  ggtitle("(f)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"SRR_gold_wglmm",suffix,".pdf",sep=""),width=4,height=4)


###############################################
## Comparing SRS to Gold 
ggplot(SRR) + 
  geom_hline(aes(yintercept=1),alpha=0.5) + 
  geom_vline(aes(xintercept=1),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(y=SRR_SRS,x=SRR_gold),colour=SRR$colChange_SRS,alpha=0.8,size=0.5) +
  ## next line is for ZOOMBOX
  #geom_rect(aes(ymin=0.92, ymax=1.08, xmin=0.92, xmax=1.08),fill=alpha("grey",0),size=0.3,colour="darkgrey")+
  theme_classic()+
  scale_x_continuous("Complete Data",limits=c(minlim,maxlim), expand = c(0.0, 0.0)) + 
  scale_y_continuous("CSRS",limits=c(minlim,maxlim), expand = c(0.0, 0.0)) +
  ggtitle("(h)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"SRR_gold_SRS",suffix,".pdf",sep=""),width=4,height=4)

###############################################
## Hist of Offset method
ggplot(data=SRR,aes(SRR_off)) + 
  geom_histogram(bins=50)+#binwidth=50) +
  ggtitle("(a)") + 
  scale_x_continuous("Frequency")+ #, expand = c(0.01, 0.0))+
  scale_y_continuous("SRR (CSCC - Offset)")+ #, expand = c(0.01, 0.0))+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"SRR_hist_off",suffix,".pdf",sep=""),width=4,height=4)

## Comparing offset to Gold 
ggplot(SRR) + 
  geom_hline(aes(yintercept=1),alpha=0.5) + 
  geom_vline(aes(xintercept=1),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(y=SRR_off,x=SRR_gold),colour=SRR$colChange_off,alpha=0.5,size=0.5) +
  ## next line is for ZOOMBOX
  #geom_rect(aes(ymin=0.92, ymax=1.08, xmin=0.92, xmax=1.08),fill=alpha("grey",0),size=0.3,colour="darkgrey")+
  theme_classic()+
  scale_x_continuous("Complete Data",limits=c(minlim,maxlim), expand = c(0.0, 0.0)) + 
  scale_y_continuous("CSCC (Offset)",limits=c(minlim,maxlim), expand = c(0.0, 0.0)) +
  ggtitle("(f)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"SRR_gold_off",suffix,".pdf",sep=""),width=4,height=4)

## same but 1x1 and in black
ggplot(SRR) + 
  geom_hline(aes(yintercept=1),alpha=0.5) + 
  geom_vline(aes(xintercept=1),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(y=SRR_off,x=SRR_gold),colour="black",alpha=0.5,size=0.5) +
  theme_classic()+
  scale_x_continuous("Complete Data",limits=c(minlim,maxlim), expand = c(0.0, 0.0)) + 
  scale_y_continuous("CSCC (Offset)",limits=c(minlim,maxlim), expand = c(0.0, 0.0)) +
  ggtitle("") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"SRR_gold_off_single",suffix,".pdf",sep=""),width=4,height=4)

###############################################
## Comparing wglmm_prop to Gold 
ggplot(SRR) + 
  geom_hline(aes(yintercept=1),alpha=0.5) + 
  geom_vline(aes(xintercept=1),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(y=SRR_wglmm_prop,x=SRR_gold),colour=SRR$colChange_wglmm_prop,alpha=0.8,size=0.5) +
  ## next line is for ZOOMBOX
  #geom_rect(aes(ymin=0.92, ymax=1.08, xmin=0.92, xmax=1.08),fill=alpha("grey",0),size=0.3,colour="darkgrey")+
  theme_classic()+
  scale_x_continuous("Complete Data",limits=c(minlim,maxlim), expand = c(0.0, 0.0)) + 
  scale_y_continuous("CSCC",limits=c(minlim,maxlim), expand = c(0.0, 0.0)) +
  ggtitle("(g)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"SRR_gold_wglmm_prop",suffix,".pdf",sep=""),width=4,height=4)

###############################################
## Comparing SRS_prop to Gold 
ggplot(SRR) + 
  geom_hline(aes(yintercept=1),alpha=0.5) + 
  geom_vline(aes(xintercept=1),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(y=SRR_SRS_prop,x=SRR_gold),colour=SRR$colChange_SRS_prop,alpha=0.8,size=0.5) +
  ## next line is for ZOOMBOX
  #geom_rect(aes(ymin=0.92, ymax=1.08, xmin=0.92, xmax=1.08),fill=alpha("grey",0),size=0.3,colour="darkgrey")+
  theme_classic()+
  scale_x_continuous("Complete Data",limits=c(minlim,maxlim), expand = c(0.0, 0.0)) + 
  scale_y_continuous("CSRS",limits=c(minlim,maxlim), expand = c(0.0, 0.0)) +
  ggtitle("(h)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"SRR_gold_SRS_prop",suffix,".pdf",sep=""),width=4,height=4)

###############################################
## Comparing offset prop to Gold 
ggplot(SRR) + 
  geom_hline(aes(yintercept=1),alpha=0.5) + 
  geom_vline(aes(xintercept=1),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(y=SRR_off_prop,x=SRR_gold),colour=SRR$colChange_off_prop,alpha=0.8,size=0.5) +
  ## next line is for ZOOMBOX
  #geom_rect(aes(ymin=0.92, ymax=1.08, xmin=0.92, xmax=1.08),fill=alpha("grey",0),size=0.3,colour="darkgrey")+
  theme_classic()+
  scale_x_continuous("Complete Data",limits=c(minlim,maxlim), expand = c(0.0, 0.0)) + 
  scale_y_continuous("CSCC (Offset)",limits=c(minlim,maxlim), expand = c(0.0, 0.0)) +
  ggtitle("(g)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"SRR_gold_off_prop",suffix,".pdf",sep=""),width=4,height=4)


####################################################################
#                        Hospital Rank Plots                       #
####################################################################

###############################################
## Comparing WO to Gold - RANKS

ggplot(data=SRR) + 
  geom_hline(aes(yintercept=locut),alpha=0.5) + 
  geom_vline(aes(xintercept=locut),alpha=0.5) +
  geom_hline(aes(yintercept=hicut),alpha=0.5) + 
  geom_vline(aes(xintercept=hicut),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(y=rank_WO,x=rank_gold),colour=SRR$colRankFract_WO,alpha=0.8,size=0.5) + 
  theme_classic()+
  scale_x_continuous("Complete Data") + 
  scale_y_continuous("Full Without SES") + 
  ggtitle("Unadjusted\n(a)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18),
        axis.text.y = element_text(margin=margin(l=-5)))
ggsave(paste("",prefix,"SRR_ranks_gold_WO",suffix,".pdf",sep=""),width=4,height=4.3)

###############################################
## Comparing WOCOMORB to Gold - RANKS

ggplot(data=SRR) + 
  geom_hline(aes(yintercept=locut),alpha=0.5) + 
  geom_vline(aes(xintercept=locut),alpha=0.5) +
  geom_hline(aes(yintercept=hicut),alpha=0.5) + 
  geom_vline(aes(xintercept=hicut),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(y=rank_WOCOMORB,x=rank_gold),colour=SRR$colRankFract_WOCOMORB,alpha=0.8,size=0.5) + 
  theme_classic()+
  scale_x_continuous("Complete Data") + 
  scale_y_continuous("Full Without Comorb") + 
  ggtitle("Unadjusted\n(a)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18),
        axis.text.y = element_text(margin=margin(l=-5)))
ggsave(paste("",prefix,"SRR_ranks_gold_WOCOMORB",suffix,".pdf",sep=""),width=4,height=4.3)



###############################################
## Comparing wglmm to Gold - RANKS
ggplot(SRR) + 
  geom_hline(aes(yintercept=locut),alpha=0.5) + 
  geom_vline(aes(xintercept=locut),alpha=0.5) +
  geom_hline(aes(yintercept=hicut),alpha=0.5) + 
  geom_vline(aes(xintercept=hicut),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(y=rank_wglmm,x=rank_gold),colour=SRR$colRankFract_wglmm,alpha=0.8,size=0.5) +
  theme_classic()+
  scale_x_continuous("Complete Data") + 
  scale_y_continuous("CSCC") +
  ggtitle("CSCC (Bal)\n(b)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18),
        axis.text.y = element_text(margin=margin(l=-5)))
ggsave(paste("",prefix,"SRR_ranks_gold_wglmm",suffix,".pdf",sep=""),width=4,height=4.3)

###############################################
## Comparing SRS to Gold - RANKS
ggplot(SRR) + 
  geom_hline(aes(yintercept=locut),alpha=0.5) + 
  geom_vline(aes(xintercept=locut),alpha=0.5) +
  geom_hline(aes(yintercept=hicut),alpha=0.5) + 
  geom_vline(aes(xintercept=hicut),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(y=rank_SRS,x=rank_gold),colour=SRR$colRankFract_SRS,alpha=0.8,size=0.5) +
  theme_classic()+
  scale_x_continuous("Complete Data") + 
  scale_y_continuous("CSRS") +
  ggtitle("CSRS (Bal)\n(d)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18),
        axis.text.y = element_text(margin=margin(l=-5)))
ggsave(paste("",prefix,"SRR_ranks_gold_SRS",suffix,".pdf",sep=""),width=4,height=4.3)

###############################################
## Comparing offset to Gold - RANKS
ggplot(SRR) + 
  geom_hline(aes(yintercept=locut),alpha=0.5) + 
  geom_vline(aes(xintercept=locut),alpha=0.5) +
  geom_hline(aes(yintercept=hicut),alpha=0.5) + 
  geom_vline(aes(xintercept=hicut),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(y=rank_off,x=rank_gold),colour=SRR$colRankFract_off,alpha=0.8,size=0.5) +
  theme_classic()+
  scale_x_continuous("Complete Data") + 
  scale_y_continuous("CSCC (Offset)") +
  ggtitle("CSCC (Bal)\n(b)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18),
        axis.text.y = element_text(margin=margin(l=-5)))
ggsave(paste("",prefix,"SRR_ranks_gold_off",suffix,".pdf",sep=""),width=4,height=4.3)


###############################################
## Comparing wglmm prop to Gold - RANKS
ggplot(SRR) + 
  geom_hline(aes(yintercept=locut),alpha=0.5) + 
  geom_vline(aes(xintercept=locut),alpha=0.5) +
  geom_hline(aes(yintercept=hicut),alpha=0.5) + 
  geom_vline(aes(xintercept=hicut),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(y=rank_wglmm_prop,x=rank_gold),colour=SRR$colRankFract_wglmm_prop,alpha=0.8,size=0.5) +
  theme_classic()+
  scale_x_continuous("Complete Data") + 
  scale_y_continuous("CSCC") +
  ggtitle("CSCC (Prop)\n(c)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18),
        axis.text.y = element_text(margin=margin(l=-5)))
ggsave(paste("",prefix,"SRR_ranks_gold_wglmm_prop",suffix,".pdf",sep=""),width=4,height=4.3)

###############################################
## Comparing SRS prop to Gold - RANKS
ggplot(SRR) + 
  geom_hline(aes(yintercept=locut),alpha=0.5) + 
  geom_vline(aes(xintercept=locut),alpha=0.5) +
  geom_hline(aes(yintercept=hicut),alpha=0.5) + 
  geom_vline(aes(xintercept=hicut),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(y=rank_SRS_prop,x=rank_gold),colour=SRR$colRankFract_SRS_prop,alpha=0.8,size=0.5) +
  theme_classic()+
  scale_x_continuous("Complete Data") + 
  scale_y_continuous("CSRS") +
  ggtitle("CSRS (Prop)\n(d)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18),
        axis.text.y = element_text(margin=margin(l=-5)))
ggsave(paste("",prefix,"SRR_ranks_gold_SRS_prop",suffix,".pdf",sep=""),width=4,height=4.3)

###############################################
## Comparing offset to Gold - RANKS
ggplot(SRR) + 
  geom_hline(aes(yintercept=locut),alpha=0.5) + 
  geom_vline(aes(xintercept=locut),alpha=0.5) +
  geom_hline(aes(yintercept=hicut),alpha=0.5) + 
  geom_vline(aes(xintercept=hicut),alpha=0.5) +
  geom_abline(aes(intercept=0,slope=1),alpha=0.5) +
  geom_point(aes(y=rank_off_prop,x=rank_gold),colour=SRR$colRankFract_off_prop,alpha=0.8,size=0.5) +
  theme_classic()+
  scale_x_continuous("Complete Data") + 
  scale_y_continuous("CSCC (Offset)") +
  ggtitle("CSCC (Prop)\n(c)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18),
        axis.text.y = element_text(margin=margin(l=-5)))
ggsave(paste("",prefix,"SRR_ranks_gold_off_prop",suffix,".pdf",sep=""),width=4,height=4.3)


####################################################################
#                        Ordered SRR Plots                         #
####################################################################
####################################
## ordered plots
SRR.ordered <- SRR[(order(SRR_gold)),]
SRR_order <- seq(1,length(SRR_gold)) 
SRR.ordered <- data.frame(cbind(SRR.ordered,SRR_order))

## WO vs GOLD
ggplot(SRR.ordered) + 
  geom_point( aes(x=SRR_order,y=SRR_WO),alpha=0.8,colour=pal[3],size=2) + 
  geom_point( aes(x=SRR_order,y=SRR_gold),alpha=0.8, colour=pal[4],size=2) +   
  #  ggtitle("Standardized Readmission Rates (WO vs Gold)") +
  theme_classic()+
  scale_x_continuous("Rank") + scale_y_continuous("RASRR",limits=c(minlim,maxlim)) +
  ggtitle("(a)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"ordered_gold_WO",suffix,".pdf",sep=""),width=4,height=4)

## WOCOMORB vs GOLD
ggplot(SRR.ordered) + 
  geom_point( aes(x=SRR_order,y=SRR_WOCOMORB),alpha=0.8,colour=pal[3],size=2) + 
  geom_point( aes(x=SRR_order,y=SRR_gold),alpha=0.8, colour=pal[4],size=2) +   
  #  ggtitle("Standardized Readmission Rates (WOCOMORB vs Gold)") +
  theme_classic()+
  scale_x_continuous("Rank") + scale_y_continuous("RASRR",limits=c(minlim,maxlim)) +
  ggtitle("(a)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"ordered_gold_WOCOMORB",suffix,".pdf",sep=""),width=4,height=4)

## wglmm vs GOLD
ggplot(SRR.ordered) + 
  geom_point( aes(x=SRR_order,y=SRR_wglmm),alpha=0.8,colour=pal[3],size=2) + 
  geom_point( aes(x=SRR_order,y=SRR_gold),alpha=0.8, colour=pal[4],size=2) +   
  #  ggtitle("Standardized Readmission Rates (wglmm vs Gold)") +
  theme_classic()+
  scale_x_continuous("Rank") + scale_y_continuous("RASRR",limits=c(minlim,maxlim)) +
  ggtitle("(b)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"ordered_gold_wglmm",suffix,".pdf",sep=""),width=4,height=4)

## SRS vs GOLD
ggplot(SRR.ordered) + 
  geom_point( aes(x=SRR_order,y=SRR_SRS),alpha=0.8,colour=pal[3],size=2) + 
  geom_point( aes(x=SRR_order,y=SRR_gold),alpha=0.8, colour=pal[4],size=2) +   
  #  ggtitle("Standardized Readmission Rates (SRS vs Gold)") +
  theme_classic()+
  scale_x_continuous("Rank") + scale_y_continuous("RASRR",limits=c(minlim,maxlim)) +
  ggtitle("(c)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"ordered_gold_SRS",suffix,".pdf",sep=""),width=4,height=4)


# offset vs GOLD
ggplot(SRR.ordered) + 
  geom_point( aes(x=SRR_order,y=SRR_off),alpha=0.8,colour=pal[3],size=2) + 
  geom_point( aes(x=SRR_order,y=SRR_gold),alpha=0.8, colour=pal[4],size=2) +   
  #  ggtitle("Standardized Readmission Rates (CSCC [Offset] vs Gold)") +
  theme_classic()+
  scale_x_continuous("Rank") + scale_y_continuous("RASRR",limits=c(minlim,maxlim)) +
  ggtitle("(c)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"ordered_gold_off",suffix,".pdf",sep=""),width=4,height=4)

## wglmm prop vs GOLD
ggplot(SRR.ordered) + 
  geom_point( aes(x=SRR_order,y=SRR_wglmm_prop),alpha=0.8,colour=pal[3],size=2) + 
  geom_point( aes(x=SRR_order,y=SRR_gold),alpha=0.8, colour=pal[4],size=2) +   
  #  ggtitle("Standardized Readmission Rates (wglmm vs Gold)") +
  theme_classic()+
  scale_x_continuous("Rank") + scale_y_continuous("RASRR",limits=c(minlim,maxlim)) +
  ggtitle("(b)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"ordered_gold_wglmm_prop",suffix,".pdf",sep=""),width=4,height=4)

## SRS prop vs GOLD
ggplot(SRR.ordered) + 
  geom_point( aes(x=SRR_order,y=SRR_SRS_prop),alpha=0.8,colour=pal[3],size=2) + 
  geom_point( aes(x=SRR_order,y=SRR_gold),alpha=0.8, colour=pal[4],size=2) +   
  #  ggtitle("Standardized Readmission Rates (SRS vs Gold)") +
  theme_classic()+
  scale_x_continuous("Rank") + scale_y_continuous("RASRR",limits=c(minlim,maxlim)) +
  ggtitle("(c)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"ordered_gold_SRS_prop",suffix,".pdf",sep=""),width=4,height=4)

# offset prop vs GOLD 
ggplot(SRR.ordered) + 
  geom_point( aes(x=SRR_order,y=SRR_off_prop),alpha=0.8,colour=pal[3],size=2) + 
  geom_point( aes(x=SRR_order,y=SRR_gold),alpha=0.8, colour=pal[4],size=2) +   
  #  ggtitle("Standardized Readmission Rates (CSCC [Offset] vs Gold)") +
  theme_classic()+
  scale_x_continuous("Rank") + scale_y_continuous("RASRR",limits=c(minlim,maxlim)) +
  ggtitle("(c)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"ordered_gold_off_prop",suffix,".pdf",sep=""),width=4,height=4)


####################################################################
#                        SES-profile Plots                         #
####################################################################

###############################################
## SES Plot Comparing WO to Gold 
SES_col_WO <- SRR$colChange_WO[SRR$colChange_WO!="lightgray"]

ggplot(data=SRR,aes(x=hosp_race,y=hosp_dual)) + 
  #  geom_hex(bins=25)+
  stat_density_2d(aes(x=hosp_race,y=hosp_dual),colour="black",bins=40,alpha=0.4)+
  #  stat_density_2d(aes(fill = (..level..)), geom="polygon",bins=50,alpha=0.5)+
  geom_point(data=SRR[SRR$colChange_WO!="lightgray",],aes(x=hosp_race,y=hosp_dual),colour=SES_col_WO,size=2.5) + 
  scale_fill_gradientn(colours=c("azure2","black"),name = "Frequency")+
  theme_classic()+
  scale_y_continuous("Percent Dual-Eligible",limits=c(0,1)) + 
  scale_x_continuous("Percent Non-White",limits=c(0,1)) +
  ggtitle("(i)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"SES_gold_WO",suffix,".pdf",sep=""),width=4,height=4)

###############################################
## SES Plot Comparing wglmm to Gold 
SES_col_wglmm <- SRR$colChange_wglmm[SRR$colChange_wglmm!="lightgray"]

ggplot(data=SRR,aes(x=hosp_race,y=hosp_dual)) + 
  #  geom_hex(bins=25)+
  stat_density_2d(aes(x=hosp_race,y=hosp_dual),colour="black",bins=40,alpha=0.4)+
  #  stat_density_2d(aes(fill = (..level..)), geom="polygon",bins=50,alpha=0.5)+
  geom_point(data=SRR[SRR$colChange_wglmm!="lightgray",],aes(x=hosp_race,y=hosp_dual),colour=SES_col_wglmm,size=2.5) + 
  scale_fill_gradientn(colours=c("azure2","black"),name = "Frequency")+
  theme_classic()+
  scale_y_continuous("Percent Dual-Eligible",limits=c(0,1)) + 
  scale_x_continuous("Percent Non-White",limits=c(0,1)) +
  ggtitle("(j)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"SES_gold_wglmm",suffix,".pdf",sep=""),width=4,height=4)


###############################################
## SES Plot Comparing SRS to Gold 
SES_col_SRS <- SRR$colChange_SRS[SRR$colChange_SRS!="lightgray"]

ggplot(data=SRR,aes(x=hosp_race,y=hosp_dual)) + 
  #  geom_hex(bins=25)+
  stat_density_2d(aes(x=hosp_race,y=hosp_dual),colour="black",bins=40,alpha=0.4)+
  #  stat_density_2d(aes(fill = (..level..)), geom="polygon",bins=50,alpha=0.5)+
  geom_point(data=SRR[SRR$colChange_SRS!="lightgray",],aes(x=hosp_race,y=hosp_dual),colour=SES_col_SRS,size=2.5) + 
  scale_fill_gradientn(colours=c("azure2","black"),name = "Frequency")+
  theme_classic()+
  scale_y_continuous("Percent Dual-Eligible",limits=c(0,1)) + 
  scale_x_continuous("Percent Non-White",limits=c(0,1)) +
  ggtitle("(l)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"SES_gold_SRS",suffix,".pdf",sep=""),width=4,height=4)


###############################################
## SES Plot Comparing off to Gold 
SES_col_off <- SRR$colChange_off[SRR$colChange_off!="lightgray"]

ggplot(data=SRR,aes(x=hosp_race,y=hosp_dual)) + 
  #  geom_hex(bins=25)+
  stat_density_2d(aes(x=hosp_race,y=hosp_dual),colour="black",bins=40,alpha=0.4)+
  #  stat_density_2d(aes(fill = (..level..)), geom="polygon",bins=50,alpha=0.5)+
  geom_point(data=SRR[SRR$colChange_off!="lightgray",],aes(x=hosp_race,y=hosp_dual),colour=SES_col_off,size=2.5) + 
  scale_fill_gradientn(colours=c("azure2","black"),name = "Frequency")+
  theme_classic()+
  scale_y_continuous("Percent Dual-Eligible",limits=c(0,1)) + 
  scale_x_continuous("Percent Non-White",limits=c(0,1)) +
  ggtitle("") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"SES_gold_off",suffix,".pdf",sep=""),width=4,height=4)


###############################################
## SES Plot Comparing wglmm_prop prop to Gold 
SES_col_wglmm_prop <- SRR$colChange_wglmm_prop[SRR$colChange_wglmm_prop!="lightgray"]

ggplot(data=SRR,aes(x=hosp_race,y=hosp_dual)) + 
  #  geom_hex(bins=25)+
  stat_density_2d(aes(x=hosp_race,y=hosp_dual),colour="black",bins=40,alpha=0.4)+
  #  stat_density_2d(aes(fill = (..level..)), geom="polygon",bins=50,alpha=0.5)+
  geom_point(data=SRR[SRR$colChange_wglmm_prop!="lightgray",],aes(x=hosp_race,y=hosp_dual),colour=SES_col_wglmm_prop,size=2.5) + 
  scale_fill_gradientn(colours=c("azure2","black"),name = "Frequency")+
  theme_classic()+
  scale_y_continuous("Percent Dual-Eligible",limits=c(0,1)) + 
  scale_x_continuous("Percent Non-White",limits=c(0,1)) +
  ggtitle("(k)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"SES_gold_wglmm_prop",suffix,".pdf",sep=""),width=4,height=4)


###############################################
## SES Plot Comparing SRS_prop prop to Gold 
SES_col_SRS_prop <- SRR$colChange_SRS_prop[SRR$colChange_SRS_prop!="lightgray"]

ggplot(data=SRR,aes(x=hosp_race,y=hosp_dual)) + 
  #  geom_hex(bins=25)+
  stat_density_2d(aes(x=hosp_race,y=hosp_dual),colour="black",bins=40,alpha=0.4)+
  #  stat_density_2d(aes(fill = (..level..)), geom="polygon",bins=50,alpha=0.5)+
  geom_point(data=SRR[SRR$colChange_SRS_prop!="lightgray",],aes(x=hosp_race,y=hosp_dual),colour=SES_col_SRS_prop,size=2.5) + 
  scale_fill_gradientn(colours=c("azure2","black"),name = "Frequency")+
  theme_classic()+
  scale_y_continuous("Percent Dual-Eligible",limits=c(0,1)) + 
  scale_x_continuous("Percent Non-White",limits=c(0,1)) +
  ggtitle("(l)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"SES_gold_SRS_prop",suffix,".pdf",sep=""),width=4,height=4)

###############################################
## SES Plot Comparing off prop to Gold 
SES_col_off_prop <- SRR$colChange_off_prop[SRR$colChange_off_prop!="lightgray"]

ggplot(data=SRR,aes(x=hosp_race,y=hosp_dual)) + 
  #  geom_hex(bins=25)+
  stat_density_2d(aes(x=hosp_race,y=hosp_dual),colour="black",bins=40,alpha=0.4)+
  #  stat_density_2d(aes(fill = (..level..)), geom="polygon",bins=50,alpha=0.5)+
  geom_point(data=SRR[SRR$colChange_off_prop!="lightgray",],aes(x=hosp_race,y=hosp_dual),colour=SES_col_off_prop,size=2.5) + 
  scale_fill_gradientn(colours=c("azure2","black"),name = "Frequency")+
  theme_classic()+
  scale_y_continuous("Percent Dual-Eligible",limits=c(0,1)) + 
  scale_x_continuous("Percent Non-White",limits=c(0,1)) +
  ggtitle("") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"SES_gold_off_prop",suffix,".pdf",sep=""),width=4,height=4)


####################################################################
#                        Comorb-profile Plots                      #
####################################################################

###############################################
## PROF Plot Comparing WOCOMORB to Gold 
PROF_col_WOCOMORB <- SRR$colChange_WOCOMORB[SRR$colChange_WOCOMORB!="lightgray"]

ggplot(data=SRR,aes(x=hosp_renal,y=hosp_protein)) + 
  #  geom_hex(bins=25)+
  stat_density_2d(aes(x=hosp_renal,y=hosp_protein),colour="black",bins=40,alpha=0.5,size=0.2)+
  #  stat_density_2d(aes(fill = (..level..)), geom="polygon",bins=50,alpha=0.5)+
  geom_point(data=SRR[SRR$colChange_WOCOMORB!="lightgray",],aes(x=hosp_renal,y=hosp_protein,colour=as.factor(PROF_col_WOCOMORB)),size=2) + 
  #scale_fill_gradientn(colours=c("azure2","black"),name = "Frequency")+
  scale_color_manual(values=c(pal[1],pal[2]),name="",labels=c("Benefit","Suffer"))+
  theme_classic()+
  scale_y_continuous("Percent Malnutrition",limits=c(0,0.6),expand=c(0,0)) + 
  scale_x_continuous("Percent Renal Failure",limits=c(0,0.6),expand=c(0,0)) +
  ggtitle("(i)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18),legend.position=c(0.15,0.95))
ggsave(paste("",prefix,"PROF_gold_WOCOMORB",suffix,".pdf",sep=""),width=4,height=4)


###############################################
## PROF Plot Comparing wglmm to Gold 
PROF_col_wglmm <- SRR$colChange_wglmm[SRR$colChange_wglmm!="lightgray"]

ggplot(data=SRR,aes(x=hosp_renal,y=hosp_protein)) + 
  #  geom_hex(bins=25)+
  stat_density_2d(aes(x=hosp_renal,y=hosp_protein),colour="black",bins=40,alpha=0.5,size=0.2)+
  #  stat_density_2d(aes(fill = (..level..)), geom="polygon",bins=50,alpha=0.5)+
  geom_point(data=SRR[SRR$colChange_wglmm!="lightgray",],aes(x=hosp_renal,y=hosp_protein,colour=PROF_col_wglmm),size=2) + 
  scale_color_manual(values=c(pal[1],pal[2]),name="",labels=c("Benefit","Suffer"))+
  theme_classic()+
  scale_y_continuous("Percent Malnutrition",limits=c(0,0.6),expand=c(0,0)) + 
  scale_x_continuous("Percent Renal Failure",limits=c(0,0.6),expand=c(0,0)) +
  ggtitle("(j)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18),legend.position=c(0.15,0.95))
ggsave(paste("",prefix,"PROF_gold_wglmm",suffix,".pdf",sep=""),width=4,height=4)


###############################################
## PROF Plot Comparing SRS to Gold 
PROF_col_SRS <- SRR$colChange_SRS[SRR$colChange_SRS!="lightgray"]

ggplot(data=SRR,aes(x=hosp_renal,y=hosp_protein)) + 
  #  geom_hex(bins=25)+
  stat_density_2d(aes(x=hosp_renal,y=hosp_protein),colour="black",bins=40,alpha=0.5,size=0.2)+
  #  stat_density_2d(aes(fill = (..level..)), geom="polygon",bins=50,alpha=0.5)+
  geom_point(data=SRR[SRR$colChange_SRS!="lightgray",],aes(x=hosp_renal,y=hosp_protein,colour=PROF_col_SRS),size=2) + 
  scale_color_manual(values=c(pal[1],pal[2]),name="",labels=c("Benefit","Suffer"))+
  theme_classic()+
  scale_y_continuous("Percent Malnutrition",limits=c(0,0.6),expand=c(0,0)) + 
  scale_x_continuous("Percent Renal Failure",limits=c(0,0.6),expand=c(0,0)) +
  ggtitle("(l)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18),legend.position=c(0.15,0.95))
ggsave(paste("",prefix,"PROF_gold_SRS",suffix,".pdf",sep=""),width=4,height=4)


###############################################
## PROF Plot Comparing off to Gold 
PROF_col_off <- SRR$colChange_off[SRR$colChange_off!="lightgray"]

ggplot(data=SRR,aes(x=hosp_renal,y=hosp_protein)) + 
  #  geom_hex(bins=25)+
  stat_density_2d(aes(x=hosp_renal,y=hosp_protein),colour="black",bins=40,alpha=0.5,size=0.2)+
  #  stat_density_2d(aes(fill = (..level..)), geom="polygon",bins=50,alpha=0.5)+
  geom_point(data=SRR[SRR$colChange_off!="lightgray",],aes(x=hosp_renal,y=hosp_protein,colour=PROF_col_off),size=2) + 
  scale_color_manual(values=c(pal[1],pal[2]),name="",labels=c("Benefit","Suffer"))+
  theme_classic()+
  scale_y_continuous("Percent Malnutrition",limits=c(0,0.6),expand=c(0,0)) + 
  scale_x_continuous("Percent Renal Failure",limits=c(0,0.6),expand=c(0,0)) +
  ggtitle("") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18),legend.position=c(0.15,0.95))
ggsave(paste("",prefix,"PROF_gold_off",suffix,".pdf",sep=""),width=4,height=4)


###############################################
## PROF Plot Comparing wglmm_prop to Gold 
PROF_col_wglmm_prop <- SRR$colChange_wglmm_prop[SRR$colChange_wglmm_prop!="lightgray"]

ggplot(data=SRR,aes(x=hosp_renal,y=hosp_protein)) + 
  #  geom_hex(bins=25)+
  stat_density_2d(aes(x=hosp_renal,y=hosp_protein),colour="black",bins=40,alpha=0.5,size=0.2)+
  #  stat_density_2d(aes(fill = (..level..)), geom="polygon",bins=50,alpha=0.5)+
  geom_point(data=SRR[SRR$colChange_wglmm_prop!="lightgray",],aes(x=hosp_renal,y=hosp_protein,colour=PROF_col_wglmm_prop),size=2) + 
  scale_color_manual(values=c(pal[1],pal[2]),name="",labels=c("Benefit","Suffer"))+
  theme_classic()+
  scale_y_continuous("Percent Malnutrition",limits=c(0,0.6),expand=c(0,0)) + 
  scale_x_continuous("Percent Renal Failure",limits=c(0,0.6),expand=c(0,0)) +
  ggtitle("(k)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18),legend.position=c(0.15,0.95))
ggsave(paste("",prefix,"PROF_gold_wglmm_prop",suffix,".pdf",sep=""),width=4,height=4)


###############################################
## PROF Plot Comparing SRS_prop to Gold 
PROF_col_SRS_prop <- SRR$colChange_SRS_prop[SRR$colChange_SRS_prop!="lightgray"]

ggplot(data=SRR,aes(x=hosp_renal,y=hosp_protein)) + 
  #  geom_hex(bins=25)+
  stat_density_2d(aes(x=hosp_renal,y=hosp_protein),colour="black",bins=40,alpha=0.5,size=0.2)+
  #  stat_density_2d(aes(fill = (..level..)), geom="polygon",bins=50,alpha=0.5)+
  geom_point(data=SRR[SRR$colChange_SRS_prop!="lightgray",],aes(x=hosp_renal,y=hosp_protein,colour=PROF_col_SRS_prop),size=2) + 
  scale_color_manual(values=c(pal[1],pal[2]),name="",labels=c("Benefit","Suffer"))+
  theme_classic()+
  scale_y_continuous("Percent Malnutrition",limits=c(0,0.6),expand=c(0,0)) + 
  scale_x_continuous("Percent Renal Failure",limits=c(0,0.6),expand=c(0,0)) +
  ggtitle("(l)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18),legend.position=c(0.15,0.95))
ggsave(paste("",prefix,"PROF_gold_SRS_prop",suffix,".pdf",sep=""),width=4,height=4)

## PROF Plot Comparing off prop to Gold 
PROF_col_off_prop <- SRR$colChange_off_prop[SRR$colChange_off_prop!="lightgray"]

ggplot(data=SRR,aes(x=hosp_renal,y=hosp_protein)) + 
  #  geom_hex(bins=25)+
  stat_density_2d(aes(x=hosp_renal,y=hosp_protein),colour="black",bins=40,alpha=0.5,size=0.2)+
  #  stat_density_2d(aes(fill = (..level..)), geom="polygon",bins=50,alpha=0.5)+
  geom_point(data=SRR[SRR$colChange_off_prop!="lightgray",],aes(x=hosp_renal,y=hosp_protein,colour=PROF_col_off_prop),size=2) + 
  scale_color_manual(values=c(pal[1],pal[2]),name="",labels=c("Benefit","Suffer"))+
  theme_classic()+
  scale_y_continuous("Percent Malnutrition",limits=c(0,0.6),expand=c(0,0)) + 
  scale_x_continuous("Percent Renal Failure",limits=c(0,0.6),expand=c(0,0)) +
  ggtitle("") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18),legend.position=c(0.15,0.95))
ggsave(paste("",prefix,"PROF_gold_off_prop",suffix,".pdf",sep=""),width=4,height=4)



####################################################################
#                     Changes in SRR and ranks                     #
####################################################################

###############################################
## Hist of SRR Changes from WO to gold
ggplot(data=SRR, aes(SRR_WO-SRR_gold)) + 
  geom_histogram(binwidth=0.003) +
  labs( x='Change in RASRR',y='Frequency')+
  ggtitle("(a)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"hist_WO",suffix,".pdf",sep=""),width=5,height=5)

###############################################
## Hist of SRR Changes from WOCOMORB to gold
ggplot(data=SRR, aes(SRR_WOCOMORB-SRR_gold)) + 
  geom_histogram(binwidth=0.003) +
  labs( x='Change in RASRR',y='Frequency')+
  ggtitle("(a)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"hist_WOCOMORB",suffix,".pdf",sep=""),width=5,height=5)


###############################################
## Hist of SRR Changes from wglmm to gold
ggplot(data=SRR, aes(SRR_wglmm -SRR_gold)) + 
  geom_histogram(binwidth=0.003) +
  labs( x='Change in RASRR',y='Frequency')+
  ggtitle("(b)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"hist_wglmm",suffix,".pdf",sep=""),width=5,height=5)


###############################################
## Hist of SRR Changes from SRS to gold
ggplot(data=SRR, aes(SRR_SRS -SRR_gold)) + 
  geom_histogram(binwidth=0.003) +
  labs( x='Change in RASRR',y='Frequency')+
  ggtitle("(c)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"hist_SRS",suffix,".pdf",sep=""),width=5,height=5)


####################################################################
##  Overlaid hist of SRR Changes 
## To highlight the asymmetry

## SRR diff
ggplot(data=SRR) + 
  geom_density(aes(x=SRR_WOCOMORB-SRR_gold),fill=pal[2],alpha=0.4) +
  geom_density(aes(x=SRR_wglmm-SRR_gold),fill=pal[1],alpha=0.4) +
  geom_density(aes(x=SRR_SRS-SRR_gold),fill=pal[3],alpha=0.4) +
  labs( x='Change in RASRR',y='Density')+
  ggtitle("(a)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"overlay_hist_change",suffix,".pdf",sep=""),width=4,height=4)

## rank diff
ggplot(data=SRR) + 
  geom_density(aes(x=rankdif_WOCOMORB_gold),fill=pal[2],alpha=0.4) +
  geom_density(aes(x=rankdif_wglmm_gold),fill=pal[1],alpha=0.4) +
  geom_density(aes(x=rankdif_SRS_gold),fill=pal[3],alpha=0.4) +
  labs( x='Change in RASRR Rank',y='Density')+
  ggtitle("(b)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("",prefix,"overlay_hist_rankchange",suffix,".pdf",sep=""),width=4,height=4)  


####################################################################
#                               Checks                             #
####################################################################


#########################
##  Histogram of weights 

## Controls
ggplot(data=weights) + 
  geom_histogram(aes(x=M0),bins=40) +
  labs( x='Weight',y='Frequency')+
  ggtitle("Controls") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("hist_M0",suffix,".pdf",sep=""),width=4,height=4)

## Cases
ggplot(data=weights) + 
  geom_histogram(aes(x=M1),bins=40) +
  labs( x='Weight',y='Frequency')+
  ggtitle("Cases") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("hist_M1",suffix,".pdf",sep=""),width=4,height=4)

#########################
##  Boxplot of weights 

weights_long <- data.frame(case_cont <- factor(rep(c("Control","Case"), each=nrow(weights))), wts = c(weights$M0,weights$M1))
ggplot(data=weights_long, aes(x=factor(case_cont), y=wts))+
  geom_boxplot()+
  labs( x='',y='Weight')+
  #ggtitle("Weights for CSCC Sample") + 
  #theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
  ggsave(paste("weights_boxplot",suffix,".pdf",sep=""),width=4,height=4)

####################################################################
##  Scatter of weights vs Predicted random effects (model check)
ggplot() + 
  geom_point(aes(x=weights$M0,y=pred_wglmm),colour=pal[4],size=2.5,alpha=0.4) + 
  geom_point(aes(x=weights$M1,y=pred_wglmm),colour=pal[3],size=2.5,alpha=0.4) + 
  theme_classic()+
  scale_y_continuous("Random Effect") + 
  scale_x_continuous("Weight") +
  #ggtitle("") + 
  #theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
  ggsave(paste("weights_pred",suffix,".pdf",sep=""),width=4,height=4)

####################################################################
##  Scatter of sizes vs Predicted random effects (model check)
ggplot() + 
  geom_point(aes(x=hosp$hosp_sizes,y=pred_wglmm),size=2.5,alpha=0.4) + 
  theme_classic()+
  scale_y_continuous("Random Effect") + 
  scale_x_continuous("Hospital Size") +
  #ggtitle("") + 
  #theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
  ggsave(paste("sizes_pred",suffix,".pdf",sep=""),width=4,height=4)

###############################################
## Weights Plot Comparing WO to Gold  
## (For model checking to see if wglmm is systematically biased in any way)
ggplot(data=SRR,aes(x=hosp_race,y=hosp_dual)) + 
  #  geom_hex(bins=25)+
  #stat_density_2d(aes(x=hosp_race,y=hosp_dual),colour="black",bins=40,alpha=0.4)+
  #  stat_density_2d(aes(fill = (..level..)), geom="polygon",bins=50,alpha=0.5)+
  geom_point(data=weights[SRR$colChange_WO!="lightgray",],aes(x=M0,y=M1),colour=SES_col_WO,size=2.5) + 
  #scale_fill_gradientn(colours=c("azure2","black"),name = "Frequency")+
  theme_classic()+
  scale_y_continuous("Case Weights") + 
  scale_x_continuous("Control Weights") +
  ggtitle("(a)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("weights_gold_WO",suffix,".pdf",sep=""),width=4,height=4)

###############################################
## Weights Plot Comparing WOCOMORB to Gold  
## (For model checking to see if wglmm is systematically biased in any way)
ggplot(data=SRR,aes(x=hosp_race,y=hosp_dual)) + 
  #  geom_hex(bins=25)+
  #stat_density_2d(aes(x=hosp_race,y=hosp_dual),colour="black",bins=40,alpha=0.4)+
  #  stat_density_2d(aes(fill = (..level..)), geom="polygon",bins=50,alpha=0.5)+
  geom_point(data=weights[SRR$colChange_WOCOMORB!="lightgray",],aes(x=M0,y=M1),colour=PROF_col_WOCOMORB,size=2.5) + 
  #scale_fill_gradientn(colours=c("azure2","black"),name = "Frequency")+
  theme_classic()+
  scale_y_continuous("Case Weights") + 
  scale_x_continuous("Control Weights") +
  ggtitle("(a)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("weights_gold_WOCOMORB",suffix,".pdf",sep=""),width=4,height=4)

###############################################
## Weights Plot Comparing wglmm to Gold 
## (For model checking to see if wglmm is systematically biased in any way)
ggplot(data=SRR,aes(x=hosp_race,y=hosp_dual)) + 
  #  geom_hex(bins=25)+
  #stat_density_2d(aes(x=hosp_race,y=hosp_dual),colour="black",bins=40,alpha=0.4)+
  #  stat_density_2d(aes(fill = (..level..)), geom="polygon",bins=50,alpha=0.5)+
  geom_point(data=weights[SRR$colChange_wglmm!="lightgray",],aes(x=M0,y=M1),colour=SES_col_wglmm,size=2.5) + 
  #scale_fill_gradientn(colours=c("azure2","black"),name = "Frequency")+
  theme_classic()+
  scale_y_continuous("Case Weights") + 
  scale_x_continuous("Control Weights") +
  ggtitle("(b)") + 
  theme(plot.title = element_text(hjust=0.5, face="bold",size=18))
ggsave(paste("weights_gold_wglmm",suffix,".pdf",sep=""),width=4,height=4)










####################################################################
#                        Regression Table                          #
####################################################################
est_gold <- fixed_eff$est_gold; est_WO <- fixed_eff$est_WO; est_WOCOMORB <- fixed_eff$est_WOCOMORB; 
SE_gold <- fixed_eff$SE_gold; SE_WO <- fixed_eff$SE_WO; SE_WOCOMORB <- fixed_eff$SE_WOCOMORB;


if(small20==TRUE){
  est_wglmm <- fixed_eff$est_wglmm_20; est_SRS<- fixed_eff$est_SRS_20
  SE_wglmm <- fixed_eff$SE_wglmm_20; SE_SRS<- fixed_eff$SE_SRS_20
  est_off <- fixed_eff$est_off_20; SE_off <- fixed_eff$SE_off_20
  est_wglmm_prop <- fixed_eff$est_wglmm_prop_20; est_SRS_prop<- fixed_eff$est_SRS_prop_20
  SE_wglmm_prop <- fixed_eff$SE_wglmm_prop_20; SE_SRS_prop<- fixed_eff$SE_SRS_prop_20
  est_off_prop <- fixed_eff$est_off_prop_20; SE_off_prop <- fixed_eff$SE_off_prop_20
}else{
  est_wglmm <- fixed_eff$est_wglmm; est_SRS<- fixed_eff$est_SRS
  SE_wglmm <- fixed_eff$SE_wglmm; SE_SRS<- fixed_eff$SE_SRS
  est_off <- fixed_eff$est_off; SE_off <- fixed_eff$SE_off ## offset method
  est_wglmm_prop <- fixed_eff$est_wglmm_prop; est_SRS_prop<- fixed_eff$est_SRS_prop
  SE_wglmm_prop <- fixed_eff$SE_wglmm_prop; SE_SRS_prop<- fixed_eff$SE_SRS_prop
  est_off_prop <- fixed_eff$est_off_prop; SE_off_prop <- fixed_eff$SE_off_prop 
}

fixed_eff <- data.frame(cbind(est_gold,est_WO,est_WOCOMORB,est_wglmm,est_SRS,est_off,est_wglmm_prop,est_SRS_prop,est_off_prop,
                              SE_gold,SE_WO,SE_WOCOMORB,SE_wglmm,SE_SRS,SE_off,SE_wglmm_prop,SE_SRS_prop,SE_off_prop))


## function to get table output:
OR_tex <- function(est,SE,round=TRUE){
  
  ## coef estimates // ci bounds for reg coef
  ests <- est
  Lci <- est-1.96*SE 
  Rci <- est+1.96*SE
  
  ## convert to OR (except sigma)
  ests[-length(ests)] <- exp(ests[-length(ests)])
  Lci[-length(Lci)] <- exp(Lci[-length(Lci)])
  Rci[-length(Rci)] <- exp(Rci[-length(Rci)])
  
  ## round
  ests <- sprintf("%.2f", round(ests,2)  )
  if(round==FALSE){
    CI <- paste(Lci,Rci)
  }else{
    CI <- paste("(",sprintf("%.2f", round(Lci,digits=2)  ),", ",sprintf("%.2f", round(Rci,digits=2)  ),")",sep="")
  }
  
  return(cbind(ests,CI))
}


## get estimates and CIs
OR_gold <- OR_tex( fixed_eff$est_gold,fixed_eff$SE_gold)
OR_WO <- OR_tex( fixed_eff$est_WO,fixed_eff$SE_WO)
OR_WOCOMORB <- OR_tex( fixed_eff$est_WOCOMORB,fixed_eff$SE_WOCOMORB)
OR_wglmm <- OR_tex( fixed_eff$est_wglmm,fixed_eff$SE_wglmm)
OR_SRS <- OR_tex( fixed_eff$est_SRS,fixed_eff$SE_SRS)
OR_off <- OR_tex( fixed_eff$est_off,fixed_eff$SE_off)
OR_wglmm_prop <- OR_tex( fixed_eff$est_wglmm_prop,fixed_eff$SE_wglmm_prop)
OR_SRS_prop <- OR_tex( fixed_eff$est_SRS_prop,fixed_eff$SE_SRS_prop)
OR_off_prop <- OR_tex( fixed_eff$est_off_prop,fixed_eff$SE_off_prop)


### Table of Estimates
table_est <- cbind(OR_gold,
                   OR_WO,OR_WOCOMORB,
                   OR_wglmm,OR_off,OR_SRS,
                   OR_wglmm_prop,OR_off_prop,OR_SRS_prop)

## clean labels
colnames(table_est) <- paste(colnames(table_est),
                             rep(c("Gold",
                                   "WO SES","WO Comorb",
                                   "WGLMM (Bal)","Offset (Bal)", "CSRS (Bal)",
                                   "WGLMM (Prop)","Offset (Prop)","CSRS (Prop)"),each=2))

rownames(table_est) <- c("Intercept","70-74","75-79","80-84","85-89","90-94","95+","Male","AA","Other","Pct Non-White","Dual Eligible","Pct Dual",
                         "4-6","7-13","14+","HomeCare","ICF/SNF","Hospice","Others","Renal Failure","Malnutrition","SD of Random Effect")

## fix NA values
table_est[is.na(table_est) | table_est=="(NA, NA)" | table_est=="NA"] <- "-"



## 10 yr categories: "75-84","85-94","95+" ## this would use ageCat0, 1, 2, 3.



###############################################
## Report results

## all models ## table for supplementary materials
stargazer(table_est,title="Regression Estimates",summary=F,digits=0)

## remove unadjusted-SES analysis (WO)
regtab <- table_est[,-c(3,4)] 

## just balanced-designs
regtab_bal <- regtab[,1:10]
stargazer(regtab_bal,title="Regression Estimates - Balanced Designs",summary=F,digits=0)

## just proportional-designs
regtab_prop <- regtab[,c(1:4,11:16)]
stargazer(regtab_prop,title="Regression Estimates - Proportional Designs",summary=F,digits=0)

## final table 3
#finaltab <- regtab[,c(1:4,5:6,11:16)]
finaltab <- regtab[,c(1:4,5:6,7:8,11:12,9:10)]
stargazer(finaltab,title="Regression Estimates",summary=F,digits=0)



###############################################
## Median Odds Ratio:

## id for the SD of random effects param
sig_id <- length(fixed_eff$est_gold)

## identify SD 
sig_gold <- c(fixed_eff$est_gold[sig_id],fixed_eff$est_gold[sig_id]-1.96*fixed_eff$SE_gold[sig_id],fixed_eff$est_gold[sig_id]+1.96*fixed_eff$SE_gold[sig_id])
sig_WO <- c(fixed_eff$est_WO[sig_id],fixed_eff$est_WO[sig_id]-1.96*fixed_eff$SE_WO[sig_id],fixed_eff$est_WO[sig_id]+1.96*fixed_eff$SE_WO[sig_id])
sig_WOCOMORB <- c(fixed_eff$est_WOCOMORB[sig_id],fixed_eff$est_WOCOMORB[sig_id]-1.96*fixed_eff$SE_WOCOMORB[sig_id],fixed_eff$est_WOCOMORB[sig_id]+1.96*fixed_eff$SE_WOCOMORB[sig_id])
sig_wglmm <- c(fixed_eff$est_wglmm[sig_id],fixed_eff$est_wglmm[sig_id]-1.96*fixed_eff$SE_wglmm[sig_id],fixed_eff$est_wglmm[sig_id]+1.96*fixed_eff$SE_wglmm[sig_id])
sig_SRS <- c(fixed_eff$est_SRS[sig_id],fixed_eff$est_SRS[sig_id]-1.96*fixed_eff$SE_SRS[sig_id],fixed_eff$est_SRS[sig_id]+1.96*fixed_eff$SE_SRS[sig_id])
sig_off <- c(fixed_eff$est_off[sig_id],fixed_eff$est_off[sig_id]-1.96*fixed_eff$SE_off[sig_id],fixed_eff$est_off[sig_id]+1.96*fixed_eff$SE_off[sig_id])
sig_wglmm_prop <- c(fixed_eff$est_wglmm_prop[sig_id],fixed_eff$est_wglmm_prop[sig_id]-1.96*fixed_eff$SE_wglmm_prop[sig_id],fixed_eff$est_wglmm_prop[sig_id]+1.96*fixed_eff$SE_wglmm_prop[sig_id])
sig_SRS_prop <- c(fixed_eff$est_SRS_prop[sig_id],fixed_eff$est_SRS_prop[sig_id]-1.96*fixed_eff$SE_SRS_prop[sig_id],fixed_eff$est_SRS_prop[sig_id]+1.96*fixed_eff$SE_SRS_prop[sig_id])
sig_off_prop <- c(fixed_eff$est_off_prop[sig_id],fixed_eff$est_off_prop[sig_id]-1.96*fixed_eff$SE_off_prop[sig_id],fixed_eff$est_off_prop[sig_id]+1.96*fixed_eff$SE_off_prop[sig_id])

## report SD estimates
round(sig_gold,3)
round(sig_WO,3)
round(sig_WOCOMORB,3)
round(sig_wglmm,3)
round(sig_SRS,3)
round(sig_off[1],3)
round(sig_wglmm_prop[1],3)
round(sig_SRS_prop[1],3)
round(sig_off_prop[1],3)

## report MOR
round(exp(sqrt(2)*sig_gold*qnorm(3/4)),3)
round(exp(sqrt(2)*sig_WO*qnorm(3/4)),3)
round(exp(sqrt(2)*sig_WOCOMORB*qnorm(3/4)),3)
round(exp(sqrt(2)*sig_wglmm*qnorm(3/4)),3)
round(exp(sqrt(2)*sig_SRS*qnorm(3/4)),3)
round(exp(sqrt(2)*sig_wglmm_prop*qnorm(3/4)),3)
round(exp(sqrt(2)*sig_SRS_prop*qnorm(3/4)),3)






####################################################################
#                    Prediction Accuracy Table                     #
####################################################################

halfcut <- length(SRR_gold)/2
loss_WO <- c(round(mean((pred_WO-pred_gold)^2),digits=5),sum(pred_WO>0 & pred_gold<=0),sum(pred_WO<=0 & pred_gold>0),
             round(mean((SRR_WO-SRR_gold)^2),digits=5),sum(SRR_WO>1 & SRR_gold<=1),sum(SRR_WO<=1 & SRR_gold>1),
             round(mean((SRR$rank_WO-SRR$rank_gold)^2),digits=0),sum(SRR$rank_WO>halfcut & SRR$rank_gold<=halfcut),sum(SRR$rank_WO<=halfcut & SRR$rank_gold>halfcut) )

loss_WOCOMORB <- c(round(mean((pred_WOCOMORB-pred_gold)^2),digits=5),sum(pred_WOCOMORB>0 & pred_gold<=0),sum(pred_WOCOMORB<=0 & pred_gold>0),
                   round(mean((SRR_WOCOMORB-SRR_gold)^2),digits=5),sum(SRR_WOCOMORB>1 & SRR_gold<=1),sum(SRR_WOCOMORB<=1 & SRR_gold>1),
                   round(mean((SRR$rank_WOCOMORB-SRR$rank_gold)^2),digits=0),sum(SRR$rank_WOCOMORB>halfcut & SRR$rank_gold<=halfcut),sum(SRR$rank_WOCOMORB<=halfcut & SRR$rank_gold>halfcut) )

loss_wglmm <- c(round(mean((pred_wglmm-pred_gold)^2),digits=5),sum(pred_wglmm>0 & pred_gold<=0),sum(pred_wglmm<=0 & pred_gold>0),
                round(mean((SRR_wglmm-SRR_gold)^2),digits=5),sum(SRR_wglmm>1 & SRR_gold<=1),sum(SRR_wglmm<=1 & SRR_gold>1),
                round(mean((SRR$rank_wglmm-SRR$rank_gold)^2),digits=0),sum(SRR$rank_wglmm>halfcut & SRR$rank_gold<=halfcut),sum(SRR$rank_wglmm<=halfcut & SRR$rank_gold>halfcut) )

loss_SRS <- c(round(mean((pred_SRS-pred_gold)^2),digits=5),sum(pred_SRS>0 & pred_gold<=0),sum(pred_SRS<=0 & pred_gold>0),
              round(mean((SRR_SRS-SRR_gold)^2),digits=5),sum(SRR_SRS>1 & SRR_gold<=1),sum(SRR_SRS<=1 & SRR_gold>1),
              round(mean((SRR$rank_SRS-SRR$rank_gold)^2),digits=0),sum(SRR$rank_SRS>halfcut & SRR$rank_gold<=halfcut),sum(SRR$rank_SRS<=halfcut & SRR$rank_gold>halfcut) )

loss_wglmm_prop <- c(round(mean((pred_wglmm_prop-pred_gold)^2),digits=5),sum(pred_wglmm_prop>0 & pred_gold<=0),sum(pred_wglmm_prop<=0 & pred_gold>0),
                     round(mean((SRR_wglmm_prop-SRR_gold)^2),digits=5),sum(SRR_wglmm_prop>1 & SRR_gold<=1),sum(SRR_wglmm_prop<=1 & SRR_gold>1),
                     round(mean((SRR$rank_wglmm_prop-SRR$rank_gold)^2),digits=0),sum(SRR$rank_wglmm_prop>halfcut & SRR$rank_gold<=halfcut),sum(SRR$rank_wglmm_prop<=halfcut & SRR$rank_gold>halfcut) )

loss_SRS_prop <- c(round(mean((pred_SRS_prop-pred_gold)^2),digits=5),sum(pred_SRS_prop>0 & pred_gold<=0),sum(pred_SRS_prop<=0 & pred_gold>0),
                   round(mean((SRR_SRS_prop-SRR_gold)^2),digits=5),sum(SRR_SRS_prop>1 & SRR_gold<=1),sum(SRR_SRS_prop<=1 & SRR_gold>1),
                   round(mean((SRR$rank_SRS_prop-SRR$rank_gold)^2),digits=0),sum(SRR$rank_SRS_prop>halfcut & SRR$rank_gold<=halfcut),sum(SRR$rank_SRS_prop<=halfcut & SRR$rank_gold>halfcut) )


loss <- c("MSEP","Misclassified > 0","Misclassified <= 0",
          "MSEP","Misclassified > 1.0","Misclassified <= 1.0",
          "MSEP","Misclassified > 50%","Misclassified <= 50%")


table_acc <- cbind(loss,loss_WO,loss_WOCOMORB,loss_wglmm,loss_SRS)
table_acc <- cbind(table_acc,loss_wglmm_prop,loss_SRS_prop)


### Table of Prediction Accuracies
stargazer(table_acc,title="Prediction Accuracy",summary=F,digits=0)


####################################################################
#                     Misclassification Table                     #
####################################################################
## make categories
SRR_cuts <- c(0,1.0,1.05,1.1,10)

cat_gold <- cut(SRR_gold,SRR_cuts,right=TRUE)
cat_WO <- cut(SRR_WO,SRR_cuts,right=TRUE)
cat_WOCOMORB <- cut(SRR_WOCOMORB,SRR_cuts,right=TRUE)
cat_wglmm <- cut(SRR_wglmm,SRR_cuts,right=TRUE)
cat_SRS <- cut(SRR_SRS,SRR_cuts,right=TRUE)
cat_wglmm_prop <- cut(SRR_wglmm_prop,SRR_cuts,right=TRUE)
cat_SRS_prop <- cut(SRR_SRS_prop,SRR_cuts,right=TRUE)

## tables visible in R
# table(cat_gold,cat_WO)
# table(cat_gold,cat_wglmm)

## individual tables
#xtable(as.matrix(table(cat_WO,cat_gold)), caption="Misclassification, Unadjusted")
#xtable(as.matrix(table(cat_wglmm,cat_gold)), caption="Misclassification, CSCC")

## Combined table (gold on top)
misclas_tab <- rbind(as.matrix(table(cat_WO,cat_gold)),
                     as.matrix(table(cat_WOCOMORB,cat_gold)),
                     as.matrix(table(cat_wglmm,cat_gold)),
                     as.matrix(table(cat_SRS,cat_gold)),
                     as.matrix(table(cat_wglmm_prop,cat_gold)),
                     as.matrix(table(cat_SRS_prop,cat_gold)))
xtable(misclas_tab, caption="Misclassification")


###############################################
## Prediction accuracy

## MSEP of bks
mean((pred_WO-pred_gold)^2)
mean((pred_WOCOMORB-pred_gold)^2)
mean((pred_wglmm-pred_gold)^2)
mean((pred_SRS-pred_gold)^2)
mean((pred_wglmm_prop-pred_gold)^2)
mean((pred_SRS_prop-pred_gold)^2)

## MSEP of SRR
mean((SRR_WO-SRR_gold)^2)
mean((SRR_WOCOMORB-SRR_gold)^2)
mean((SRR_wglmm-SRR_gold)^2)
mean((SRR_SRS-SRR_gold)^2)
mean((SRR_wglmm_prop-SRR_gold)^2)
mean((SRR_SRS_prop-SRR_gold)^2)

###############################################
## Correlations

## correlation of bks
cor(pred_WO,pred_gold)
cor(pred_WOCOMORB,pred_gold)
cor(pred_wglmm,pred_gold)
cor(pred_SRS,pred_gold)
cor(pred_wglmm_prop,pred_gold)
cor(pred_SRS_prop,pred_gold)

## correlation of SRRs
cor(SRR_WO,SRR_gold)
cor(SRR_WOCOMORB,SRR_gold)
cor(SRR_wglmm,SRR_gold)
cor(SRR_SRS,SRR_gold)
cor(SRR_wglmm_prop,SRR_gold)
cor(SRR_SRS_prop,SRR_gold)


## spearman rank correlation of bks
cor(pred_WO,pred_gold,method="spearman")
cor(pred_WOCOMORB,pred_gold,method="spearman")
cor(pred_wglmm,pred_gold,method="spearman")
cor(pred_SRS,pred_gold,method="spearman")
cor(pred_wglmm_prop,pred_gold,method="spearman")
cor(pred_SRS_prop,pred_gold,method="spearman")

## spearman rank correlation of SRRs
cor(SRR_WO,SRR_gold,method="spearman")
cor(SRR_WOCOMORB,SRR_gold,method="spearman")
cor(SRR_wglmm,SRR_gold,method="spearman")
cor(SRR_SRS,SRR_gold,method="spearman")
cor(SRR_wglmm_prop,SRR_gold,method="spearman")
cor(SRR_SRS_prop,SRR_gold,method="spearman")


##############################################
#              Comparing Ranks               #
##############################################
IQR(SRR$rankdif_WO_gold)
IQR(SRR$rankdif_WOCOMORB_gold)
IQR(SRR$rankdif_wglmm_gold)
IQR(SRR$rankdif_SRS_gold)
IQR(SRR$rankdif_wglmm_prop_gold)
IQR(SRR$rankdif_SRS_prop_gold)


sd(SRR$rankdif_WO_gold)
sd(SRR$rankdif_WOCOMORB_gold)
sd(SRR$rankdif_wglmm_gold)
sd(SRR$rankdif_SRS_gold)
sd(SRR$rankdif_wglmm_prop_gold)
sd(SRR$rankdif_SRS_prop_gold)



