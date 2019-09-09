####################################################################
#                                                                  #
#                       Case Study Model Fit                       #
#                                                                  #
####################################################################
#                              Dec 2017                            #
####################################################################

library(stargazer)
library(glmmML)
library(gaussquad)
library(ggplot2)
library(wesanderson)
source('draw_samples.R')
source('wglmm_test.R')
source('compute_SRR.R')
set.seed(9999)



####################################################################
#                            Load Data                             #
####################################################################
#source("process_data.R")


####################################################################
#                            Draw Samples                          #
####################################################################
## returns CSCC sample of size n1k=50
sampleCSCC <- data.frame(draw_CSCC_each(clean_full,100)) 
sampleSRS <- data.frame(draw_SRS_each(clean_full,100))
sampleCSCC_prop <- data.frame(draw_CSCCpropto_each(clean_full,0.316055))  ## slightly different fraction to get same sample size
sampleSRS_prop <- data.frame(draw_SRSpropto_each(clean_full,0.316015))  ## slightly different fraction to get same sample size

sampleCSCC_20 <- data.frame(draw_CSCC_each(clean_full,40)) 
sampleSRS_20 <- data.frame(draw_SRS_each(clean_full,40))
sampleCSCC_prop_20 <- data.frame(draw_CSCCpropto_each(clean_full,0.13874))  ## slightly different fraction to get same sample size
sampleSRS_prop_20 <- data.frame(draw_SRSpropto_each(clean_full,0.138815))  ## slightly different fraction to get same sample size


####################################################################
#               Collect Weights for Later Analysis                 #
####################################################################
##   CSCC Weights   ##

## data for population
fullY <- clean_full$event
full_clstr <- clean_full$clinic2
full_tibble <- as_tibble(cbind(fullY,full_clstr))
## summary
N_yk <- full_tibble %>% group_by(full_clstr) %>% summarise(N0k=sum(1-fullY),N1k=sum(fullY))

## balanced CSCC sample
## data for CSCC sample
Y <- sampleCSCC$event
clstr <- sampleCSCC$clinic2
sample_tibble <- as_tibble(cbind(Y,clstr))
## summary
n_yk <- sample_tibble %>% group_by(clstr) %>% summarise(n0k=sum(1-Y),n1k=sum(Y))

## weights
M0 <- c(N_yk$N0k/n_yk$n0k)  
M1 <- c(N_yk$N1k/n_yk$n1k)

## proportional CSCC sample
## data for CSCC sample
Y <- sampleCSCC_prop$event
clstr <- sampleCSCC_prop$clinic2
sample_tibble <- as_tibble(cbind(Y,clstr))
## summary
n_yk <- sample_tibble %>% group_by(clstr) %>% summarise(n0k=sum(1-Y),n1k=sum(Y))

## weights
M0_prop <- c(N_yk$N0k/n_yk$n0k)  
M1_prop <- c(N_yk$N1k/n_yk$n1k)


## sample nyk=20
## data for CSCC sample
Y <- sampleCSCC_20$event
clstr <- sampleCSCC_20$clinic2
sample_tibble <- as_tibble(cbind(Y,clstr))
## summary
n_yk <- sample_tibble %>% group_by(clstr) %>% summarise(n0k=sum(1-Y),n1k=sum(Y))

## weights
M0_20 <- c(N_yk$N0k/n_yk$n0k)  
M1_20 <- c(N_yk$N1k/n_yk$n1k)

## proportional CSCC sample
## data for CSCC sample
Y <- sampleCSCC_prop_20$event
clstr <- sampleCSCC_prop_20$clinic2
sample_tibble <- as_tibble(cbind(Y,clstr))
## summary
n_yk <- sample_tibble %>% group_by(clstr) %>% summarise(n0k=sum(1-Y),n1k=sum(Y))

## weights
M0_prop_20 <- c(N_yk$N0k/n_yk$n0k)  
M1_prop_20 <- c(N_yk$N1k/n_yk$n1k)



####################################################################
#                            Model Fits                            #
####################################################################
####################################################################
### Gold Standard
g_gold <- glmmML(y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition,family=binomial,data=clean_full,cluster=clinic2)
est_gold <- c(g_gold$coefficients,g_gold$sigma)
SE_gold <- c(g_gold$coef.sd,g_gold$sigma.sd)
pred_gold <- g_gold$posterior.modes
SRR_gold <- compute_SRR(g_gold$coefficients,g_gold$sigma,pred_gold,clean_full,'y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition','clinic2')

###################################################################
## w/o SES  (race and dual )
g_WO <- glmmML(y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+prop_nonwhite_std+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition,family=binomial,data=clean_full,cluster=clinic2)
est_WO <- c(g_WO$coefficients,g_WO$sigma)
SE_WO <- c(g_WO$coef.sd,g_WO$sigma.sd)
pred_WO <- g_WO$posterior.modes
SRR_WO <- compute_SRR(g_WO$coefficients,g_WO$sigma,pred_WO,clean_full,'y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+prop_nonwhite_std+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition','clinic2')

###################################################################
## w/o COMORB
g_WOCOMORB <- glmmML(y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4,family=binomial,data=clean_full,cluster=clinic2)
est_WOCOMORB <- c(g_WOCOMORB$coefficients,g_WOCOMORB$sigma)
SE_WOCOMORB <- c(g_WOCOMORB$coef.sd,g_WOCOMORB$sigma.sd)
pred_WOCOMORB <- g_WOCOMORB$posterior.modes
SRR_WOCOMORB <- compute_SRR(g_WOCOMORB$coefficients,g_WOCOMORB$sigma,pred_WOCOMORB,clean_full,'y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4','clinic2')


####################################################################
### SRS
g_SRS <- glmmML(y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition,family=binomial,data=sampleSRS,cluster=clinic2)
est_SRS <- c(g_SRS$coefficients,g_SRS$sigma)
SE_SRS <- c(g_SRS$coef.sd,g_SRS$sigma.sd)
pred_SRS <- g_SRS$posterior.modes
SRR_SRS <- compute_SRR(g_SRS$coefficients,g_SRS$sigma,pred_SRS,sampleSRS,'y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition','clinic2')


####################################################################
### wGLMM
sample <- data.frame(sampleCSCC)
## make observation-level weights
## assign to cases and controls separately
M0_ki <- M0[sample$clinic2]
M1_ki <- M1[sample$clinic2]
M_ki <- c(M1_ki*sample$event+M0_ki*(1-sample$event) )

g_wglmm <- glmmML(y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition,family=binomial,cluster=clinic2,weights=M_ki,data=sample)
est_wglmm <- c(g_wglmm$coefficients,g_wglmm$sigma)
robust <- wglmm_robust(clean_full,sample,'(y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition)',testid='clinic2',est_wglmm)
pred_wglmm <- g_wglmm$posterior.modes
SRR_wglmm <- compute_SRR(g_wglmm$coefficients,g_wglmm$sigma,pred_wglmm,sample,'y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition','clinic2',Mki=M_ki)
cheese <- robust[[1]]
invI <- robust[[2]]
SE_wglmm <- sqrt(diag(invI%*%cheese%*%invI)) 

####################################################################
### offset
g_off <- glmmML(y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition,family=binomial,cluster=clinic2,offset=log(M0_ki/M1_ki),data=sample)
est_off <- c(g_off$coefficients,g_off$sigma)
SE_off <- c(g_off$coef.sd,g_off$sigma.sd)
pred_off <- g_off$posterior.modes
SRR_off <- compute_SRR(g_off$coefficients,g_off$sigma,pred_off,sample,'y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition','clinic2',offset=log(M0_ki/M1_ki))



####################################################################
### SRS_prop -- proportional sampling
g_SRS_prop <- glmmML(y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition,family=binomial,data=sampleSRS_prop,cluster=clinic2)
est_SRS_prop <- c(g_SRS_prop$coefficients,g_SRS_prop$sigma)
SE_SRS_prop <- c(g_SRS_prop$coef.sd,g_SRS_prop$sigma.sd)
pred_SRS_prop <- g_SRS_prop$posterior.modes
SRR_SRS_prop <- compute_SRR(g_SRS_prop$coefficients,g_SRS_prop$sigma,pred_SRS_prop,sampleSRS_prop,'y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition','clinic2')


####################################################################
### wglmm -- proportional sampling
sample <- data.frame(sampleCSCC_prop)
## make observation-level weights
## assign to cases and controls separately
M0_prop_ki <- M0_prop[sample$clinic2]
M1_prop_ki <- M1_prop[sample$clinic2]
M_prop_ki <- c(M1_prop_ki*sample$event+M0_prop_ki*(1-sample$event) )

g_wglmm_prop <- glmmML(y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition,family=binomial,cluster=clinic2,weights=M_prop_ki,data=sample)
est_wglmm_prop <- c(g_wglmm_prop$coefficients,g_wglmm_prop$sigma)
robust <- wglmm_robust(clean_full,sample,'(y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition)',testid='clinic2',est_wglmm_prop)
pred_wglmm_prop <- g_wglmm_prop$posterior.modes
SRR_wglmm_prop <- compute_SRR(g_wglmm_prop$coefficients,g_wglmm_prop$sigma,pred_wglmm_prop,sample,'y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition','clinic2',Mki=M_prop_ki)
cheese <- robust[[1]]
invI <- robust[[2]]
SE_wglmm_prop <- sqrt(diag(invI%*%cheese%*%invI)) 

####################################################################
### offset -- proportional sampling
g_off_prop <- glmmML(y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition,family=binomial,cluster=clinic2,offset=log(M0_prop_ki/M1_prop_ki),data=sample)
est_off_prop <- c(g_off_prop$coefficients,g_off_prop$sigma)
SE_off_prop <- c(g_off_prop$coef.sd,g_off_prop$sigma.sd)
pred_off_prop <- g_off_prop$posterior.modes
SRR_off_prop <- compute_SRR(g_off_prop$coefficients,g_off_prop$sigma,pred_off_prop,sample,'y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition','clinic2',offset=log(M0_prop_ki/M1_prop_ki))





####################################################################
### SRS 20
g_SRS_20 <- glmmML(y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition,family=binomial,data=sampleSRS_20,cluster=clinic2)
est_SRS_20 <- c(g_SRS_20$coefficients,g_SRS_20$sigma)
SE_SRS_20 <- c(g_SRS_20$coef.sd,g_SRS_20$sigma.sd)
pred_SRS_20 <- g_SRS_20$posterior.modes
SRR_SRS_20 <- compute_SRR(g_SRS_20$coefficients,g_SRS_20$sigma,pred_SRS_20,sampleSRS_20,'y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition','clinic2')


####################################################################
### wGLMM 20
sample <- data.frame(sampleCSCC_20)
## make observation-level weights
## assign to cases and controls separately
M0_20_ki <- M0_20[sample$clinic2]
M1_20_ki <- M1_20[sample$clinic2]
M_20_ki <- c(M1_20_ki*sample$event+M0_20_ki*(1-sample$event) )

g_wglmm_20 <- glmmML(y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition,family=binomial,cluster=clinic2,weights=M_20_ki,data=sample)
est_wglmm_20 <- c(g_wglmm_20$coefficients,g_wglmm_20$sigma)
robust <- wglmm_robust(clean_full,sample,'(y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition)',testid='clinic2',est_wglmm_20)
pred_wglmm_20 <- g_wglmm_20$posterior.modes
SRR_wglmm_20 <- compute_SRR(g_wglmm_20$coefficients,g_wglmm_20$sigma,pred_wglmm_20,sample,'y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition','clinic2',Mki=M_20_ki)
cheese <- robust[[1]]
invI <- robust[[2]]
SE_wglmm_20 <- sqrt(diag(invI%*%cheese%*%invI)) 

####################################################################
### offset 20
g_off_20 <- glmmML(y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition,family=binomial,cluster=clinic2,offset=log(M0_20_ki/M1_20_ki),data=sample)
est_off_20 <- c(g_off_20$coefficients,g_off_20$sigma)
SE_off_20 <- c(g_off_20$coef.sd,g_off_20$sigma.sd)
pred_off_20 <- g_off_20$posterior.modes
SRR_off_20 <- compute_SRR(g_off_20$coefficients,g_off_20$sigma,pred_off_20,sample,'y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition','clinic2',offset=log(M0_20_ki/M1_20_ki))



####################################################################
### SRS 20 -- proportional
g_SRS_prop_20 <- glmmML(y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition,family=binomial,data=sampleSRS_prop_20,cluster=clinic2)
est_SRS_prop_20 <- c(g_SRS_prop_20$coefficients,g_SRS_prop_20$sigma)
SE_SRS_prop_20 <- c(g_SRS_prop_20$coef.sd,g_SRS_prop_20$sigma.sd)
pred_SRS_prop_20 <- g_SRS_prop_20$posterior.modes
SRR_SRS_prop_20 <- compute_SRR(g_SRS_prop_20$coefficients,g_SRS_prop_20$sigma,pred_SRS_prop_20,sampleSRS_prop_20,'y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition','clinic2')


####################################################################
### wGLMM 20 -- proportional
sample <- data.frame(sampleCSCC_prop_20)
## make observation-level weights
## assign to cases and controls separately
M0_prop_20_ki <- M0_prop_20[sample$clinic2]
M1_prop_20_ki <- M1_prop_20[sample$clinic2]
M_prop_20_ki <- c(M1_prop_20_ki*sample$event+M0_prop_20_ki*(1-sample$event) )

g_wglmm_prop_20 <- glmmML(y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition,family=binomial,cluster=clinic2,weights=M_prop_20_ki,data=sample)
est_wglmm_prop_20 <- c(g_wglmm_prop_20$coefficients,g_wglmm_prop_20$sigma)
robust <- wglmm_robust(clean_full,sample,'(y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition)',testid='clinic2',est_wglmm_prop_20)
pred_wglmm_prop_20 <- g_wglmm_prop_20$posterior.modes
SRR_wglmm_prop_20 <- compute_SRR(g_wglmm_prop_20$coefficients,g_wglmm_prop_20$sigma,pred_wglmm_prop_20,sample,'y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition','clinic2',Mki=M_prop_20_ki)
cheese <- robust[[1]]
invI <- robust[[2]]
SE_wglmm_prop_20 <- sqrt(diag(invI%*%cheese%*%invI)) 

####################################################################
### offset 20 -- proportional
g_off_prop_20 <- glmmML(y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition,family=binomial,cluster=clinic2,offset=log(M0_prop_20_ki/M1_prop_20_ki),data=sample)
est_off_prop_20 <- c(g_off_prop_20$coefficients,g_off_prop_20$sigma)
SE_off_prop_20 <- c(g_off_prop_20$coef.sd,g_off_prop_20$sigma.sd)
pred_off_prop_20 <- g_off_prop_20$posterior.modes
SRR_off_prop_20 <- compute_SRR(g_off_prop_20$coefficients,g_off_prop_20$sigma,pred_off_prop_20,sample,'y1~ageCat1_5yr+ageCat2_5yr+ageCat3_5yr+ageCat4_5yr+ageCat5_5yr+ageCat6_5yr+male+raceCat1+raceCat2+prop_nonwhite_std+dual+prop_dual_std+LOSCat1+LOSCat2+LOSCat3+disCat1+disCat2+disCat3+disCat4+Renal_failure+Protein_calorie_malnutrition','clinic2',offset=log(M0_prop_20_ki/M1_prop_20_ki))










####################################################################
#                              Output                              #
####################################################################
## output version where WO means without SES 
## also output version where WO means with SES but without comorb (full has both)

## Output Sampling Weights
weights <- cbind(M0,M1,M0_prop,M1_prop,M0_20,M1_20,M0_prop_20,M1_prop_20)
write.table(weights,"weights_results.txt",col.names=TRUE)


## Output Estimates and SEs
## add NAs for length so estimates can be compared
est_WO <- c(est_WO[1:8],rep(NA,2),est_WO[9],NA,est_WO[10:length(est_WO)])
SE_WO <- c(SE_WO[1:8],rep(NA,2),SE_WO[9],NA,SE_WO[10:length(SE_WO)])
est_WOCOMORB <- c(est_WOCOMORB[1:20],rep(NA,2),est_WOCOMORB[21])
SE_WOCOMORB <- c(SE_WOCOMORB[1:20],rep(NA,2),SE_WOCOMORB[21])
fixed_eff <- cbind(est_gold,SE_gold,
              est_WO,SE_WO,
              est_WOCOMORB,SE_WOCOMORB,
              est_wglmm,SE_wglmm,
              est_SRS,SE_SRS,
              est_off,SE_off,
              est_wglmm_prop,SE_wglmm_prop,
              est_SRS_prop,SE_SRS_prop,
              est_off_prop,SE_off_prop,
              est_wglmm_20,SE_wglmm_20,
              est_SRS_20,SE_SRS_20,
              est_off_20,SE_off_20,
              est_wglmm_prop_20,SE_wglmm_prop_20,
              est_SRS_prop_20,SE_SRS_prop_20,
              est_off_prop_20,SE_off_prop_20)
write.table(fixed_eff,"fixed_eff_results.txt",col.names=TRUE)


## Output random effects
pred <- cbind(pred_gold,pred_WO,pred_WOCOMORB,
              pred_wglmm,pred_SRS,pred_off,
              pred_wglmm_prop,pred_SRS_prop,pred_off_prop,
              pred_wglmm_prop_20,pred_SRS_prop_20,pred_off_prop_20,
              pred_wglmm_20,pred_SRS_20,pred_off_20)
write.table(pred,"pred_results.txt",col.names=TRUE)


## Output SRRs
SRR <- cbind(SRR_gold,SRR_WO,SRR_WOCOMORB,
             SRR_wglmm,SRR_SRS,SRR_off,
             SRR_wglmm_prop,SRR_SRS_prop,SRR_off_prop,
             SRR_wglmm_20,SRR_SRS_20,SRR_off_20,
             SRR_wglmm_prop_20,SRR_SRS_prop_20,SRR_off_prop_20)
write.table(SRR,"SRR_results.txt",col.names=TRUE)




