####################################################################
#                                                                  #
#                         Clean/Process Data                       #
#                                                                  #
####################################################################
#                              Feb 2017                            #
####################################################################
## CHF (CHF)
## Readmission (y1)
## Entire Country


library(tidyverse)

## exclude hospitals with fewer than 25 admissions
## according to CMS rules
exclude_less_25 <- TRUE

## all co-morbidities
allcomorbid <- TRUE 

####################################################################
#                        Load Full Data                            #
####################################################################
##
if(allcomorbid==TRUE){
  load("cleanCHF_09_28_2015_allcomorbid.dat")
} else{
  load("cleanCHF_09_28_2015.dat")
}

#clean <- cleanAMI
clean <- cleanCHF

## Set the outcomes
## * y1: readmission within tau days of discharge
## * y2: death within tau days of discharge
##
tau <- 30
clean$y1 <- as.numeric(clean$T1 < tau)
clean$y2 <- as.numeric(clean$T2 < tau)


##
#summary(clean$age)
myCat <- c(75, 85, 95)
clean$ageCat <- 0
for(i in 1:length(myCat)) clean$ageCat[clean$age >= myCat[i]] <- i

##
clean$ageCat0 <- as.numeric(clean$ageCat == 0)
clean$ageCat1 <- as.numeric(clean$ageCat == 1)
clean$ageCat2 <- as.numeric(clean$ageCat == 2)
clean$ageCat3 <- as.numeric(clean$ageCat == 3)

clean$ageCat0_5yr <- as.numeric(clean$age<70)
clean$ageCat1_5yr <- as.numeric(clean$age<75 & clean$age>=70)
clean$ageCat2_5yr <- as.numeric(clean$age<80 & clean$age>=75)
clean$ageCat3_5yr <- as.numeric(clean$age<85 & clean$age>=80)
clean$ageCat4_5yr <- as.numeric(clean$age<90 & clean$age>=85)
clean$ageCat5_5yr <- as.numeric(clean$age<95 & clean$age>=90)
clean$ageCat6_5yr <- as.numeric(clean$age>=95)

##
myCat <- c(4,7,14)
clean$LOSCat <- 0
for(i in 1:length(myCat)) clean$LOSCat[clean$LOS >= myCat[i]] <- i

##
clean$LOSCat0 <- as.numeric(clean$LOSCat == 0)
clean$LOSCat1 <- as.numeric(clean$LOSCat == 1)
clean$LOSCat2 <- as.numeric(clean$LOSCat == 2)
clean$LOSCat3 <- as.numeric(clean$LOSCat == 3)

##
clean$raceCat0 <- as.numeric(clean$race == 0)
clean$raceCat1 <- as.numeric(clean$race == 1)
clean$raceCat2 <- as.numeric(clean$race == 2)

## Discharge: "Home" is the referent
##
clean$disCat0 <- as.numeric(clean$discharge == "Home")
clean$disCat1 <- as.numeric(clean$discharge == "HomeCare")
clean$disCat2 <- as.numeric(clean$discharge == "ICFSNF")
clean$disCat3 <- as.numeric(clean$discharge == "Hospice")
clean$disCat4 <- as.numeric(is.element(clean$discharge, c("Approved swing bed", "Inpatient care", "Long-term care", "Others", "Rehabilitation")))

clean$disCat <- as.factor(0*clean$disCat0+1*clean$disCat1+2*clean$disCat2+3*clean$disCat3+4*clean$disCat4)


## final data (using full data)
clean_full <- clean

####################################################################
#                         Clean Clusters                           #
####################################################################
# which outcome - readmission
clean_full$event <- clean_full$y1
clean_full$clinic <- clean_full$hospital
clean_full$clinic2 <- as.numeric(factor(clean_full$clinic))
clean_full$clinic <- clean_full$clinic2
clean_full <- clean_full[(order(clean_full$clinic2)),]


####################################################################
#             Exclude Tiny Hospital (as per CMS rules)             #
####################################################################
if (exclude_less_25==TRUE){
  clean_full <- as_tibble(clean_full)
  clean_full <- clean_full %>% group_by(clinic2) %>% filter(n() >= 25)
  clean_full <- data.frame(clean_full)
  clean_full$clinic2 <- as.numeric(factor(clean_full$clinic))
  clean_full$clinic <- clean_full$clinic2
}





####################################################################
#                        Cluster Level Data                        #
####################################################################
hosp_level <- as_tibble(clean_full)
hosp <- hosp_level %>% 
  group_by(clinic2) %>%
  summarise(
    hosp_sizes=n(),
    hosp_dual=mean(dual),
    hosp_race=mean(race!=0),
    hosp_rate=mean(event),
    hosp_age=mean(age),
    hosp_gender=mean(male),
    hosp_state=first(state),
    hosp_renal=mean(Renal_failure),
    hosp_protein=mean(Protein_calorie_malnutrition))



## create cluster level covariates in patient level data
clean_full$prop_dual <- hosp$hosp_dual[clean_full$clinic2]
clean_full$prop_nonwhite <- hosp$hosp_race[clean_full$clinic2]
clean_full$hosp_size <- hosp$hosp_sizes[clean_full$clinic2]

## standardize (mean-center and make contrast 10%)
clean_full$prop_dual_std <- 10*(clean_full$prop_dual) 
clean_full$prop_nonwhite_std <- 10*(clean_full$prop_nonwhite) 


rm(list=c("clean","cleanCHF","hosp_level"))

