###########################################
## Master file to run simulations A, B, and C

## includes code for case study 
## will not run without data file


## set TRUE to run case study analysis (not available without data)
case_study=FALSE

### run simulations
## sim A
source('1simA_150.R')
source('1simA_200.R')
source('1simA_300.R')
source('12simA_results.R') ## report results
## sim B
source('2simB.R')
source('22simB_results.R') ## report results
## sim C
source('3simC.R')
source('32simC_results.R') ## report results


## run case study (requires data)
if(case_study==TRUE){
  
  source('case1_datacleaning.R')
  source('case2_modelfits.R')
  source('case3_reportresults.R')
  
}