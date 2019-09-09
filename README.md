# hospODS

R code to reproduce methodology and simulations from "Outcome Dependent Sampling in Cluster Correlated Data with Application to Hospital Profiling".

## Simulations

#### 1simA_150.R, 1simA_200.R, 1simA_300.R
Run simulation A with K=150,200 and 300 clusters.

#### 2simB.R
Run simulation B for predictive random effects and RASRRs.

#### 3simC.R
Run simulation C for efficiency comparisons.

#### 12simA_results.R, 22simB_results.R, 32simC_results.R
Produces tables and plots summarizing results from simulations A,B and C, respectively. 

#### samp_sizes.txt
Hospital volumes for generating simulated data.

#### 0masterfile.R
Run all simulations.

## Data Analysis
Code for analyzing CMS data application. Code will not run, because data is not publicly available.

#### case1_datacleaning.R
Clean CMS data for analysis.

#### case2_modelfits.R
Fit readmissions models to CMS data.

#### case3_reportresults.R
Summarize results of analyses via tables and plots.

## Functions
Generic functions to generate simulated data, draw samples, fit models, and compute RASRRs.

#### draw_samples.R
Functions to draw samples of various designs from full populations.

#### model_fits.R
Functions to compute weights and fit models to CSCC samples via PML or Offset method.

#### compute_SRR.R
Function to compute RASRRs (risk adjusted standardized readmissions ratio) based on a given model fit.

#### wglmm_test.R
Function to compute components of sandwich-form standard error for PML method. 
