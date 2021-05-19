# Bayesian continuous marker adaptive trial design
This repository provides R code to implement the Bayesian adaptive enrichment design with a continuous biomarker.


### Code related to simulation study
The code related to simulation study is located in the folder ./simulations. <br />

The function Sim_Shell_Exponential.R is the top-level function for users to <br />
1) simulate clinical trial data with an exponentially distributed time-to-event endpoint T and one continuous biomarker whose prognostic and predictive effects are determined by the user, by calling simdata.R. <br />
2) implement our proposed Bayesian continuous marker (BCM) design, by calling Bayes_adaptive_design.R which in turn calls several other R functions defined in this folder. <br />
3) implement the comparator frequentist dichotomizing marker (FDM) design, by calling Freq_adaptive_design.R which in turn calls Freq_detect_marker.R. <br />

The R code to run each simulation scenario described in the paper by calling Sim_Shell_Exponential.R can be found in the file ./simulations/one_interim/settingx.R when one interim analysis is planned (x=1,2,...10), and ./simulations/two_interim/settingx.R when two interim analyses are planned (x=1,2,4,9).   <br />
 
The R code to tabulate the simulation results for each scenario, as presented in Web Table 1, can be found in the file ./simulations/one_interim/settingx_table.R and ./simulations/two_interim/settingx_table.R. <br />




### Code related to application to a real trial dataset
The code related to application to a real trial dataset is located in the folder ./application.  <br />

The R code file code_data_application.R describes how to run the BCM design on a given trial dataset by calling BCM_data_application.R defined in the folder. Right censoring of time-to-event data is allowed. The code in this folder accommodates multiple continuous biomarkers whose prognostic and predictive effects are assumed to be independent and additive. In practice, however, usually only one or at most a few biomarkers should be considered.  <br />



### Dependency on other R packages
The following R packages need to be installed before running our proposed BCM design: splines, MASS, survival, Matrix, coda.
