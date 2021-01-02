# Continuous marker adaptive trial design
This repository provides the R code to implement our proposed Bayesian adaptive design with continuous biomarkers.


### Code in the main folder 
The function Sim_Shell_Exponential.R is the top-level function for users to <br />
1) simulate clinical trial data with a continuous biomarker and an exponentially distributed time-to-event endpoint T, by calling simdata.R. <br />
2) implement our proposed Bayesian continuous marker (BCM) design, by calling Bayes_adaptive_design.R which in turn calls several other R functions defined in this folder. <br />
3) implement the comparator frequentist dichotomizing marker (FDM) design, by calling Freq_adaptive_design.R which in turn calls Freq_detect_marker.R.



### Code to reproduce simulation studies in the paper_code subfolder
1) The R code to run each scenario described in the paper by calling Sim_Shell_Exponential.R is saved in the file settingx_simulation.R, where x=1,2,...,10. <br />
2) The R code to tabulate the simulation results for each scenario, as presented in Table S1, is saved in the file settingx_table.R, where x=1,2,...,10. <br />
3) The R code to produce Figure 2 and 3 in the paper based on the simulation results is saved respectively in Figure2.R and Figure3.R. <br />
4) The R code to run each scenario described in the supplement by calling Sim_Shell_Exponential.R is saved in the file settingx_K_2_simulation.R, where x=1,2,4,9. <br />
5) The R code to tabulate the simulation results in the supplement for each scenario, as presented in Table S2, is saved in the file settingx_K_2_table.R, where x=1,2,4,9.


