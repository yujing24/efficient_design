# Designing Efficient Hybrid and Single-Arm Trials: External Control Borrowing and Sample Size Calculation

R codes for simulation studies and data analysis in this paper. 


## File structure

* `simulation`: Folder containing R code for the simulation studies.
  * `1. f_sample_size.R`: Core functions for sample size calculations across different trial design methods, including RCT-only designs (difference-in-means estimator, RCT-only AIPW estimator), hybrid designs, and single-arm designs. This file also implements the sample size calculations for each scenario described in the paper.
  * `2. data_gen_function.R`: Functions for generating simulated data to evaluate performance (power and type I error) of different design methods.
  * `3. est_function.R`: Estimation and assessment functions for all design methods, including the z-test and RCT-only AIPW estimator for randomized trials, ACW estimator for hybrid designs, and EIF-based estimator for single-arm trials.
  * `4. simu_function.R`: Simulation wrapper functions for running each scenario in parallel.
  * `5. f_run_simu.R`: Main script for running the simulations in parallel for each scenario.
  * `6. f_sim_plot.R`: Code for generating plots based on the simulation assessment results.
  * `results`: Folder containing the simulation output files used by `6. f_sim_plot.R` to produce figures.
  * `plots`: Folder containing plots generated from `6. f_sim_plot.R`.
  
* `case_study`: Folder containing code and data for the case study.
  * `case_study.R`: R script for the case study, including synthetic EC dataset generation, sample size calculation, and treatment effect estimation.
  * `chapter15_example.sas7bdat`: Original dataset used in the case study.
  * `sufficient_EC.RData`: Synthetic external control dataset (EC1) used in the case study, with a sample size of 1,000.
  * `insufficient_EC.RData`: Synthetic external control dataset (EC2) used in the case study, with a sample size of 60.

