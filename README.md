# Master Thesis: Analysis of Unsupervised Learning Models for Multiblock Neuroimaging Data: Robustness to Data-Contamination
The code written for this thesis is split into two parts: Python code and R code. The python code is used to read the fMRI brain imaging data and uses the Harvard-Oxford atlas to transform this data to data for specific regions of interest in the brain, and writes the data into a .csv file which can be used with the R code included in this repository. Python is also used to implement the joint ICA method.
The python files are:
  - Read_all_data.py is used for reading the regional homogeneity (REHO), local functional connectivity density (LFCD) and amplitude of low frequency fluctuations (ALFF) fMRI data
  - jointICA_functions.py contains functions used to apply the joint ICA method to data
  - jointICA.py applies joint ICA to the selected datasets

The R code provides implementations of the DISCO-SCA and JIVE models as described in the paper. To use the data, ensure the file paths for the .csv files containing the data is changed to your own directory path. The R code provides the following functionality:
  - Simulation_Study.R creates multiple realisations of simulated data through a DGP where characteristics and dependencies for the original multiblock data are maintained
  - Contamination_Loops.R iterates through the realisations of simulated data and creates noise- and outlier-contaminated datasets for differing levels of contamination
  - DISCO-SCA_model_loops.R loops through all real, simulated and contaminated data and creates a DISCO-SCA model for every dataset
  - DISCO-SCA_result_loops.R reads the created DISCO-SCA models and computes some results and performance measures
  - DISCO-SCA_pearson_tucker_scores.R computes the Tucker's Congruence and Pearson Correlation Coefficients for all DISCO-SCA models
  - JIVE_model_loops.R loops through all real, simulated and contaminated data and creates a DISCO-SCA model for every dataset
  - JIVE_result_loops.R reads the created DISCO-SCA models and computes some results and performance measures
  - JIVE_pearson_tucker_scores.R computes the Tucker's Congruence and Pearson Correlation Coefficients for all DISCO-SCA models
