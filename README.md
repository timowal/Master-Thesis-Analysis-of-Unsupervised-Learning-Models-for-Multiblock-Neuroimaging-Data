# Master Thesis: Analysis of Unsupervised Learning Models for Multiblock Neuroimaging Data: Robustness to Data-Contamination
The python code is used to read the fMRI brain imaging data and uses the Harvard-Oxford atlas to transform this data to data for specific regions of interest in the brain, and writes the data into a .csv file which can be used with the R code included in this repository.
The Python files are:
  - readREHO.py is used for reading regional homogeneity fMRI data
  - readLFCD.py is used for reading local functional connectivity density data
  - readALFF.py is used for reading amplitude of low frequency fluctuations data

The R code provides implementations of the DISCO-SCA and JIVE models as described in the paper. To use the data, ensure the file paths for the .csv files containing the data is changed to your own directory path. The R code provides the following functionality:
  - DISCO_SCA.R creates DISCO-SCA models for the data and provides some plots and performance measures
  - Evaluation_DISCO.R provides further performance measures for the DISCO-SCA models
  - JIVE.R creates JIVE models for the data and provides some plots and performance measures
  - Simulation_Study.R creates data through a data generation process where characteristics and dependencies for the original multiblock data are maintained
  - Noise_evaluation.R adds differing levels of noise-contamination to the simulated data
  - Outlier_evaluation.R adds differing levels of outlier-contamination to the simulated data
