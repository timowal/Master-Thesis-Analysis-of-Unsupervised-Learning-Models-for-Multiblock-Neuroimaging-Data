### Code For Obtaining DISCO-SCA models results
# Author: Timo Waling (581706tw)


## Install and load packages
install.packages("multiblock")
install.packages("dplyr")
install.packages("tidyr")
install.packages("actuar")
install.packages("ggplot2")

library(multiblock)
library(dplyr)
library(tidyr)
library(actuar)
library(ggplot2)


# Set project directory
project_directory <- "F:/Studie/Thesis/Thesis R Project/"


# Function to compute the individual, common and residual variation for the DISCO-SCA model
explained_variation_components <- function(disco_sca_model) {
  # Extract the needed information from the model
  component_structure <- disco_sca_model$DISCOsca$comdist[[1]]
  explained_var_values <- disco_sca_model$explvar
  
  # Initiate a matrix in which we will store results
  result_matrix <- matrix(0, nrow = 3, ncol = 3)
  rownames(result_matrix) <- c("REHO", "LFCD", "ALFF")
  colnames(result_matrix) <- c("Individual", "Common", "Residual")
  
  # Go through the model's components
  for (i in 1:4) {
    reho_value <- explained_var_values[, i][1] * 0.01
    lfcd_value <- explained_var_values[, i][2] * 0.01
    alff_value <- explained_var_values[, i][3] * 0.01
    
    # Check if this component is an individual or a common component
    if (sum(component_structure[, i]) > 1) { # Common component
      result_matrix[1, 2] <- result_matrix[1, 2] + reho_value
      result_matrix[2, 2] <- result_matrix[2, 2] + lfcd_value
      result_matrix[3, 2] <- result_matrix[3, 2] + alff_value
    } else { # Individual component
      result_matrix[1, 1] <- result_matrix[1, 1] + reho_value
      result_matrix[2, 1] <- result_matrix[2, 1] + lfcd_value
      result_matrix[3, 1] <- result_matrix[3, 1] + alff_value
    }
  }
  
  # Compute the explained variation for the column of residual values
  result_matrix[, 3] <- 1 - (result_matrix[, 1] + result_matrix[, 2])
  
  return(result_matrix)
}


# Function to visualize the results for a single DISCO-SCA model
DISCO_SCA_single_results <- function(disco_sca_model) {
  # Obtain the explained variation values for the DISCO-SCA model, afterward turn it into dataframe
  result_matrix <- explained_variation_components(disco_sca_model)
  result_data <- as.data.frame(as.table(result_matrix))
  
  
  # Ensure the created plot has the correct order of components and data blocks
  colnames(result_data) <- c("Modality", "Component", "Value")
  result_data$Modality <- factor(result_data$Modality, levels = c("REHO", "LFCD", "ALFF"))
  result_data$Component <- factor(result_data$Component, levels = c("Residual", "Individual", "Common"))
  
  # Plot the results
  ggplot(result_data, aes(x = Modality, y = Value, fill = Component)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    scale_fill_manual(values = c("Residual" = "#3f3f3f", "Individual" = "#83c8e4", "Common" = "#2275af")) +
    labs(x = "Data Block", y = "Proportion of Explained Variation", title = "") +
    theme_minimal() +
    theme(axis.text = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 14, face = "bold"),
          legend.title = element_blank(), legend.text = element_text(size = 12))
  
  # Add a column to show the total explained variation as well
  result_matrix <- cbind(result_matrix, Total = 1 - result_matrix[, "Residual"])
  
  # Then also print the results so we can use the numerical values
  print(round(result_matrix, 3))
}

# Function to take DISCO-SCA models for multiple realisations of the simulated data and visualize the results
DISCO_SCA_aggregate_results <- function(n_simulation_datasets, model_directory, file_type) {
  # Initiate array to save the matrices for each realisation
  explvar_all_simulations <- array(0, dim = c(3, 3, n_simulation_datasets), 
                                   dimnames = list(c("REHO", "LFCD", "ALFF"),
                                                   c("Individual", "Common", "Residual"), NULL))
  
  # Create a matrix to save the average values over the realisations
  average_explvar_all <- matrix(0, nrow = 3, ncol = 3)
  rownames(average_explvar_all) <- c("REHO", "LFCD", "ALFF")
  colnames(average_explvar_all) <- c("Individual", "Common", "Residual")
  
  # Loop through the realisations and save the explained variation matrices
  for (i in 1:n_simulation_datasets) {
    filename <- file.path(model_directory, paste0(file_type, i, ".rds"))
    disco_sca_model <- readRDS(filename)
    
    # Extract the explained variation
    result_matrix <- explained_variation_components(disco_sca_model)
    
    # Save the matrix of values for this current realisation
    explvar_all_simulations[, , i] <- result_matrix
    
    average_explvar_all <- average_explvar_all + result_matrix
  }
  
  # Divide to get the average over the simulations, and turn it into a dataframe for visualization
  average_explvar_all <- average_explvar_all / n_simulation_datasets
  average_explvar_data <- as.data.frame(as.table(average_explvar_all))
  
  # Ensure the created plot has the correct order of components and data blocks
  colnames(average_explvar_data) <- c("Modality", "Component", "Value")
  average_explvar_data$Modality <- factor(average_explvar_data$Modality,
                                          levels = c("REHO", "LFCD", "ALFF"))
  
  average_explvar_data$Component <- factor(average_explvar_data$Component,
                                           levels = c("Residual", "Individual", "Common"))
  

  print(ggplot(average_explvar_data, aes(x = Modality, y = Value, fill = Component)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    scale_fill_manual(values = c("Residual" = "#3f3f3f", "Individual" = "#83c8e4", "Common" = "#2275af")) +
    labs(x = "Data Block", y = "Proportion of Explained Variation", title = "") +
    theme_minimal() +
    theme(axis.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 14, face = "bold"),
          legend.title = element_blank(), legend.text = element_text(size = 12)))
  
  # Add a column to show the total explained variation as well
  average_explvar_all <- cbind(average_explvar_all, Total = 1 - average_explvar_all[, "Residual"])
  
  # Then also print the results so we can use the numerical values
  print(round(average_explvar_all, 3))
}

##### DISCO-SCA Results for the Original Data #####
# Load the DISCO-SCA models for the real data
filename_control <- file.path(project_directory, "Results/Real_data/disco_sca_real_control.rds")
filename_asd <- file.path(project_directory, "Results/Real_data/disco_sca_real_asd.rds")
filename_combined <- file.path(project_directory, "Results/Real_data/disco_sca_real_combined.rds")

disco_sca_real_control <- readRDS(filename_control)
disco_sca_real_asd <- readRDS(filename_asd)
disco_sca_real_combined <- readRDS(filename_combined)

# Visualize the results
DISCO_SCA_single_results(disco_sca_real_control)
DISCO_SCA_single_results(disco_sca_real_asd)
DISCO_SCA_single_results(disco_sca_real_combined)


##### Aggregate DISCO-SCA Results for Simulated Data without Contamination #####
model_directory <- paste0(project_directory, "Results/Simulated_Data")
file_type <- "disco_sca_simulated_"
n_simulation_datasets <- 10

DISCO_SCA_aggregate_results(n_simulation_datasets, model_directory, file_type)




##### Aggregate DISCO-SCA Results for Simulated Data with Noise Contamination #####
# Noise Level 01
model_directory <- paste0(project_directory, "Results/Noise_Results/Level_01")
file_type <- "disco_sca_noise_01_simulation_"
n_simulation_datasets <- 10

DISCO_SCA_aggregate_results(n_simulation_datasets, model_directory, file_type)


# Noise Level 05
model_directory <- paste0(project_directory, "Results/Noise_Results/Level_05")
file_type <- "disco_sca_noise_05_simulation_"
n_simulation_datasets <- 10

DISCO_SCA_aggregate_results(n_simulation_datasets, model_directory, file_type)


# Noise Level 10
model_directory <- paste0(project_directory, "Results/Noise_Results/Level_10")
file_type <- "disco_sca_noise_10_simulation_"
n_simulation_datasets <- 10

DISCO_SCA_aggregate_results(n_simulation_datasets, model_directory, file_type)


# Noise Level 25
model_directory <- paste0(project_directory, "Results/Noise_Results/Level_25")
file_type <- "disco_sca_noise_25_simulation_"
n_simulation_datasets <- 10

DISCO_SCA_aggregate_results(n_simulation_datasets, model_directory, file_type)


# Noise Level 50
model_directory <- paste0(project_directory, "Results/Noise_Results/Level_50")
file_type <- "disco_sca_noise_50_simulation_"
n_simulation_datasets <- 10

DISCO_SCA_aggregate_results(n_simulation_datasets, model_directory, file_type)





##### Aggregate JIVE Results for Simulated Data with Outlier Contamination #####
# Outlier Level 01
model_directory <- paste0(project_directory, "Results/Outlier_Results/Level_01")
file_type <- "disco_sca_outlier_01_simulation_"
n_simulation_datasets <- 10

DISCO_SCA_aggregate_results(n_simulation_datasets, model_directory, file_type)


# Outlier Level 05
model_directory <- paste0(project_directory, "Results/Outlier_Results/Level_05")
file_type <- "disco_sca_outlier_05_simulation_"
n_simulation_datasets <- 10

DISCO_SCA_aggregate_results(n_simulation_datasets, model_directory, file_type)


# Outlier Level 10
model_directory <- paste0(project_directory, "Results/Outlier_Results/Level_10")
file_type <- "disco_sca_outlier_10_simulation_"
n_simulation_datasets <- 10

DISCO_SCA_aggregate_results(n_simulation_datasets, model_directory, file_type)


# Outlier Level 25
model_directory <- paste0(project_directory, "Results/Outlier_Results/Level_25")
file_type <- "disco_sca_outlier_25_simulation_"
n_simulation_datasets <- 10

DISCO_SCA_aggregate_results(n_simulation_datasets, model_directory, file_type)


# Outlier Level 50
model_directory <- paste0(project_directory, "Results/Outlier_Results/Level_50")
file_type <- "disco_sca_outlier_50_simulation_"
n_simulation_datasets <- 10

DISCO_SCA_aggregate_results(n_simulation_datasets, model_directory, file_type)
