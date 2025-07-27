### Code For Obtaining JIVE models results
# Author: Timo Waling (581706tw)


## Install and load packages
install.packages("r.jive")
install.packages("dplyr")
install.packages("tidyr")
install.packages("actuar")
install.packages("ggplot2")

library(r.jive)
library(dplyr)
library(tidyr)
library(actuar)
library(ggplot2)


# Set project directory
project_directory <- "F:/Studie/Thesis/Thesis R Project/"

# Function for extracting wanted explained variation values from JIVE models
JIVE_variance_explained <- function(jive_result) {
  # Extract data and components
  data <- jive_result$data
  common <- jive_result$joint
  individual <- jive_result$individual
  
  # Compute variation for each data block and component
  compute_var <- function(x) sum((x - mean(x))^2) / length(x)
  
  # Function to calculate explained variance for each data block
  calc_variance <- function(data, common, individual) {
    residual <- data - common - individual
    total_var <- compute_var(data)
    common_var <- compute_var(common)
    individual_var <- compute_var(individual)
    residual_var <- compute_var(residual)
    
    c(common = common_var / total_var, individual = individual_var / total_var, residual = residual_var / total_var)
  }
  
  # Apply for our three data blocks
  variance_explained <- data.frame(
    Modality = names(data),
    t(sapply(1:3, function(i) calc_variance(data[[i]], common[[i]], individual[[i]])))
  )
  
  return(variance_explained)
}

barplot_JIVE <- function(values) {
  data <- data.frame(
    Modality = rep(c("REHO", "LFCD", "ALFF"), each = 3),
    Value = values,
    Component = rep(c("Common", "Individual", "Residual"), times = 3)
  )
  
  # Ensure that the Feature column is an ordered factor
  data$Modality <- factor(data$Modality, levels = c("REHO", "LFCD", "ALFF"))
  
  # Ensure that the Component column is an ordered factor with "Residual" on top
  data$Component <- factor(data$Component, levels = c("Residual", "Individual", "Common"))
  
  # Create the plot
  ggplot(data, aes(x = Modality, y = Value, fill = Component)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    scale_fill_manual(values = c("Residual" = "#3f3f3f", "Individual" = "#83c8e4", "Common" = "#2275af")) +
    labs(x = "Data Block", y = "Proportion of Explained Variance", title = "") +
    theme_minimal() +
    theme(axis.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 14, face = "bold"),
          legend.title = element_blank(), legend.text = element_text(size = 12))
}

JIVE_single_results <- function(jive_model) {
  # Extract variance components from the JIVE model
  var_result <- JIVE_variance_explained(jive_model)
  
  # Order the values
  ordered_values <- c(
    var_result$joint[var_result$Modality == "reho"],
    var_result$individual[var_result$Modality == "reho"],
    var_result$residual[var_result$Modality == "reho"],
    
    var_result$joint[var_result$Modality == "lfcd"],
    var_result$individual[var_result$Modality == "lfcd"],
    var_result$residual[var_result$Modality == "lfcd"],
    
    var_result$joint[var_result$Modality == "alff"],
    var_result$individual[var_result$Modality == "alff"],
    var_result$residual[var_result$Modality == "alff"]
  )
  
  
  ordered_values_matrix <- matrix(ordered_values, nrow = 3, byrow = TRUE,
                             dimnames = list(c("REHO", "LFCD", "ALFF"), c("Common", "Individual", "Residual")))
    
  ordered_values_matrix <- cbind(ordered_values_matrix, Total = 1 - ordered_values_matrix[, "Residual"])
  
  print(round(ordered_values_matrix, 3))
  
  # Plot the barplot using our existing function
  barplot_JIVE(ordered_values)
}

JIVE_aggregate_results <- function(n_simulation_datasets, model_directory, file_type) {
  # Store the results here
  explained_var_table <- data.frame()
  
  # Loop through models
  for (i in 1:n_simulation_datasets) {
    # Load the JIVE model
    filename <- file.path(model_directory, paste0(file_type, i, ".rds"))
    jive_model <- readRDS(filename)
    
    # Use the function to compute explained variance for all data blocks
    var_result <- JIVE_variance_explained(jive_model)
    
    # Add the model number and Total variance explained
    var_result$Model <- i
    var_result$total <- 1 - var_result$residual
    
    # Reorder columns for consistency
    var_result <- var_result[, c("Model", "Modality", "joint", "individual", "residual", "total")]
    
    # Append to main results table
    explained_var_table <- rbind(explained_var_table, var_result)
  }
  
  # Split up REHO, LFCD and ALFF and obtain the mean values per component type for the data blocsk
  reho_df <- subset(explained_var_table, Modality == "reho")
  alff_df <- subset(explained_var_table, Modality == "alff")
  lfcd_df <- subset(explained_var_table, Modality == "lfcd")
  
  reho_joint_average <- mean(reho_df$joint)
  reho_individual_average <- mean(reho_df$individual)
  reho_residual_average <- mean(reho_df$residual)
  
  lfcd_joint_average <- mean(lfcd_df$joint)
  lfcd_individual_average <- mean(lfcd_df$individual)
  lfcd_residual_average <- mean(lfcd_df$residual)
  
  alff_joint_average <- mean(alff_df$joint)
  alff_individual_average <- mean(alff_df$individual)
  alff_residual_average <- mean(alff_df$residual)
  
  average_variances <- c(reho_joint_average, reho_individual_average, reho_residual_average,
                         lfcd_joint_average, lfcd_individual_average, lfcd_residual_average,
                         alff_joint_average, alff_individual_average, alff_residual_average)
  
  names(average_variances) <- c(
    "REHO_Joint", "REHO_Individual", "REHO_Residual",
    "LFCD_Joint", "LFCD_Individual", "LFCD_Residual",
    "ALFF_Joint", "ALFF_Individual", "ALFF_Residual"
  )
  
  # Create a matrix to print results
  result_matrix <- matrix(average_variances, nrow = 3, byrow = TRUE)
  rownames(result_matrix) <- c("REHO", "LFCD", "ALFF")
  colnames(result_matrix) <- c("Common", "Individual", "Residual")
  
  # Add a fourth column which is the Total
  result_matrix <- cbind(result_matrix, Total = 1 - result_matrix[, "Residual"])
  print(round(result_matrix, 3))
  
  print(barplot_JIVE(average_variances))
}

##### JIVE Results for Real Data #####
# Load the JIVE models for the real data
filename_control <- file.path(project_directory, "Results/Real_data/jive_real_control.rds")
filename_asd <- file.path(project_directory, "Results/Real_data/jive_real_asd.rds")
filename_combined <- file.path(project_directory, "Results/Real_data/jive_real_combined.rds")

jive_real_control <- readRDS(filename_control)
jive_real_asd <- readRDS(filename_asd)
jive_real_combined <- readRDS(filename_combined)

JIVE_single_results(jive_real_control)
JIVE_single_results(jive_real_asd)
JIVE_single_results(jive_real_combined)





##### Aggregate JIVE Results for Simulated Data without Contamination #####
model_directory <- paste0(project_directory, "Results/Simulated_Data")
file_type <- "jive_simulated_"
n_simulation_datasets <- 10

JIVE_aggregate_results(n_simulation_datasets, model_directory, file_type)




##### Aggregate JIVE Results for Simulated Data with Noise Contamination #####
# Noise Level 01
model_directory <- paste0(project_directory, "Results/Noise_Results/Level_01")
file_type <- "jive_noise_01_simulation_"
n_simulation_datasets <- 10

JIVE_aggregate_results(n_simulation_datasets, model_directory, file_type)


# Noise Level 05
model_directory <- paste0(project_directory, "Results/Noise_Results/Level_05")
file_type <- "jive_noise_05_simulation_"
n_simulation_datasets <- 10

JIVE_aggregate_results(n_simulation_datasets, model_directory, file_type)


# Noise Level 10
model_directory <- paste0(project_directory, "Results/Noise_Results/Level_10")
file_type <- "jive_noise_10_simulation_"
n_simulation_datasets <- 10

JIVE_aggregate_results(n_simulation_datasets, model_directory, file_type)


# Noise Level 25
model_directory <- paste0(project_directory, "Results/Noise_Results/Level_25")
file_type <- "jive_noise_25_simulation_"
n_simulation_datasets <- 10

JIVE_aggregate_results(n_simulation_datasets, model_directory, file_type)


# Noise Level 50
model_directory <- paste0(project_directory, "Results/Noise_Results/Level_50")
file_type <- "jive_noise_50_simulation_"
n_simulation_datasets <- 10

JIVE_aggregate_results(n_simulation_datasets, model_directory, file_type)





##### Aggregate JIVE Results for Simulated Data with Outlier Contamination #####
# Outlier Level 01
model_directory <- paste0(project_directory, "Results/Outlier_Results/Level_01")
file_type <- "jive_outlier_01_simulation_"
n_simulation_datasets <- 10

JIVE_aggregate_results(n_simulation_datasets, model_directory, file_type)


# Outlier Level 05
model_directory <- paste0(project_directory, "Results/Outlier_Results/Level_05")
file_type <- "jive_outlier_05_simulation_"
n_simulation_datasets <- 10

JIVE_aggregate_results(n_simulation_datasets, model_directory, file_type)


# Outlier Level 10
model_directory <- paste0(project_directory, "Results/Outlier_Results/Level_10")
file_type <- "jive_outlier_10_simulation_"
n_simulation_datasets <- 10

JIVE_aggregate_results(n_simulation_datasets, model_directory, file_type)


# Outlier Level 25
model_directory <- paste0(project_directory, "Results/Outlier_Results/Level_25")
file_type <- "jive_outlier_25_simulation_"
n_simulation_datasets <- 10

JIVE_aggregate_results(n_simulation_datasets, model_directory, file_type)


# Outlier Level 50
model_directory <- paste0(project_directory, "Results/Outlier_Results/Level_50")
file_type <- "jive_outlier_50_simulation_"
n_simulation_datasets <- 10

JIVE_aggregate_results(n_simulation_datasets, model_directory, file_type)
