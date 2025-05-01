### Code For Running Distinctive and Common structures - Simultaneous Components Analysis (DISCO-SCA)
# Author: Timo Waling (581706tw)

# Install relevant packages if not already installed
install.packages("multiblock")
install.packages("neuroim")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("psych")
# Load relevant packages
library(multiblock)
library(neuroim)
library(dplyr)
library(ggplot2)
library(tidyr)
library(psych)

#### Functions used
# Function that fits the input data into a DISCO-SCA model
create_DISCOSCA_model <- function(data_reho, data_lfcd, data_alff, n_comp) {
  # Scale data blocks and put them into one list
  data_list <- list(scale(as.matrix(data_reho)),
                    scale(as.matrix(data_lfcd)),
                    scale(as.matrix(data_alff)))
  # Create the DISCO-SCA model with the given data and given number of components
  result <- disco(data_list, ncomp = n_comp)
  return(result)
}
# Function to plot the explained variance data
plot_explvar_block_component_DISCO <- function(discosca_result) {
  # Get the explained variance data
  explvar <- discosca_result$explvar
  # Sum the variance explained for each block (across all components)
  total_explained_variance <- rowSums(explvar)
  # Calculate the unexplained variance for each block
  unexplained_variance <- 100 - total_explained_variance
  # Create a data frame with the total explained and unexplained variance
  result_df <- data.frame(Block = c("REHO", "LFCD", "ALFF"),
                          TotalVarianceExplained = total_explained_variance,
                          UnexplainedVariance = unexplained_variance)
  # Calculate the averages
  average_explained <- mean(result_df$TotalVarianceExplained)
  average_unexplained <- mean(result_df$UnexplainedVariance)
  # Add an average row for the table (not for the plot)
  average_row <- data.frame(Block = "AVERAGE",
                            TotalVarianceExplained = average_explained,
                            UnexplainedVariance = average_unexplained)
  # Append the average row to the result_df
  result_df_with_avg <- rbind(result_df, average_row)
  print(result_df_with_avg)
  # Tidy format for explained variance per component
  explvar_df <- as.data.frame(explvar)
  colnames(explvar_df) <- colnames(discosca_result$scores)
  explvar_df$Block <- c("REHO", "LFCD", "ALFF")
  explvar_df$Block <- factor(explvar_df$Block,
                             levels = c("REHO", "LFCD", "ALFF"))
  explvar_long <- pivot_longer(explvar_df, cols = -Block,
                               names_to = "Component",
                               values_to = "VarianceExplained")
  # Plot the explained variance by component and block
  plot_expl_var <- ggplot(explvar_long, aes(x = Component,
                                            y = VarianceExplained,
                                            fill = Block)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "COMPONENT", y = "EXPLAINED VARIANCE(%)") +
    scale_fill_manual(values = c("REHO" = "steelblue",
                                 "LFCD" = "lightgreen",
                                 "ALFF" = "coral1")) +
    theme_minimal() + theme(axis.text = element_text(face = "bold"),
                            axis.title = element_text(face = "bold"),
                            legend.title = element_text(face = "bold"),
                            legend.text = element_text(face = "bold"))
  print(plot_expl_var)
}

# Function for computing the Tucker Congruence Coefficient
# between two DISCO-SCA models
compute_tucker_coef <- function(discosca_1, discosca_2) {
  # Get the score matrices from the DISCO-SCA models
  scores_1 <- discosca_1$scores
  scores_2 <- discosca_2$scores
  # Create a matrix to store the Tucker Congruence Coefficients
  tcc_matrix <- matrix(0, ncol(scores_1), ncol(scores_2))
  # Compute TCC for each component pair
  # (each component of model 1 with each component of model 2)
  for (i in 1:ncol(scores_1)) {
    for (j in 1:ncol(scores_2)) {
      # Extract the corresponding components
      component_1 <- scores_1[, i]
      component_2 <- scores_2[, j]
      # Calculate the dot product of the components
      dot_product <- sum(component_1 * component_2)
      # Calculate the norms of the components
      norm_1 <- sqrt(sum(component_1^2))
      norm_2 <- sqrt(sum(component_2^2))
      # Compute the Tucker Congruence Coefficient for the component pair
      tcc_matrix[i, j] <- abs(dot_product) / (norm_1 * norm_2)
    }
  }
  return(tcc_matrix)
}

# Function that computes Pearson correlations
compute_pearson_corr <- function(discosca_result_1, discosca_result_2, n_comp) {
  # Create an empty list to store individual results
  pearson_results <- list()
  # Loop through components
  for (i in 1:n_comp) {
    # Extract the i-th component loadings for both models
    loading_1 <- discosca_result_1$loadings[, i]
    loading_2 <- discosca_result_2$loadings[, i]
    # Compute Pearson correlation for the i-th component
    pearson_correlation <- cor(loading_1, loading_2)
    # Store results
    pearson_results[[paste0("Comp", i)]] <- pearson_correlation
  }
  # Compute the full Pearson correlation matrix for all components
  pearson_matrix <- cor(discosca_result_1$loadings, discosca_result_2$loadings)
  # Return the results as a list
  return(list(individual = pearson_results, matrix = pearson_matrix))
}

#### Load in the original data
# Read the REHO data
reho_data_control <- read.csv("F:/Studie/Thesis/Thesis R Project/data/reho_data/reho_atlas_Control_group.csv", header = FALSE)
reho_data_asd <- read.csv("F:/Studie/Thesis/Thesis R Project/data/reho_data/reho_atlas_ASD_group.csv", header = FALSE)
reho_data_combined <- read.csv("F:/Studie/Thesis/Thesis R Project/data/reho_data/reho_atlas_combined.csv", header = FALSE)
reho_data_combined <- reho_data_combined[-1, ] #Remove first row of 0's
reho_data_control <- reho_data_control[-1, ]
reho_data_asd <- reho_data_asd[-1, ]

# Read the LFCD data
lfcd_data_control <- read.csv("F:/Studie/Thesis/Thesis R Project/data/lfcd_data/lfcd_atlas_Control_group.csv", header = FALSE)
lfcd_data_asd <- read.csv("F:/Studie/Thesis/Thesis R Project/data/lfcd_data/lfcd_atlas_ASD_group.csv", header = FALSE)
lfcd_data_combined<- read.csv("F:/Studie/Thesis/Thesis R Project/data/lfcd_data/lfcd_atlas_combined.csv", header = FALSE)
lfcd_data_control <- lfcd_data_control[-1, ] #Remove first row of 0's
lfcd_data_asd <- lfcd_data_asd[-1, ]
lfcd_data_combined <- lfcd_data_combined[-1, ]

# Read the ALFF data
alff_data_control <- read.csv("F:/Studie/Thesis/Thesis R Project/data/alff_data/alff_atlas_Control_group.csv", header = FALSE)
alff_data_asd <- read.csv("F:/Studie/Thesis/Thesis R Project/data/alff_data/alff_atlas_ASD_group.csv", header = FALSE)
alff_data_combined <- read.csv("F:/Studie/Thesis/Thesis R Project/data/alff_data/alff_atlas_combined.csv", header = FALSE)
alff_data_control <- alff_data_control[-1, ] #Remove first row of 0's
alff_data_asd <- alff_data_asd[-1, ]
alff_data_combined <- alff_data_combined[-1, ]

#### Load in the Simulated  data (see: Simulation_Study.R for the generation of this data)
simulated_data_control <- read.csv("F:/Studie/Thesis/Thesis R Project/data/simulated_data/simulated_data_Control_group.csv", header = FALSE)
simulated_data_asd <- read.csv("F:/Studie/Thesis/Thesis R Project/data/simulated_data/simulated_data_ASD_group.csv", header = FALSE)
simulated_data_combined <- read.csv("F:/Studie/Thesis/Thesis R Project/data/simulated_data/simulated_data_combined.csv", header = FALSE)
simulated_data_control <- simulated_data_control[-1,] #Remove first row of 0's
simulated_data_asd <- simulated_data_asd[-1,]
simulated_data_combined <- simulated_data_combined[-1,]

# Simulation data is still Characters instead of numeric, fix it
simulated_data_control[] <- lapply(simulated_data_control, function(x) as.numeric(as.character(x)))
simulated_data_asd[] <- lapply(simulated_data_asd, function(x) as.numeric(as.character(x)))
simulated_data_combined[] <- lapply(simulated_data_combined, function(x) as.numeric(as.character(x)))

# Extract the REHO, LFCD and ALFF data from the simulated sets of the Control group
simulated_reho_control <- simulated_data_control[, 1:111]
simulated_lfcd_control <- simulated_data_control[, 112:222]
simulated_alff_control <- simulated_data_control[, 223:333]

# Extract the REHO, LFCD and ALFF data from the simulated sets of the ASD group
simulated_reho_asd <- simulated_data_asd[, 1:111]
simulated_lfcd_asd <- simulated_data_asd[, 112:222]
simulated_alff_asd <- simulated_data_asd[, 223:333]

# Extract the REHO, LFCD and ALFF data from the simulated sets of the Combined group
simulated_reho_combined <- simulated_data_combined[, 1:111]
simulated_lfcd_combined <- simulated_data_combined[, 112:222]
simulated_alff_combined <- simulated_data_combined[, 223:333]

#### Create DISCO-SCA models for the different groups
# Set the number of components to be used for DISCO-SCA
n_comp = 3
# Perform DISCO-SCA for Original Control group data:
discosca_result_original_control <- create_DISCOSCA_model(reho_data_control,
                                                          lfcd_data_control,
                                                          alff_data_control,
                                                          n_comp)
# Perform DISCO-SCA for Original ASD group data:
discosca_result_original_asd <- create_DISCOSCA_model(reho_data_asd,
                                                      lfcd_data_asd,
                                                      alff_data_asd,
                                                      n_comp)
# Perform DISCO-SCA for Original Combined data:
discosca_result_original_combined <- create_DISCOSCA_model(reho_data_combined,
                                                           lfcd_data_combined,
                                                           alff_data_combined,
                                                           n_comp)
# Perform DISCO-SCA for Simulated Control group data:
discosca_result_simulated_control <- create_DISCOSCA_model(
                                          simulated_reho_control,
                                          simulated_lfcd_control,
                                          simulated_alff_control,
                                          n_comp)
# Perform DISCO-SCA for Simulated ASD group data:
discosca_result_simulated_asd <- create_DISCOSCA_model(simulated_reho_asd,
                                                       simulated_lfcd_asd,
                                                       simulated_alff_asd,
                                                       n_comp)
# Perform DISCO-SCA for Simulated ASD group data:
discosca_result_simulated_combined <- create_DISCOSCA_model(
                                          simulated_reho_combined,
                                          simulated_lfcd_combined,
                                          simulated_alff_combined,
                                          n_comp)
#### Run the functions to show results
# Show results about explanatory variance
plot_explvar_block_component_DISCO(discosca_result_original_control)
plot_explvar_block_component_DISCO(discosca_result_original_asd)
plot_explvar_block_component_DISCO(discosca_result_simulated_control)
plot_explvar_block_component_DISCO(discosca_result_simulated_asd)
plot_explvar_block_component_DISCO(discosca_result_original_combined)
plot_explvar_block_component_DISCO(discosca_result_simulated_combined)
#### NOTE: For this, first get noise-influenced DISCOSCA models from the Noise_evaluation.R script!
# Noise-Influenced data:
plot_explvar_block_component_DISCO(discosca_result_noise_combined_01)
plot_explvar_block_component_DISCO(discosca_result_noise_combined_05)
plot_explvar_block_component_DISCO(discosca_result_noise_combined_10)
plot_explvar_block_component_DISCO(discosca_result_noise_combined_25)
plot_explvar_block_component_DISCO(discosca_result_noise_combined_50)
#### NOTE: For this, first get outlier-influenced DISCOSCA models from the Outlier_evaluation.R script!
# Outlier-Influenced data:
plot_explvar_block_component_DISCO(discosca_result_outlier_combined_01)
plot_explvar_block_component_DISCO(discosca_result_outlier_combined_05)
plot_explvar_block_component_DISCO(discosca_result_outlier_combined_10)
plot_explvar_block_component_DISCO(discosca_result_outlier_combined_25)
plot_explvar_block_component_DISCO(discosca_result_outlier_combined_50)

#### Evaluation of DISCO-SCA Results
## Compute Tucker Congruence Coefficients
# Comparison Simulated data to noise influenced data
print(compute_tucker_coef(discosca_result_simulated_combined,
                          discosca_result_noise_combined_01))
print(compute_tucker_coef(discosca_result_simulated_combined,
                          discosca_result_noise_combined_05))
print(compute_tucker_coef(discosca_result_simulated_combined,
                          discosca_result_noise_combined_10))
print(compute_tucker_coef(discosca_result_simulated_combined,
                          discosca_result_noise_combined_25))
print(compute_tucker_coef(discosca_result_simulated_combined,
                          discosca_result_noise_combined_50))
# Comparison Simulated data to outlier influenced data
print(compute_tucker_coef(discosca_result_simulated_combined,
                          discosca_result_outlier_combined_01))
print(compute_tucker_coef(discosca_result_simulated_combined,
                          discosca_result_outlier_combined_05))
print(compute_tucker_coef(discosca_result_simulated_combined,
                          discosca_result_outlier_combined_10))
print(compute_tucker_coef(discosca_result_simulated_combined,
                          discosca_result_outlier_combined_25))
print(compute_tucker_coef(discosca_result_simulated_combined,
                          discosca_result_outlier_combined_50))
## Compute Pearson scores
# Comparison Simulated data to noise influenced data
print(compute_pearson_score(discosca_result_simulated_combined,
                            discosca_result_noise_combined_01, 3))
print(compute_pearson_score(discosca_result_simulated_combined,
                            discosca_result_noise_combined_05, 3))
print(compute_pearson_score(discosca_result_simulated_combined,
                            discosca_result_noise_combined_10, 3))
print(compute_pearson_score(discosca_result_simulated_combined,
                            discosca_result_noise_combined_25, 3))
print(compute_pearson_score(discosca_result_simulated_combined,
                            discosca_result_noise_combined_50, 3))
# Comparison Simulated data to outlier influenced data
print(compute_pearson_score(discosca_result_simulated_combined,
                            discosca_result_outlier_combined_01, 3))
print(compute_pearson_score(discosca_result_simulated_combined,
                            discosca_result_outlier_combined_05, 3))
print(compute_pearson_score(discosca_result_simulated_combined,
                            discosca_result_outlier_combined_10, 3))
print(compute_pearson_score(discosca_result_simulated_combined,
                            discosca_result_outlier_combined_25, 3))
print(compute_pearson_score(discosca_result_simulated_combined,
                            discosca_result_outlier_combined_50, 3))
