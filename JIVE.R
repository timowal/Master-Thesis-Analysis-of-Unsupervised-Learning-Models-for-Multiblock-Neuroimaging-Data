### Code for running Joint and Individual Variation Explained (JIVE)
# Author: Timo Waling (581706tw)


# Install relevant packages if not already installed
install.packages("multiblock")
install.packages("neuroim")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("r.jive")
install.packages("reshape2")
# Load relevant packages
library(multiblock)
library(neuroim)
library(dplyr)
library(ggplot2)
library(r.jive)
library(reshape2)

#### Functions used
# Function to calculate variance explained
calculate_variance_explained <- function(matrix, total_variance) {
  svd_matrix <- svd(matrix)
  sum(svd_matrix$d^2) / total_variance * 100
}

# Function to plot the variance explained by Joint and Individual components
plot_explvar_block_component_jive <- function(jive_result) {
  # Extract joint and individual scores
  joint_scores <- jive_result$joint
  individual_scores <- jive_result$individual
  # Convert lists to matrices
  joint_scores_matrices <- lapply(joint_scores, as.matrix)
  individual_scores_matrices <- lapply(individual_scores, as.matrix)
  # Combine matrices for total variance calculation
  joint_scores_combined <- do.call(cbind,
                                   joint_scores_matrices)
  individual_scores_combined <- do.call(cbind,
                                        individual_scores_matrices)
  combined_all_components <- cbind(joint_scores_combined,
                                   individual_scores_combined)
  # Calculate total variance
  total_variance <- sum(svd(combined_all_components)$d^2)
  # Calculate explained variances
  joint_variance_explained <- sapply(joint_scores_matrices,
                                     calculate_variance_explained,
                                     total_variance)
  individual_variance_explained <- sapply(individual_scores_matrices,
                                          calculate_variance_explained,
                                          total_variance)
  # Create data frame for plotting
  explained_df <- data.frame(
    Component = c(paste0("Joint Comp", seq_along(joint_variance_explained)),
                  paste0("Indiv Comp", seq_along(individual_variance_explained))),
    Variance = c(joint_variance_explained, individual_variance_explained),
    Type = rep(c("Joint", "Individual"), c(length(joint_variance_explained),
                                           length(individual_variance_explained))))
  # Plot
  ggplot(explained_df, aes(x = Component, y = Variance, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Component", y = "Total Variance Explained (%)") +
    theme_minimal() +
    theme(axis.text = element_text(face = "bold"),
          axis.title = element_text(face = "bold"))
}

# Function to give values of variance distributed
JIVE_variance_explained <- function(jive_result) {
  # Extract data and component scores
  data <- jive_result$data
  joint <- jive_result$joint
  individual <- jive_result$individual
  # Create empty matrix to store results
  results <- matrix(NA, nrow = length(data), ncol = 3)
  for (i in seq_along(data)) {
    d <- data[[i]]
    j <- joint[[i]]
    ind <- individual[[i]]
    res <- d - j - ind
    # Calculate variances
    total_var <- mean((d - mean(d))^2)
    joint_var <- mean((j - mean(j))^2)
    individual_var <- mean((ind - mean(ind))^2)
    residual_var <- mean((res - mean(res))^2)
    # Store explained variances    
    results[i, ] <- c(joint_var / total_var,
                      individual_var / total_var,
                      residual_var / total_var)
  }
  # Convert to data frame with appropriate column names
  variance_explained <- data.frame(Modality = names(data),
                                   joint = results[, 1],
                                   individual = results[, 2],
                                   residual = results[, 3])
  return(variance_explained)
}

# Function to compute Pearson scores for JIVE models
JIVE_pearson_score <- function(model_1, model_2) {
  # Extract joint components from both models
  joint_model_1 <- model_1$joint
  joint_model_2 <- model_2$joint
  # Ensure that both are matrices
  joint_model_1 <- lapply(joint_model_1, as.matrix)
  joint_model_2 <- lapply(joint_model_2, as.matrix)
  # Create an empty vector to store the Pearson correlations
  result <- numeric(length = length(joint_model_1))
  # Compute Pearson correlation for each component (across rows for each component)
  for (i in 1:length(joint_model_1)) {
    # Flatten the matrices and calculate the correlation (treating the data as vectors)
    result[i] <- cor(as.vector(joint_model_1[[i]]), as.vector(joint_model_2[[i]]))
  }
  return(result)
}

# Function for computing Tucker Congruence Coefficients for JIVE models
JIVE_tucker_congruence <- function(model_1, model_2) {
  # Extract joint components from both models
  joint_model_1 <- model_1$joint
  joint_model_2 <- model_2$joint
  # Ensure that both are matrices
  joint_model_1 <- lapply(joint_model_1, as.matrix)
  joint_model_2 <- lapply(joint_model_2, as.matrix)
  # Create an empty vector to store the Tucker Congruence Coefficients
  tucker_corr <- numeric(length = length(joint_model_1))
  # Compute Tucker Congruence for each component
  for (i in 1:length(joint_model_1)) {
    # Flatten the matrices into vectors
    vec_1 <- as.vector(joint_model_1[[i]])
    vec_2 <- as.vector(joint_model_2[[i]])
    # Compute Tucker Congruence Coefficient
    numerator <- sum(vec_1 * vec_2)^2
    denominator <- sum(vec_1^2) * sum(vec_2^2)
    # Store the result
    tucker_corr[i] <- numerator / denominator
  }
  return(tucker_corr)
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

#### Creating datablocks for JIVE input 
# Transpose and scale data, and put it into combined datablock frame
data_blocks_original_control <- list(reho = t(reho_data_control), lfcd = t(lfcd_data_control), 
                             alff = t(alff_data_control))
data_blocks_original_asd <- list(reho = t(reho_data_asd), lfcd = t(lfcd_data_asd), 
                                alff = t(alff_data_asd))
data_blocks_original_combined <- list(reho = t(reho_data_combined), lfcd = t(lfcd_data_combined), 
                                 alff = t(alff_data_combined))
data_blocks_simulated_control <- list(reho = t(simulated_reho_control), lfcd = t(simulated_lfcd_control), 
                                     alff = t(simulated_alff_control))
data_blocks_simulated_asd <- list(reho = t(simulated_reho_asd), lfcd = t(simulated_lfcd_asd), 
                                      alff = t(simulated_alff_asd))
data_blocks_simulated_combined <- list(reho = t(simulated_reho_combined), lfcd = t(simulated_lfcd_combined), 
                                  alff = t(simulated_alff_combined))

#### Create JIVE models for the different groups
jive_original_control <- jive(data_blocks_original_control,
                            method = "given", rankJ = 3, rankA = c(2, 2, 2))
jive_original_asd <- jive(data_blocks_original_asd,
                            method = "given", rankJ = 3, rankA = c(2, 2, 2))
jive_original_combined <- jive(data_blocks_original_combined,
                            method = "given", rankJ = 3, rankA = c(2, 2, 2))
jive_simulated_control <- jive(data_blocks_simulated_control,
                            method = "given", rankJ = 3, rankA = c(2, 2, 2))
jive_simulated_asd <- jive(data_blocks_simulated_asd,
                            method = "given", rankJ = 3, rankA = c(2, 2, 2))
jive_simulated_combined <- jive(data_blocks_simulated_combined,
                            method = "given", rankJ = 3, rankA = c(2, 2, 2))

# Plot the results from JIVE
plot(jive_original_control)
plot(jive_original_asd)
plot(jive_simulated_control)
plot(jive_simulated_asd)
plot(jive_original_combined)
plot(jive_simulated_combined)

# Plot the fractions of explained variance for JIVE models
JIVE_variance_explained(jive_original_control)
JIVE_variance_explained(jive_original_asd)
JIVE_variance_explained(jive_simulated_control)
JIVE_variance_explained(jive_simulated_asd)
JIVE_variance_explained(jive_original_combined)
JIVE_variance_explained(jive_simulated_combined)

# Example usage with jive_original_asd:
plot_explvar_block_component_jive(jive_original_control)
plot_explvar_block_component_jive(jive_original_asd)
plot_explvar_block_component_jive(jive_simulated_control)
plot_explvar_block_component_jive(jive_simulated_asd)
plot_explvar_block_component_jive(jive_original_combined)
plot_explvar_block_component_jive(jive_simulated_combined)

# Plot the average loadings per data block for JIVE
plot_ROI_loadings(jive_original_control)
plot_ROI_loadings(jive_original_asd)
plot_ROI_loadings(jive_simulated_control)
plot_ROI_loadings(jive_simulated_asd)
plot_ROI_loadings(jive_original_combined)
plot_ROI_loadings(jive_simulated_combined)



#### NOTE: For this, first get noise-influenced JIVE models from the Noise_evaluation.R script!
# Plot the results of the noise-influenced JIVE models
plot(jive_simulated_combined)
plot(jive_noise_combined_01)
plot(jive_noise_combined_05)
plot(jive_noise_combined_10)
plot(jive_noise_combined_25)
plot(jive_noise_combined_50)

# Plot the fractions of explained variance for noise-contaminated JIVE models
JIVE_variance_explained(jive_simulated_combined)
JIVE_variance_explained(jive_noise_combined_01)
JIVE_variance_explained(jive_noise_combined_05)
JIVE_variance_explained(jive_noise_combined_10)
JIVE_variance_explained(jive_noise_combined_25)
JIVE_variance_explained(jive_noise_combined_50)

# Compute Pearson scores for noise-contaminated data
JIVE_pearson_score(jive_simulated_combined, jive_noise_combined_01)
JIVE_pearson_score(jive_simulated_combined, jive_noise_combined_05)
JIVE_pearson_score(jive_simulated_combined, jive_noise_combined_10)
JIVE_pearson_score(jive_simulated_combined, jive_noise_combined_25)
JIVE_pearson_score(jive_simulated_combined, jive_noise_combined_50)

# Compute Tucker congruence coefficients for noise-contaminated data
JIVE_tucker_congruence(jive_simulated_combined, jive_noise_combined_01)
JIVE_tucker_congruence(jive_simulated_combined, jive_noise_combined_05)
JIVE_tucker_congruence(jive_simulated_combined, jive_noise_combined_10)
JIVE_tucker_congruence(jive_simulated_combined, jive_noise_combined_25)
JIVE_tucker_congruence(jive_simulated_combined, jive_noise_combined_50)


#### NOTE: For this, first get outlier-influenced JIVE models from the Outlier_evaluation.R script!
# Plot the results of the outlier-influenced JIVE models
plot(jive_simulated_combined)
plot(jive_outlier_combined_01)
plot(jive_outlier_combined_05)
plot(jive_outlier_combined_10)
plot(jive_outlier_combined_25)
plot(jive_outlier_combined_50)

# Plot the fractions of explained variance for noise-contaminated JIVE models
JIVE_variance_explained(jive_simulated_combined)
JIVE_variance_explained(jive_outlier_combined_01)
JIVE_variance_explained(jive_outlier_combined_05)
JIVE_variance_explained(jive_outlier_combined_10)
JIVE_variance_explained(jive_outlier_combined_25)
JIVE_variance_explained(jive_outlier_combined_50)

# Compute Pearson scores for outlier-contaminated data
JIVE_pearson_score(jive_simulated_combined, jive_outlier_combined_01)
JIVE_pearson_score(jive_simulated_combined, jive_outlier_combined_05)
JIVE_pearson_score(jive_simulated_combined, jive_outlier_combined_10)
JIVE_pearson_score(jive_simulated_combined, jive_outlier_combined_25)
JIVE_pearson_score(jive_simulated_combined, jive_outlier_combined_50)

# Compute Tucker congruence coefficients for outlier-contaminated data
JIVE_tucker_congruence(jive_simulated_combined, jive_outlier_combined_01)
JIVE_tucker_congruence(jive_simulated_combined, jive_outlier_combined_05)
JIVE_tucker_congruence(jive_simulated_combined, jive_outlier_combined_10)
JIVE_tucker_congruence(jive_simulated_combined, jive_outlier_combined_25)
JIVE_tucker_congruence(jive_simulated_combined, jive_outlier_combined_50)

