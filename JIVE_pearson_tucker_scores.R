### Code for obtaining Tucker congruence and Pearson correlation coefficients for JIVE models
# Author: Timo Waling (581706tw)


## Install and load packages
install.packages("r.jive")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("tidyr")
install.packages("tools")
install.packages("clue")

library(r.jive)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tools)
library(clue)


# Set project directory
project_directory <- "F:/Studie/Thesis/Thesis R Project/"

# Function to compute Pearson scores for JIVE models
JIVE_pearson_score <- function(model_1, model_2) {
  # Extract joint components from both JIVE models
  joint_model_1 <- model_1$joint
  joint_model_2 <- model_2$joint
  
  # Ensure that both are matrices
  joint_model_1 <- lapply(joint_model_1, as.matrix)
  joint_model_2 <- lapply(joint_model_2, as.matrix)
  
  # Initialize a matrix to store the Pearson correlations
  n_comp <- length(joint_model_1)
  pearson_scores_matrix <- matrix(0, n_comp, n_comp)
  
  # Compute Pearson correlation coefficients for all possible component pairs
  for (i in 1:n_comp) {
    for (j in 1:n_comp) {
      pearson_scores_matrix[i, j] <- abs(cor(as.vector(joint_model_1[[i]]),
                                             as.vector(joint_model_2[[j]])))
    }
  }
  
  # Find the best component pairings
  best_pairings <- solve_LSAP(pearson_scores_matrix, maximum = TRUE)
  best_component_scores <- pearson_scores_matrix[cbind(1:n_comp, best_pairings)]
  
  # Return the mean of the best pairings
  return(mean(best_component_scores))
}


# Function to compute Tucker Congruence Coefficient for JIVE models
JIVE_tucker_congruence <- function(model_1, model_2) {
  # Extract common components from both JIVE models
  joint_model_1 <- model_1$joint
  joint_model_2 <- model_2$joint
  
  # Ensure that both are matrices
  joint_model_1 <- lapply(joint_model_1, as.matrix)
  joint_model_2 <- lapply(joint_model_2, as.matrix)
  
  # Initialize matrix to store Tucker congruence between each pair of components
  n_comp <- length(joint_model_1)
  tucker_scores_matrix <- matrix(0, nrow = n_comp, ncol = n_comp)
  
  # Compute Tucker's congruence coefficients for all possible component pairs
  for (i in 1:n_comp) {
    vec_1 <- as.vector(joint_model_1[[i]])
    for (j in 1:n_comp) {
      vec_2 <- as.vector(joint_model_2[[j]])
      numerator <- sum(vec_1 * vec_2)^2
      denominator <- sum(vec_1^2) * sum(vec_2^2)
      tucker_scores_matrix[i, j] <- numerator / denominator
    }
  }
  
  # Find the best component pairings
  best_pairings <- solve_LSAP(tucker_scores_matrix, maximum = TRUE)
  best_component_scores <- tucker_scores_matrix[cbind(1:n_comp, best_pairings)]
  
  # Return the mean of the best pairings
  return(mean(best_component_scores))
}

# Function to loop for the multiple realisations of the simulated/contaminated data JIVE models
JIVE_score_loops <- function(n_simulation_datasets, model_directory, file_type) {
  tucker_result_list <- numeric(10)
  pearson_result_list <- numeric(10)
  
  # Loop over all combinations of uncontaminated simulated data and contaminated data JIVE models
  for (i in 1:n_simulation_datasets) {
    filename_uncontaminated <- file.path(project_directory, paste0("Results/Simulated_data/jive_simulated_", i, ".rds"))
    jive_model_uncontaminated <- readRDS(filename_uncontaminated)
    
    filename_contaminated <- file.path(model_directory, paste0(file_type, i, ".rds"))
    jive_model_contaminated <- readRDS(filename_contaminated)
    
    current_tucker <- JIVE_tucker_congruence(jive_model_uncontaminated, jive_model_contaminated)
    tucker_result_list[i] <- current_tucker
    
    current_pearson <- JIVE_pearson_score(jive_model_uncontaminated, jive_model_contaminated)
    pearson_result_list[i] <- current_pearson
  }
  
  # Return the results
  return(list(tucker = tucker_result_list, pearson = pearson_result_list))
}

# Function used in the plotting function below this one
confidence_interval_info <- function(data) {
  data_mean <- mean(data)
  data_standard_dev <- sd(data)
  
  # Obtain confidence intervals for the scores
  error_range <- qt(0.975, df = length(data) - 1) * data_standard_dev / sqrt(length(data))
  lower_bound <- data_mean - error_range
  upper_bound <- data_mean + error_range
  
  # Return the information we need to make the confidence interval error bars
  return(c(mean = data_mean, lower = lower_bound, upper= upper_bound))
}

# Function to visualize the results
plot_scores <- function(scores_matrix) {
  # Obtain confidence interval information for these scores
  score_info <- apply(scores_matrix, 1, confidence_interval_info)
  
  # Transpose the confidence interval information and ensure it is a data frame
  score_info <- t(score_info)
  score_info <- as.data.frame(score_info)
  
  # Order the columns of the scores for the order we want to have in the plot
  score_info$Condition <- rownames(scores_matrix)
  score_info <- score_info[, c("Condition", "mean", "lower", "upper")]
  
  # Separate the scores for Tucker's congruence and Pearson correlation coefficients
  tucker_rows <- grep("Tucker", score_info$Condition, value = TRUE)
  pearson_rows <- grep("Pearson", score_info$Condition, value = TRUE)
  
  # For the plot, we plot the Tucker coefficients on the left and the Pearson coefficients on the right
  score_order <- c(tucker_rows, pearson_rows)
  
  score_info$Condition <- factor(score_info$Condition, levels = score_order)
  
  # We add a vertical line in the plot to separate the Tucker and Pearson scores
  vertical_line_position <- length(tucker_rows) + 0.5
  
  ggplot(score_info, aes(x = Condition, y = mean)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
    geom_vline(xintercept = vertical_line_position, linetype = "dashed", color = "gray", size = 0.8) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(), axis.title.y = element_blank())
}



##### Results for Noise Contamination #####
# Initiate a matrix in which we will store coefficients
n_simulation_datasets <- 10
scores_matrix_noise <- matrix(0, nrow = 10, ncol = 10)
rownames(scores_matrix_noise) <- c("1% Noise Tucker", "5% Noise Tucker", "10% Noise Tucker",
                       "25% Noise Tucker", "50% Noise Tucker", "1% Noise Pearson",
                       "5% Noise Pearson", "10% Noise Pearson", "25% Noise Pearson",
                       "50% Noise Pearson")

## Use the function to fill the scores matrix
# 1% Noise Level
model_directory <- paste0(project_directory, "Results/Noise_Results/Level_01")
file_type <- "jive_noise_01_simulation_"
scores <- JIVE_score_loops(n_simulation_datasets, model_directory, file_type)
scores_matrix_noise[1, ] <- scores$tucker
scores_matrix_noise[6, ] <- scores$pearson

# 5% Noise Level
model_directory <- paste0(project_directory, "Results/Noise_Results/Level_05")
file_type <- "jive_noise_05_simulation_"
scores <- JIVE_score_loops(n_simulation_datasets, model_directory, file_type)
scores_matrix_noise[2, ] <- scores$tucker
scores_matrix_noise[7, ] <- scores$pearson

# 10% Noise Level
model_directory <- paste0(project_directory, "Results/Noise_Results/Level_10")
file_type <- "jive_noise_10_simulation_"
scores <- JIVE_score_loops(n_simulation_datasets, model_directory, file_type)
scores_matrix_noise[3, ] <- scores$tucker
scores_matrix_noise[8, ] <- scores$pearson

# 25% Noise Level
model_directory <- paste0(project_directory, "Results/Noise_Results/Level_25")
file_type <- "jive_noise_25_simulation_"
scores <- JIVE_score_loops(n_simulation_datasets, model_directory, file_type)
scores_matrix_noise[4, ] <- scores$tucker
scores_matrix_noise[9, ] <- scores$pearson

# 50% Noise Level
model_directory <- paste0(project_directory, "Results/Noise_Results/Level_50")
file_type <- "jive_noise_50_simulation_"
scores <- JIVE_score_loops(n_simulation_datasets, model_directory, file_type)
scores_matrix_noise[5, ] <- scores$tucker
scores_matrix_noise[10, ] <- scores$pearson





##### Results for Outlier Contamination #####
# Initiate a matrix in which we will store coefficients
n_simulation_datasets <- 10
scores_matrix_outlier <- matrix(0, nrow = 10, ncol = 10)
rownames(scores_matrix_outlier) <- c("1% Outlier Tucker", "5% Outlier Tucker", "10% Outlier Tucker",
                             "25% Outlier Tucker", "50% Outlier Tucker", "1% Outlier Pearson",
                             "5% Outlier Pearson", "10% Outlier Pearson", "25% Outlier Pearson",
                             "50% Outlier Pearson")

## Use the function to fill the scores matrix
# 1% Outlier Level
model_directory <- paste0(project_directory, "Results/Outlier_Results/Level_01")
file_type <- "jive_outlier_01_simulation_"
scores <- JIVE_score_loops(n_simulation_datasets, model_directory, file_type)
scores_matrix_outlier[1, ] <- scores$tucker
scores_matrix_outlier[6, ] <- scores$pearson

# 5% Outlier Level
model_directory <- paste0(project_directory, "Results/Outlier_Results/Level_05")
file_type <- "jive_outlier_05_simulation_"
scores <- JIVE_score_loops(n_simulation_datasets, model_directory, file_type)
scores_matrix_outlier[2, ] <- scores$tucker
scores_matrix_outlier[7, ] <- scores$pearson

# 10% Outlier Level
model_directory <- paste0(project_directory, "Results/Outlier_Results/Level_10")
file_type <- "jive_outlier_10_simulation_"
scores <- JIVE_score_loops(n_simulation_datasets, model_directory, file_type)
scores_matrix_outlier[3, ] <- scores$tucker
scores_matrix_outlier[8, ] <- scores$pearson

# 25% Outlier Level
model_directory <- paste0(project_directory, "Results/Outlier_Results/Level_25")
file_type <- "jive_outlier_25_simulation_"
scores <- JIVE_score_loops(n_simulation_datasets, model_directory, file_type)
scores_matrix_outlier[4, ] <- scores$tucker
scores_matrix_outlier[9, ] <- scores$pearson

# 50% Outlier Level
model_directory <- paste0(project_directory, "Results/Outlier_Results/Level_50")
file_type <- "jive_outlier_50_simulation_"
scores <- JIVE_score_loops(n_simulation_datasets, model_directory, file_type)
scores_matrix_outlier[5, ] <- scores$tucker
scores_matrix_outlier[10, ] <- scores$pearson

# Print the means of the rows of the scores matrices
print("Row Means: JIVE Noise Contamination:")
print(rowMeans(scores_matrix_noise))

print("Row Means: JIVE Outlier Contamination:")
print(rowMeans(scores_matrix_outlier))

# Visualize the results
plot_scores(scores_matrix_noise)
plot_scores(scores_matrix_outlier)
