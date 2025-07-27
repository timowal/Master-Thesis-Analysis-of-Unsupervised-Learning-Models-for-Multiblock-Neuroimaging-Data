### Code For Adding Noise and Outlier Contamination To Simulated Data for Multiple Simulated Data Realizations
# Author: Timo Waling (581706tw)


## Install packages and load them:
install.packages("gridExtra")
install.packages("pheatmap")
install.packages("grid")

library(gridExtra)
library(pheatmap)
library(grid)


## Load in the original data
# Loading in REHO data
data_reho_combined <- read.csv("F:/Studie/Thesis/Thesis R Project/data/reho_data/reho_atlas_combined.csv", header = FALSE)
data_reho_combined <- data_reho_combined[-1, ] #Remove first row of 0's

# Loading in LFCD data
data_lfcd_combined <- read.csv("F:/Studie/Thesis/Thesis R Project/data/lfcd_data/lfcd_atlas_combined.csv", header = FALSE)
data_lfcd_combined <- data_lfcd_combined[-1, ] #Remove first row of 0's

# Loading in ALFF data
data_alff_combined <- read.csv("F:/Studie/Thesis/Thesis R Project/data/alff_data/alff_atlas_combined.csv", header = FALSE)
data_alff_combined <- data_alff_combined[-1, ] #Remove first row of 0's


##### Noise Contamination #####
# Create width levels based on the original data, used for noise contamination
width_reho = (max(data_reho_combined) - min(data_reho_combined))
width_lfcd = (max(data_lfcd_combined) - min(data_lfcd_combined))
width_alff = (max(data_alff_combined) - min(data_alff_combined))

### Function to add noise based on the width of data
apply_noise_to_data <- function(data, width, noise_level) {
  data_copy <- data
  
  # Simulate random noise based on width and noise level
  noise <- matrix(rnorm(length(data_copy), mean = 0, sd = width * noise_level),
                  nrow = nrow(data_copy), ncol = ncol(data_copy))
  
  # Add the noise to the data
  result <- data_copy + noise
  return(result)
}


### Loop through the simulated datasets and add different levels of noise contamination
set.seed(2025)
n_simulation_datasets <- 10
data_directory <- "F:/Studie/Thesis/Thesis R Project/Simulated_data/Uncontaminated"

for (i in 1:n_simulation_datasets) {
  cat("Processing dataset", i, "...\n")
  
  filename <- file.path(data_directory, paste0("simulated_data_", i, ".csv"))
  simulated_data <- read.csv(filename, header = FALSE)[-1, ]
  
  # Ensure the data is read as numeric data
  simulated_data[] <- lapply(simulated_data, function(x) as.numeric(as.character(x)))
  
  # Extract the REHO, LFCD and ALFF data from the simulated sets
  simulated_reho <- simulated_data[, 1:111]
  simulated_lfcd <- simulated_data[, 112:222]
  simulated_alff <- simulated_data[, 223:333]
  
  ## Add noise to the data for different noise levels
  # For simulated REHO data:
  noise_reho_01 <- apply_noise_to_data(simulated_reho, width_reho, 0.01)
  noise_reho_05 <- apply_noise_to_data(simulated_reho, width_reho, 0.05)
  noise_reho_10 <- apply_noise_to_data(simulated_reho, width_reho, 0.1)
  noise_reho_25 <- apply_noise_to_data(simulated_reho, width_reho, 0.25)
  noise_reho_50 <- apply_noise_to_data(simulated_reho, width_reho, 0.5)
  
  # For simulated LFCD data:
  noise_lfcd_01 <- apply_noise_to_data(simulated_lfcd, width_lfcd, 0.01)
  noise_lfcd_05 <- apply_noise_to_data(simulated_lfcd, width_lfcd, 0.05)
  noise_lfcd_10 <- apply_noise_to_data(simulated_lfcd, width_lfcd, 0.1)
  noise_lfcd_25 <- apply_noise_to_data(simulated_lfcd, width_lfcd, 0.25)
  noise_lfcd_50 <- apply_noise_to_data(simulated_lfcd, width_lfcd, 0.5)
  
  # For simulated ALFF data:
  noise_alff_01 <- apply_noise_to_data(simulated_alff, width_alff, 0.01)
  noise_alff_05 <- apply_noise_to_data(simulated_alff, width_alff, 0.05)
  noise_alff_10 <- apply_noise_to_data(simulated_alff, width_alff, 0.1)
  noise_alff_25 <- apply_noise_to_data(simulated_alff, width_alff, 0.25)
  noise_alff_50 <- apply_noise_to_data(simulated_alff, width_alff, 0.5)
  
  
  # Plot the correlation heatmaps for the noise contaminated data
  cor_hm_0 <- pheatmap(cor(cbind(simulated_reho, simulated_lfcd, simulated_alff)),
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       labels_row = "", labels_col = "", main = "0% Noise", silent = TRUE)
  
  cor_hm_01 <- pheatmap(cor(cbind(noise_reho_01, noise_lfcd_01, noise_alff_01)),
                        cluster_rows = FALSE, cluster_cols = FALSE,
                        labels_row = "", labels_col = "", main = "1% Noise", silent = TRUE)
  
  cor_hm_05 <- pheatmap(cor(cbind(noise_reho_05, noise_lfcd_05, noise_alff_05)),
                        cluster_rows = FALSE, cluster_cols = FALSE,
                        labels_row = "", labels_col = "", main = "5% Noise", silent = TRUE)
  
  cor_hm_10 <- pheatmap(cor(cbind(noise_reho_10, noise_lfcd_10, noise_alff_10)),
                        cluster_rows = FALSE, cluster_cols = FALSE,
                        labels_row = "", labels_col = "", main = "10% Noise", silent = TRUE)
  
  cor_hm_25 <- pheatmap(cor(cbind(noise_reho_25, noise_lfcd_25, noise_alff_25)),
                        cluster_rows = FALSE, cluster_cols = FALSE,
                        labels_row = "", labels_col = "", main = "25% Noise", silent = TRUE)
  
  cor_hm_50 <- pheatmap(cor(cbind(noise_reho_50, noise_lfcd_50, noise_alff_50)),
                        cluster_rows = FALSE, cluster_cols = FALSE,
                        labels_row = "", labels_col = "", main = "50% Noise", silent = TRUE)
  
  # Plot the heatmaps
  grid.arrange(cor_hm_0$gtable, cor_hm_01$gtable, cor_hm_05$gtable,
               cor_hm_10$gtable, cor_hm_25$gtable, cor_hm_50$gtable, ncol = 3)
  
  ## Save the Contaminated data of this simulation dataset for all noise levels
  # Noise level 01
  save_directory <- "F:/Studie/Thesis/Thesis R Project/Simulated_data/Noise/Level_01"
  combined_noise_01 <- cbind(noise_reho_01, noise_lfcd_01, noise_alff_01)
  filename <- file.path(save_directory, paste0("noise_data_01_simulation_", i, ".csv"))
  write.csv(combined_noise_01, file = filename, row.names = FALSE)
  
  # Noise level 05
  save_directory <- "F:/Studie/Thesis/Thesis R Project/Simulated_data/Noise/Level_05"
  combined_noise_05 <- cbind(noise_reho_05, noise_lfcd_05, noise_alff_05)
  filename <- file.path(save_directory, paste0("noise_data_05_simulation_", i, ".csv"))
  write.csv(combined_noise_05, file = filename, row.names = FALSE)
  
  # Noise level 10
  save_directory <- "F:/Studie/Thesis/Thesis R Project/Simulated_data/Noise/Level_10"
  combined_noise_10 <- cbind(noise_reho_10, noise_lfcd_10, noise_alff_10)
  filename <- file.path(save_directory, paste0("noise_data_10_simulation_", i, ".csv"))
  write.csv(combined_noise_10, file = filename, row.names = FALSE)
  
  # Noise level 25
  save_directory <- "F:/Studie/Thesis/Thesis R Project/Simulated_data/Noise/Level_25"
  combined_noise_25 <- cbind(noise_reho_25, noise_lfcd_25, noise_alff_25)
  filename <- file.path(save_directory, paste0("noise_data_25_simulation_", i, ".csv"))
  write.csv(combined_noise_25, file = filename, row.names = FALSE)
  
  # Noise level 50
  save_directory <- "F:/Studie/Thesis/Thesis R Project/Simulated_data/Noise/Level_50"
  combined_noise_50 <- cbind(noise_reho_50, noise_lfcd_50, noise_alff_50)
  filename <- file.path(save_directory, paste0("noise_data_50_simulation_", i, ".csv"))
  write.csv(combined_noise_50, file = filename, row.names = FALSE)
}




##### Outlier Contamination #####
### Function for Adding Outlier Contamination
# Function to add outliers
apply_outlier_to_data <- function(data, fraction = 0.05, outlier_strength = 10) {
  data_copy <- as.matrix(data)
  
  n_rows <- nrow(data_copy)
  n_cols <- ncol(data_copy)
  sd_from_data <- mean(apply(data_copy, 2, sd))
  
  # Pick a fraction of rows to contaminate
  n_outlier_rows <- ceiling(n_rows * fraction)
  outlier_row_indices <- sample(seq_len(n_rows), size = n_outlier_rows, replace = FALSE)
  
  # Loop over columns
  for (i in seq_len(n_cols)) {
    # Add strong outliers to selected rows
    contamination <- rnorm(length(outlier_row_indices), mean = 0, sd = outlier_strength * sd_from_data)
    
    data_copy[outlier_row_indices, i] <- data_copy[outlier_row_indices, i] + contamination
  }
  
  return(data_copy)
}


### Loop through the simulated datasets and add different levels of outlier contamination
set.seed(2025)
n_simulation_datasets <- 10
data_directory <- "F:/Studie/Thesis/Thesis R Project/Simulated_data/Uncontaminated"

for (i in 1:n_simulation_datasets) {
  cat("Processing dataset", i, "...\n")
  
  filename <- file.path(data_directory, paste0("simulated_data_", i, ".csv"))
  simulated_data <- read.csv(filename, header = FALSE)[-1, ]
  
  # Ensure the read data is numeric
  simulated_data[] <- lapply(simulated_data, function(x) as.numeric(as.character(x)))
  
  # Extract the REHO, LFCD and ALFF data from the simulated sets
  simulated_reho <- simulated_data[, 1:111]
  simulated_lfcd <- simulated_data[, 112:222]
  simulated_alff <- simulated_data[, 223:333]
  
  ## Add outlier to the data for different outlier levels
  # For simulated REHO data:
  outlier_reho_01 <- apply_outlier_to_data(simulated_reho, 0.01, 7)
  outlier_reho_05 <- apply_outlier_to_data(simulated_reho, 0.05, 7)
  outlier_reho_10 <- apply_outlier_to_data(simulated_reho, 0.1, 7)
  outlier_reho_25 <- apply_outlier_to_data(simulated_reho, 0.25, 7)
  outlier_reho_50 <- apply_outlier_to_data(simulated_reho, 0.5, 7)
  
  # For simulated LFCD data:
  outlier_lfcd_01 <- apply_outlier_to_data(simulated_lfcd, 0.01, 7)
  outlier_lfcd_05 <- apply_outlier_to_data(simulated_lfcd, 0.05, 7)
  outlier_lfcd_10 <- apply_outlier_to_data(simulated_lfcd, 0.1, 7)
  outlier_lfcd_25 <- apply_outlier_to_data(simulated_lfcd, 0.25, 7)
  outlier_lfcd_50 <- apply_outlier_to_data(simulated_lfcd, 0.5, 7)
  
  # For simulated ALFF data:
  outlier_alff_01 <- apply_outlier_to_data(simulated_alff, 0.01, 7)
  outlier_alff_05 <- apply_outlier_to_data(simulated_alff, 0.05, 7)
  outlier_alff_10 <- apply_outlier_to_data(simulated_alff, 0.1, 7)
  outlier_alff_25 <- apply_outlier_to_data(simulated_alff, 0.25, 7)
  outlier_alff_50 <- apply_outlier_to_data(simulated_alff, 0.5, 7)
  
  
  # Plot the correlation heatmaps for the outlier contaminated data
  cor_hm_out_0 <- pheatmap(cor(cbind(simulated_reho, simulated_lfcd, simulated_alff)),
                           cluster_rows = FALSE, cluster_cols = FALSE,
                           labels_row = "", labels_col = "", main = "0% outlier", silent = TRUE)
  
  cor_hm_out_01 <- pheatmap(cor(cbind(outlier_reho_01, outlier_lfcd_01, outlier_alff_01)),
                            cluster_rows = FALSE, cluster_cols = FALSE,
                            labels_row = "", labels_col = "", main = "1% outlier", silent = TRUE)
  
  cor_hm_out_05 <- pheatmap(cor(cbind(outlier_reho_05, outlier_lfcd_05, outlier_alff_05)),
                            cluster_rows = FALSE, cluster_cols = FALSE,
                            labels_row = "", labels_col = "", main = "5% outlier", silent = TRUE)
  
  cor_hm_out_10 <- pheatmap(cor(cbind(outlier_reho_10, outlier_lfcd_10, outlier_alff_10)),
                            cluster_rows = FALSE, cluster_cols = FALSE,
                            labels_row = "", labels_col = "", main = "10% outlier", silent = TRUE)
  
  cor_hm_out_25 <- pheatmap(cor(cbind(outlier_reho_25, outlier_lfcd_25, outlier_alff_25)),
                            cluster_rows = FALSE, cluster_cols = FALSE,
                            labels_row = "", labels_col = "", main = "25% outlier", silent = TRUE)
  
  cor_hm_out_50 <- pheatmap(cor(cbind(outlier_reho_50, outlier_lfcd_50, outlier_alff_50)),
                            cluster_rows = FALSE, cluster_cols = FALSE,
                            labels_row = "", labels_col = "", main = "50% outlier", silent = TRUE)
  
  # Plot the heatmaps
  grid.arrange(cor_hm_out_0$gtable, cor_hm_out_01$gtable, cor_hm_out_05$gtable,
               cor_hm_out_10$gtable, cor_hm_out_25$gtable, cor_hm_out_50$gtable, ncol = 3)
  
  ## Save the Contaminated data of this simulation dataset for all outlier levels
  # Outlier level 01
  save_directory <- "F:/Studie/Thesis/Thesis R Project/Simulated_data/Outlier/Level_01"
  combined_outlier_01 <- cbind(outlier_reho_01, outlier_lfcd_01, outlier_alff_01)
  filename <- file.path(save_directory, paste0("outlier_data_01_simulation_", i, ".csv"))
  write.csv(combined_outlier_01, file = filename, row.names = FALSE)
  
  # Outlier level 05
  save_directory <- "F:/Studie/Thesis/Thesis R Project/Simulated_data/Outlier/Level_05"
  combined_outlier_05 <- cbind(outlier_reho_05, outlier_lfcd_05, outlier_alff_05)
  filename <- file.path(save_directory, paste0("outlier_data_05_simulation_", i, ".csv"))
  write.csv(combined_outlier_05, file = filename, row.names = FALSE)
  
  # Outlier level 10
  save_directory <- "F:/Studie/Thesis/Thesis R Project/Simulated_data/Outlier/Level_10"
  combined_outlier_10 <- cbind(outlier_reho_10, outlier_lfcd_10, outlier_alff_10)
  filename <- file.path(save_directory, paste0("outlier_data_10_simulation_", i, ".csv"))
  write.csv(combined_outlier_10, file = filename, row.names = FALSE)
  
  # Outlier level 25
  save_directory <- "F:/Studie/Thesis/Thesis R Project/Simulated_data/Outlier/Level_25"
  combined_outlier_25 <- cbind(outlier_reho_25, outlier_lfcd_25, outlier_alff_25)
  filename <- file.path(save_directory, paste0("outlier_data_25_simulation_", i, ".csv"))
  write.csv(combined_outlier_25, file = filename, row.names = FALSE)
  
  # Outlier level 50
  save_directory <- "F:/Studie/Thesis/Thesis R Project/Simulated_data/Outlier/Level_50"
  combined_outlier_50 <- cbind(outlier_reho_50, outlier_lfcd_50, outlier_alff_50)
  filename <- file.path(save_directory, paste0("outlier_data_50_simulation_", i, ".csv"))
  write.csv(combined_outlier_50, file = filename, row.names = FALSE)
}
