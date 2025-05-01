### Code For Adding Outliers To Simulated Data for Quantitative Analysis
# Author: Timo Waling (581706tw)

# Install relevant packages if not already installed
install.packages("gridExtra")
install.packages("pheatmap")
# Load relevant packages
library(gridExtra)
library(pheatmap)


#### Functions used
# Function to add outliers
apply_outlier_to_data <- function(data, fraction = 0.05, outlier_strength = 10) {
  # Create a copy in order to preserve the original input data and ensure data is a matrix
  data_copy <- as.matrix(data)
  # Calculate statistics from data
  n_rows <- nrow(data_copy)
  n_cols <- ncol(data_copy)
  sd_from_data <- mean(apply(data_copy, 2, sd))
  # Pick a fraction of rows to contaminate
  n_outlier_rows <- ceiling(n_rows * fraction)
  outlier_row_indices <- sample(seq_len(n_rows), size = n_outlier_rows, replace = FALSE)
  # Loop over columns
  for (col_idx in seq_len(n_cols)) {
    # Add strong outliers to selected rows
    noise <- rnorm(length(outlier_row_indices), mean = 0, sd = outlier_strength * sd_from_data)
    data_copy[outlier_row_indices, col_idx] <- data_copy[outlier_row_indices, col_idx] + noise
  }
  return(data_copy)
}
# Create a function that will fit the input data into a DISCO-SCA model
create_DISCOSCA_model <- function(data_reho, data_lfcd, data_alff, n_comp) {
  # Scale data blocks and put them into one list
  data_list <- list(scale(as.matrix(data_reho)),
                    scale(as.matrix(data_lfcd)),
                    scale(as.matrix(data_alff)))
  
  # Create the DISCO-SCA model with the given data and given number of components
  result <- disco(data_list, ncomp = n_comp)
  
  # Return the results of the DISCO-SCA model
  return(result)
}


#### Load in the original data
# Loading in Reho data
data_reho_combined <- read.csv("F:/Studie/Thesis/Thesis R Project/data/reho_data/reho_atlas_combined.csv", header = FALSE)
data_reho_combined <- data_reho_combined[-1, ] #Remove first row of 0's

# Loading in Lfcd data
data_lfcd_combined <- read.csv("F:/Studie/Thesis/Thesis R Project/data/lfcd_data/lfcd_atlas_combined.csv", header = FALSE)
data_lfcd_combined <- data_lfcd_combined[-1, ] #Remove first row of 0's

# Loading in ALFF data
data_alff_combined <- read.csv("F:/Studie/Thesis/Thesis R Project/data/alff_data/alff_atlas_combined.csv", header = FALSE)
data_alff_combined <- data_alff_combined[-1, ] #Remove first row of 0's

#### Load in the Simulated  data (see: Simulation_Study.R for the generation of this data)
simulated_data_combined <- read.csv("F:/Studie/Thesis/Thesis R Project/data/simulated_data/simulated_data_combined.csv", header = FALSE)
simulated_data_combined <- simulated_data_combined[-1,]

# Simulation data is still Characters instead of numeric, fix it
simulated_data_combined[] <- lapply(simulated_data_combined, function(x) as.numeric(as.character(x)))

# Extract the Reho, Lfcd and ALFF data from the simulated sets
simulated_reho_combined <- simulated_data_combined[, 1:111]
simulated_lfcd_combined <- simulated_data_combined[, 112:222]
simulated_alff_combined <- simulated_data_combined[, 223:333]

#### Add outlier to the data for different outlier levels (percentages)
# For simulated Reho data:
outlier_reho_combined_01 <- apply_outlier_to_data(simulated_reho_combined, 0.01, 7)
outlier_reho_combined_05 <- apply_outlier_to_data(simulated_reho_combined, 0.05, 7)
outlier_reho_combined_10 <- apply_outlier_to_data(simulated_reho_combined, 0.1, 7)
outlier_reho_combined_25 <- apply_outlier_to_data(simulated_reho_combined, 0.25, 7)
outlier_reho_combined_50 <- apply_outlier_to_data(simulated_reho_combined, 0.5, 7)

# For simulated Lfcd data:
outlier_lfcd_combined_01 <- apply_outlier_to_data(simulated_lfcd_combined, 0.01, 7)
outlier_lfcd_combined_05 <- apply_outlier_to_data(simulated_lfcd_combined, 0.05, 7)
outlier_lfcd_combined_10 <- apply_outlier_to_data(simulated_lfcd_combined, 0.1, 7)
outlier_lfcd_combined_25 <- apply_outlier_to_data(simulated_lfcd_combined, 0.25, 7)
outlier_lfcd_combined_50 <- apply_outlier_to_data(simulated_lfcd_combined, 0.5, 7)

# For simulated ALFF data:
outlier_alff_combined_01 <- apply_outlier_to_data(simulated_alff_combined, 0.01, 7)
outlier_alff_combined_05 <- apply_outlier_to_data(simulated_alff_combined, 0.05, 7)
outlier_alff_combined_10 <- apply_outlier_to_data(simulated_alff_combined, 0.1, 7)
outlier_alff_combined_25 <- apply_outlier_to_data(simulated_alff_combined, 0.25, 7)
outlier_alff_combined_50 <- apply_outlier_to_data(simulated_alff_combined, 0.5, 7)

#### Plot contaminated data
# Plot the correlation heatmaps for the outlier-influenced data
cor_hm_out_combined_0 <- pheatmap(cor(cbind(simulated_reho_combined, simulated_lfcd_combined, simulated_alff_combined)),
                              cluster_rows = FALSE, cluster_cols = FALSE,
                              labels_row = "", labels_col = "", main = "0% outlier", silent = TRUE)
cor_hm_out_combined_01 <- pheatmap(cor(cbind(outlier_reho_combined_01, outlier_lfcd_combined_01, outlier_alff_combined_01)),
                               cluster_rows = FALSE, cluster_cols = FALSE,
                               labels_row = "", labels_col = "", main = "1% outlier", silent = TRUE)
cor_hm_out_combined_05 <- pheatmap(cor(cbind(outlier_reho_combined_05, outlier_lfcd_combined_05, outlier_alff_combined_05)),
                               cluster_rows = FALSE, cluster_cols = FALSE,
                               labels_row = "", labels_col = "", main = "5% outlier", silent = TRUE)
cor_hm_out_combined_10 <- pheatmap(cor(cbind(outlier_reho_combined_10, outlier_lfcd_combined_10, outlier_alff_combined_10)),
                               cluster_rows = FALSE, cluster_cols = FALSE,
                               labels_row = "", labels_col = "", main = "10% outlier", silent = TRUE)
cor_hm_out_combined_25 <- pheatmap(cor(cbind(outlier_reho_combined_25, outlier_lfcd_combined_25, outlier_alff_combined_25)),
                               cluster_rows = FALSE, cluster_cols = FALSE,
                               labels_row = "", labels_col = "", main = "25% outlier", silent = TRUE)
cor_hm_out_combined_50 <- pheatmap(cor(cbind(outlier_reho_combined_50, outlier_lfcd_combined_50, outlier_alff_combined_50)),
                               cluster_rows = FALSE, cluster_cols = FALSE,
                               labels_row = "", labels_col = "", main = "50% outlier", silent = TRUE)

# Plot the heatmaps in a 3 by 3 grid for easy comparison
grid.arrange(cor_hm_out_combined_0$gtable, cor_hm_out_combined_01$gtable, cor_hm_out_combined_05$gtable,
             cor_hm_out_combined_10$gtable, cor_hm_out_combined_25$gtable, cor_hm_out_combined_50$gtable, ncol = 3)

#### DISCO-SCA
# Set the number of components to be used for DISCO-SCA
n_comp = 3
# Perform DISCO-SCA for the outlier contaminated data
discosca_result_outlier_combined_01 <- create_DISCOSCA_model(outlier_reho_combined_01, outlier_lfcd_combined_01,
                                                           outlier_alff_combined_01, n_comp)
discosca_result_outlier_combined_05 <- create_DISCOSCA_model(outlier_reho_combined_05, outlier_lfcd_combined_05,
                                                           outlier_alff_combined_05, n_comp)
discosca_result_outlier_combined_10 <- create_DISCOSCA_model(outlier_reho_combined_10, outlier_lfcd_combined_10,
                                                           outlier_alff_combined_10, n_comp)
discosca_result_outlier_combined_25 <- create_DISCOSCA_model(outlier_reho_combined_25, outlier_lfcd_combined_25,
                                                           outlier_alff_combined_25, n_comp)
discosca_result_outlier_combined_50 <- create_DISCOSCA_model(outlier_reho_combined_50, outlier_lfcd_combined_50,
                                                           outlier_alff_combined_50, n_comp)

#### JIVE
# Create datablocks that will be used to make JIVE models
data_blocks_outlier_combined_01 <- list(reho = t(outlier_reho_combined_01), lfcd = t(outlier_lfcd_combined_01), 
                                      alff = t(outlier_alff_combined_01))
data_blocks_outlier_combined_05 <- list(reho = t(outlier_reho_combined_05), lfcd = t(outlier_lfcd_combined_05), 
                                      alff = t(outlier_alff_combined_05))
data_blocks_outlier_combined_10 <- list(reho = t(outlier_reho_combined_10), lfcd = t(outlier_lfcd_combined_10), 
                                      alff = t(outlier_alff_combined_10))
data_blocks_outlier_combined_25 <- list(reho = t(outlier_reho_combined_25), lfcd = t(outlier_lfcd_combined_25), 
                                      alff = t(outlier_alff_combined_25))
data_blocks_outlier_combined_50 <- list(reho = t(outlier_reho_combined_50), lfcd = t(outlier_lfcd_combined_50), 
                                      alff = t(outlier_alff_combined_50))
# Create JIVE models on the outlier data:
jive_outlier_combined_01 <- jive(data_blocks_outlier_combined_01, method = "given", rankJ = 3, rankA = c(2, 2, 2))
jive_outlier_combined_05 <- jive(data_blocks_outlier_combined_05, method = "given", rankJ = 3, rankA = c(2, 2, 2))
jive_outlier_combined_10 <- jive(data_blocks_outlier_combined_10, method = "given", rankJ = 3, rankA = c(2, 2, 2))
jive_outlier_combined_25 <- jive(data_blocks_outlier_combined_25, method = "given", rankJ = 3, rankA = c(2, 2, 2))
jive_outlier_combined_50 <- jive(data_blocks_outlier_combined_50, method = "given", rankJ = 3, rankA = c(2, 2, 2))


