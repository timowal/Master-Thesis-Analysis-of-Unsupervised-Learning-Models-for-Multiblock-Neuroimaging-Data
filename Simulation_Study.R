### Code For Generating Data For Simulation Study
# Author: Timo Waling (581706tw)


# Install relevant packages if not already installed
install.packages("ggplot2")
install.packages("stats")
install.packages("fitdistrplus")
install.packages("copula")
install.packages("dplyr")
install.packages("pheatmap")
install.packages("VineCopula")
install.packages("statmod")
# Load relevant packages
library(ggplot2)
library(stats)
library(fitdistrplus)
library(copula)
library(dplyr)
library(pheatmap)
library(VineCopula)
library(statmod)

#### Functions used
# Function that takes data of the data block for one ROI, and fits the data to a
# Gamma, Exponential, and Weibull distribution in order to find
# which model best fits the data block
fit_data_to_distributions <- function(data) {
  #### Gamma distribution:
  # Fit the data onto a Gamma distribution:
  fitted_gamma_model <- fitdist(data, "gamma", method = "mle")
  # Get log-likelihood, AIC and BIC:
  fitted_gamma_logl <- fitted_gamma_model$loglik
  fitted_gamma_aic <- fitted_gamma_model$aic
  fitted_gamma_bic <- fitted_gamma_model$bic
  
  #### Exponential distribution
  fitted_exponential_model <- fitdist(data, "exp", method = "mle")
  # Get log-likelihood, AIC and BIC:
  fitted_exponential_logl <- fitted_exponential_model$loglik
  fitted_exponential_aic <- fitted_exponential_model$aic
  fitted_exponential_bic <- fitted_exponential_model$bic
  
  #### Weibull distribution
  # Fit the data onto a Weibull distribution:
  fitted_weibull_model <- fitdist(data, "weibull", method = "mle")
  # Get log-likelihood, AIC and BIC:
  fitted_weibull_logl <- fitted_weibull_model$loglik
  fitted_weibull_aic <- fitted_weibull_model$aic
  fitted_weibull_bic <- fitted_weibull_model$bic
  
  # Concatenate the log-likelihood, AIC and BIC values
  values_logl <- c(fitted_gamma_logl, fitted_exponential_logl, fitted_weibull_logl)
  values_aic <- c(fitted_gamma_aic, fitted_exponential_aic, fitted_weibull_aic)
  values_bic <- c(fitted_gamma_bic, fitted_exponential_bic, fitted_weibull_bic)
  
  # Obtain the index corresponding to the distribution with the best value
  best_logl_index <- which.max(values_logl) # For log-likelihood, higher values are preferred
  best_aic_index <- which.min(values_aic) # For AIC, lower values are preferred
  best_bic_index <- which.min(values_bic) # For BIC, lower values are preferred
  
  # Return the indices of the best fitting distributions
  result <- list(best_logl_index = best_logl_index, best_aic_index = best_aic_index, best_bic_index = best_bic_index)
  return(result)
}

# Function that loops through the ROIs to allow user to evaluate which distribution best describes said ROI
investigate_ROIs_for_data <- function(data) {
  # Create vectors to save best distribution choices per criterion
  counter_logl_best <- c("Gamma" = 0, "Exponential" = 0, "Weibull" = 0)
  counter_aic_best <- c("Gamma" = 0, "Exponential" = 0, "Weibull" = 0)
  counter_bic_best <- c("Gamma" = 0, "Exponential" = 0, "Weibull" = 0)
  
  # Loop through all 111 ROIs and fit the data to the distributions
  for (i in 1:111) {
    #Get the data for the i-th ROI
    investigate_ROI <- data[, i]
    # Get the log-likelihood, AIC and BIC values for fitting this ROI data onto the distributions
    fitted_distribution_scores <- fit_data_to_distributions(investigate_ROI)
    # Check which distribution had the highest log-likelihood value
    best_logl_dist <- c("Gamma", "Exponential", "Weibull")[fitted_distribution_scores$best_logl_index]
    best_aic_dist <- c("Gamma", "Exponential", "Weibull")[fitted_distribution_scores$best_aic_index]
    best_bic_dist <- c("Gamma", "Exponential", "Weibull")[fitted_distribution_scores$best_bic_index]
    # Update the counters for how many times a certain distribution provided the best fit for the ROI
    counter_logl_best[best_logl_dist] <- counter_logl_best[best_logl_dist] + 1
    counter_aic_best[best_aic_dist] <- counter_aic_best[best_aic_dist] + 1
    counter_bic_best[best_bic_dist] <- counter_bic_best[best_bic_dist] + 1
  }
  # After going through all the ROIs, print the result for which distribution was best, and how often
  table_of_results <- data.frame(Distribution = c("Gamma", "Exponential", "Weibull"),
                                 amount_logl = counter_logl_best,
                                 amount_aic = counter_aic_best,
                                 amount_bic = counter_bic_best)
  print(table_of_results)
}
#### Functions for simulating data

# Function calculates certain required data characteristics necessary for transformations
fitted_marginals_for_data <- function(data, distribution_type = "gamma") {
  # Create an empty list to return the marginal distributions of the data
  result <- list()
  # Iterate through the columns (ROIs) for the data and fit the ROI data to the given distribution
  for (i in 1:ncol(data)) {
    ROI_fit <- fitdist(data[, i], distribution_type, method = "mle")
    result[[i]] <- ROI_fit
  }
  return(result)
}

# Function that turns the simulated uniform pseudo data back to the original marginals
transform_pseudo_to_original_marginals <- function(simulated_data, marginals, distribution_type = "gamma") {
  # Create an empty matrix to return the transformed data in
  result <- matrix(NA, nrow(simulated_data), length(marginals))
  # Loop through the ROIs
  for (i in 1:ncol(simulated_data)) {
    #### Manner of transformation is dependent on distribution type
    # For Gamma distribution:
    if (distribution_type == "gamma") {
      param_shape <- marginals[[i]]$estimate[1] # Shape parameter for the Gamma function
      param_rate <- marginals[[i]]$estimate[2] # Rate parameter for the Gamma function
      # Transformation using qgamma
      result[, i] <- qgamma(simulated_data[, i], shape = param_shape, rate = param_rate)
    }
    # For Exponential distribution:
    if (distribution_type == "exp") {
      param_rate <- marginals[[i]]$estimate[1] # Rate parameter for Exponential distribution
      # Transformation using qlnorm
      result[, i] <- qexp(simulated_data[, i], rate = param_rate)
    }
    # For Weibull distribution:
    if (distribution_type == "weibull") {
      param_shape <- marginals[[i]]$estimate[1] # Shape parameter for the Weibull function
      param_scale <- marginals[[i]]$estimate[2] # Scale parameter for the Weibull function
      # Transformation using qweibull
      result[, i] <- qweibull(simulated_data[, i], shape = param_shape, scale = param_scale)
    }
  }
  return(result)
}

# Function to check correlation for original and simulated data for a certain ROI
check_correlation_ROI <- function(original_data, simulated_data) {
  # Check correlation for both simulated and the original data
  correlation_original <- cor(original_data)
  correlation_simulated <- cor(simulated_data)
  # Print correlation results
  print("Correlation Matrix")
  print("Correlation for this ROI for Original Data:")
  print(correlation_original)
  print("Correlation for this ROI for Simulated Data:")
  print(correlation_simulated)
}

# Use a seed for consistency in results
set.seed(2025)

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

#### Distributions require strictly positive values, so replace 0 with 1*10^(-18)
# For Reho
data_reho_combined[data_reho_combined == 0] <- 1e-18
# For Lfcd
data_lfcd_combined[data_lfcd_combined == 0] <- 1e-18
# For ALFF
data_alff_combined[data_alff_combined == 0] <- 1e-18


#### Data insights
# Plot histograms to see distribution of a single ROI
# For REHO
hist(data_reho_combined[[1]], breaks=25, main = "Histogram of Reho data for ROI 1", xlab = "Value")
# For LFCD
hist(data_lfcd_combined[[1]], breaks=25, main = "Histogram of Lfcd data for ROI 1", xlab = "Value")
# For ALFF
hist(data_alff_combined[[1]], breaks=25, main = "Histogram of ALFF data for ROI 1", xlab = "Value")

#### Determine which distribution to use for the Data Generating Process
# Determine the best models for the ROIs of the data blocks
investigate_ROIs_for_data(data_reho_combined)
investigate_ROIs_for_data(data_lfcd_combined)
investigate_ROIs_for_data(data_alff_combined)

#### Simulating data
# Select the amount of samples to generate:
n_samples <- 1000
# Set the distribution types to base simulation on
dist_type_reho <- "weibull" 
dist_type_lfcd <- "gamma" 
dist_type_alff <- "weibull" 
# Copy the data so the original data is not altered
data_reho <- data_reho_combined
data_lfcd <- data_lfcd_combined
data_alff <- data_alff_combined
# Fit the data to the corresponding distribution type:
marginals_reho_combined <- fitted_marginals_for_data(data_reho, dist_type_reho)
marginals_lfcd_combined <- fitted_marginals_for_data(data_lfcd, dist_type_lfcd)
marginals_alff_combined <- fitted_marginals_for_data(data_alff, dist_type_alff)

# Use (vine) copulas to capture the common and individual structure for the data blocks
# First, create pseudo-observations of the data:
pseudo_reho <- pobs(data_reho)
pseudo_lfcd <- pobs(data_lfcd)
pseudo_alff <- pobs(data_alff)

# Feed the pseudo-observations to this vine copula function to find underlying structure of data blocks
vine_copula_model <- RVineStructureSelect(cbind(pseudo_reho, pseudo_lfcd, pseudo_alff),
                                          familyset = c(1), 
                                          selectioncrit = "logLik", indeptest = TRUE)
# Use vine copula structure to retain structural dependencies and create simulated data
simulated_pseudo_data <- RVineSim(n_samples, vine_copula_model)
# Transform simulated data back to original marginal values corresponding to original data
simulated_reho_combined <- transform_pseudo_to_original_marginals(simulated_pseudo_data[, 1:ncol(data_reho)],
                                                                marginals_reho_combined, dist_type_reho)
simulated_lfcd_combined <- transform_pseudo_to_original_marginals(simulated_pseudo_data[, (ncol(data_reho)+1):(ncol(data_reho) + ncol(data_lfcd))],
                                                                marginals_lfcd_combined, dist_type_lfcd)
simulated_alff_combined <- transform_pseudo_to_original_marginals(simulated_pseudo_data[, (ncol(data_reho) + ncol(data_lfcd) + 1):(ncol(data_reho) + ncol(data_lfcd) + ncol(data_alff))],
                                                                marginals_alff_combined, dist_type_alff)

# Unfortunately, ROI 83 consistently surpasses expected values during simulation for ALFF
# To fix, ROI 83 is reduced so its mean is more sensible
mean_ROI_83 <- mean(simulated_alff_combined[, 83])
mean_other_ROI <- mean(colMeans(simulated_alff_combined)[-83])
simulated_alff_combined[, 83] <- (simulated_alff_combined[, 83] / (mean_ROI_83 / mean_other_ROI))

#### Plot correlation heatmaps to check whether simulated data has underlying block structure
# Plot the correlation heatmap for the original combined data
pheatmap(cor(cbind(data_reho_combined, data_lfcd_combined, data_alff_combined)), 
         cluster_rows = FALSE, cluster_cols = FALSE,
         labels_row = "", labels_col = "")
# Plot the correlation heatmap for the simulated Combined data
pheatmap(cor(cbind(simulated_reho_combined, simulated_lfcd_combined, simulated_alff_combined)),
         cluster_rows = FALSE, cluster_cols = FALSE,
         labels_row = "", labels_col = "")

#### Concatenate the simulated data and save it as a csv file
# Directory to save the file in:
save_directory <- "F:/Studie/Thesis/Thesis R Project/Simulated_data"
combined_simulated_combined <- cbind(simulated_reho_combined, simulated_lfcd_combined, simulated_alff_combined)
# Include number of samples in the filename
filename <- file.path(save_directory, "simulated_data_combined.csv")
# Write the  data to the CSV file
write.csv(combined_simulated_combined, file = filename, row.names = FALSE)

