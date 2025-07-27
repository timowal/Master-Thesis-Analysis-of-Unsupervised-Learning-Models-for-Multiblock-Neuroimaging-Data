### Code For Generating Data For Simulation Study
# Author: Timo Waling (581706tw)


## Install packages and load them:
install.packages("ggplot2")
install.packages("stats")
install.packages("fitdistrplus")
install.packages("copula")
install.packages("dplyr")
install.packages("pheatmap")
install.packages("VineCopula")
install.packages("statmod")
install.packages("actuar")
install.packages("mixtools")
install.packages("multiblock")
install.packages("r.jive")
install.packages("parallel")
install.packages("rlist")
install.packages("clue")

library(ggplot2)
library(stats)
library(fitdistrplus)
library(copula)
library(dplyr)
library(pheatmap)
library(VineCopula)
library(statmod)
library(actuar)
library(mixtools)
library(multiblock)
library(r.jive)
library(parallel)
library(rlist)
library(clue)

## Use a seed for consistency in results, and set project directory
set.seed(2025)
project_directory <- "F:/Studie/Thesis/Thesis R Project/"

### Loading and cleaning data
## First, load in all the original data of the patients
# Loading in REHO data
data_reho_combined <- read.csv(paste0(project_directory, "data/reho_data/reho_atlas_combined.csv"),
                               header = FALSE)
data_reho_combined <- data_reho_combined[-1, ] #Remove first row of 0's


# Loading in LFCD data
data_lfcd_combined <- read.csv(paste0(project_directory, "data/lfcd_data/lfcd_atlas_combined.csv"),
                               header = FALSE)
data_lfcd_combined <- data_lfcd_combined[-1, ] #Remove first row of 0's


# Loading in ALFF data
data_alff_combined <- read.csv(paste0(project_directory, "data/alff_data/alff_atlas_combined.csv"),
                               header = FALSE)
data_alff_combined <- data_alff_combined[-1, ] #Remove first row of 0's



## We will work with distributions that require strictly positive values, so replace 0 with 1*10^(-18)
# For REHO
data_reho_combined[data_reho_combined == 0] <- 1e-18

# For LFCD
data_lfcd_combined[data_lfcd_combined == 0] <- 1e-18

# For ALFF
data_alff_combined[data_alff_combined == 0] <- 1e-18



### Test which distribution to use for the Data Generating Process in Simulation Study

## Create a function that takes data of one data block for one ROI, and fits the data to a
## Gamma, Exponential, Weibull, Log-Normal, Log-Logistic and Bimodal distribution in order to find
## Which distribution best fits the data type

fit_data_to_distributions <- function(data) {
  ## Gamma distribution:
  # Fit the data onto a Gamma distribution:
  fitted_gamma_model <- fitdist(data, "gamma", method = "mle")

  # Get log-likelihood, AIC and BIC:
  fitted_gamma_logl <- fitted_gamma_model$loglik
  fitted_gamma_aic <- fitted_gamma_model$aic
  fitted_gamma_bic <- fitted_gamma_model$bic
  
  # Get the average value of 10 Kolmogorov-Smirnov tests for Gamma distribution
  ks_stats <- numeric(10)
  for (i in 1:10) {
    simulated_data <- rgamma(length(data), shape = fitted_gamma_model$estimate["shape"],
                             rate = fitted_gamma_model$estimate["rate"])
    
    ks_stats[i] <- ks.test(data, simulated_data, exact = FALSE)$statistic
  }
  fitted_gamma_ks <- mean(ks_stats)
  
  
  ## Exponential distribution
  fitted_exponential_model <- fitdist(data, "exp", method = "mle")
  
  # Get log-likelihood, AIC and BIC:
  fitted_exponential_logl <- fitted_exponential_model$loglik
  fitted_exponential_aic <- fitted_exponential_model$aic
  fitted_exponential_bic <- fitted_exponential_model$bic
  
  # Get the average value of 10 Kolmogorov-Smirnov tests for Exponential distribution
  ks_stats <- numeric(10)
  for (i in 1:10) {
    simulated_data <- rexp(length(data), rate = fitted_exponential_model$estimate["rate"])
    
    ks_stats[i] <- ks.test(data, simulated_data, exact = FALSE)$statistic
  }
  fitted_exponential_ks <- mean(ks_stats)
  
  
  ## Weibull distribution
  # Fit the data onto a Weibull distribution:
  fitted_weibull_model <- fitdist(data, "weibull", method = "mle")
  
  # Get log-likelihood, AIC and BIC:
  fitted_weibull_logl <- fitted_weibull_model$loglik
  fitted_weibull_aic <- fitted_weibull_model$aic
  fitted_weibull_bic <- fitted_weibull_model$bic
  
  # Get the average value of 10 Kolmogorov-Smirnov tests for Weibull distribution
  ks_stats <- numeric(10)
  for (i in 1:10) {
    simulated_data <- rweibull(length(data), shape = fitted_weibull_model$estimate["shape"],
                               scale = fitted_weibull_model$estimate["scale"])
    
    ks_stats[i] <- ks.test(data, simulated_data, exact = FALSE)$statistic
  }
  fitted_weibull_ks <- mean(ks_stats)
  
  
  
  ## Log-Normal distribution
  # Fit the data onto a Log-Normal distribution:
  fitted_lognormal_model <- fitdist(data, "lnorm", method = "mle")
  
  # Get log-likelihood, AIC and BIC:
  fitted_lognormal_logl <- fitted_lognormal_model$loglik
  fitted_lognormal_aic <- fitted_lognormal_model$aic
  fitted_lognormal_bic <- fitted_lognormal_model$bic
  
  # Get the average value of 10 Kolmogorov-Smirnov tests for Log-Normal distribution
  ks_stats <- numeric(10)
  for (i in 1:10) {
    simulated_data <- rlnorm(length(data), meanlog = fitted_lognormal_model$estimate["meanlog"],
                             sdlog = fitted_lognormal_model$estimate["sdlog"])
    
    ks_stats[i] <- ks.test(data, simulated_data, exact = FALSE)$statistic
  }
  fitted_lognormal_ks <- mean(ks_stats)
  
  
  
  ## Log-Logistic distribution
  # Fit the data onto a Log-Logistic distribution:
  fitted_loglogistic_model <- fitdist(data, "llogis", method = "mle")
  
  # Get log-likelihood, AIC and BIC:
  fitted_loglogistic_logl <- fitted_loglogistic_model$loglik
  fitted_loglogistic_aic <- fitted_loglogistic_model$aic
  fitted_loglogistic_bic <- fitted_loglogistic_model$bic
  
  # Get the average value of 10 Kolmogorov-Smirnov tests for Log-Logistic distribution
  ks_stats <- numeric(10)
  for (i in 1:10) {
    simulated_data <- rllogis(length(data), shape = fitted_loglogistic_model$estimate["shape"],
                              scale = fitted_loglogistic_model$estimate["scale"])
    
    ks_stats[i] <- ks.test(data, simulated_data, exact = FALSE)$statistic
  }
  fitted_loglogistic_ks <- mean(ks_stats)
  
  
  
  ## Bimodal Mixture Distribution
  # First subtract the mean of the data from the data, to make the data centered around 0
  data_mean <- mean(data)
  centered_data <- data - data_mean
  
  # Fit the data onto a Bimodal distribution:
  fitted_bimodal_model <- normalmixEM(centered_data, k=2, maxit=5000, epsilon=1e-8)
  
  # Get log-likelihood, AIC and BIC:
  fitted_bimodal_logl <- fitted_bimodal_model$loglik
  fitted_bimodal_aic <- -2 * fitted_bimodal_logl + 2 * 6
  fitted_bimodal_bic <- -2 * fitted_bimodal_logl + log(length(data)) * 6
  
  # Get the average value of 10 Kolmogorov-Smirnov tests for Bimodal mixture model
  ks_stats <- numeric(10)
  for (i in 1:10) {
    simulated_data <- numeric(length(data))
    
    for (j in 1:length(data)) {
      # Choose component 1 or 2 based on mixing probabilities
      comp <- sample(1:2, size = 1, prob = fitted_bimodal_model$lambda)
      
      # Draw from corresponding normal distribution and shift mean back to that of original data
      drawn_value <- rnorm(1, mean = fitted_bimodal_model$mu[comp], 
                           sd = fitted_bimodal_model$sigma[comp])
      
      drawn_value <- drawn_value + data_mean
      
      # Truncate negative values
      if (drawn_value < 0) {
        simulated_data[j] <- 0
      } else {
        simulated_data[j] <- drawn_value
      }
    }
    
    # Run KS test and store statistic
    ks_stats[i] <- ks.test(data, simulated_data, exact = FALSE)$statistic
  }
  fitted_bimodal_ks <- mean(ks_stats)
  
  
  
  # Concatenate the log-likelihood, AIC and BIC values
  values_logl <- c(fitted_gamma_logl, fitted_exponential_logl, fitted_weibull_logl,
                   fitted_lognormal_logl, fitted_loglogistic_logl, fitted_bimodal_logl)
  
  values_aic <- c(fitted_gamma_aic, fitted_exponential_aic, fitted_weibull_aic,
                  fitted_lognormal_aic, fitted_loglogistic_aic, fitted_bimodal_aic)
  
  values_bic <- c(fitted_gamma_bic, fitted_exponential_bic, fitted_weibull_bic,
                  fitted_lognormal_bic, fitted_loglogistic_bic, fitted_bimodal_bic)
  
  values_ks <- c(fitted_gamma_ks, fitted_exponential_ks, fitted_weibull_ks,
                  fitted_lognormal_ks, fitted_loglogistic_ks, fitted_bimodal_ks)
  
  # Obtain the index corresponding to the distribution with the best value
  best_logl_index <- which.max(values_logl) # For log-likelihood, higher values are preferred
  best_aic_index <- which.min(values_aic) # For AIC, lower values are preferred
  best_bic_index <- which.min(values_bic) # For BIC, lower values are preferred
  best_ks_index <- which.min(values_ks) # For KS, lower values are preferred
  
  # Return the indices of the best fitting distributions
  result <- list(best_logl_index = best_logl_index, best_aic_index = best_aic_index,
                 best_bic_index = best_bic_index, best_ks_index = best_ks_index)
  
  return(result)
}

## Create a function that loops through the ROIs to evaluate which distribution best describes it
investigate_ROIs_for_data <- function(data) {
  ## Distribution would best describe the data generating process behind the data 
  counter_logl_best <- c("Gamma" = 0, "Exponential" = 0, "Weibull" = 0,"Log-Normal" = 0,
                         "Log-Logistic" = 0, "Bimodal Mixture" = 0)
  
  counter_aic_best <- c("Gamma" = 0, "Exponential" = 0, "Weibull" = 0, "Log-Normal" = 0,
                        "Log-Logistic" = 0, "Bimodal Mixture" = 0)
  
  counter_bic_best <- c("Gamma" = 0, "Exponential" = 0, "Weibull" = 0, "Log-Normal" = 0,
                        "Log-Logistic" = 0, "Bimodal Mixture" = 0)
  
  counter_ks_best <- c("Gamma" = 0, "Exponential" = 0, "Weibull" = 0, "Log-Normal" = 0,
                       "Log-Logistic" = 0, "Bimodal Mixture" = 0)
  
  
  # Loop through the 111 ROIs and fit the data to the distributions
  for (i in 1:111) {
    investigate_ROI <- data[, i] #Get the data for the i-th ROI
    
    # Get the log-likelihood, AIC and BIC values for fitting this ROI data onto the distributions
    fitted_distribution_scores <- fit_data_to_distributions(investigate_ROI)
    
    # Check which distribution had the highest log-likelihood value
    best_logl_dist <- c("Gamma", "Exponential", "Weibull", "Log-Normal", "Log-Logistic",
                        "Bimodal Mixture")[fitted_distribution_scores$best_logl_index]
    
    # Check which distribution had the lowest AIC value
    best_aic_dist <- c("Gamma", "Exponential", "Weibull", "Log-Normal", "Log-Logistic",
                       "Bimodal Mixture")[fitted_distribution_scores$best_aic_index]
    
    # Check which distribution had the lowest BIC value
    best_bic_dist <- c("Gamma", "Exponential", "Weibull", "Log-Normal", "Log-Logistic",
                       "Bimodal Mixture")[fitted_distribution_scores$best_bic_index]
    
    # Check which distribution had the lowest KS test value
    best_ks_dist <- c("Gamma", "Exponential", "Weibull", "Log-Normal", "Log-Logistic",
                      "Bimodal Mixture")[fitted_distribution_scores$best_ks_index]
    
    # Update the counters for how many times a certain distribution provided the best fit for the ROI
    counter_logl_best[best_logl_dist] <- counter_logl_best[best_logl_dist] + 1
    counter_aic_best[best_aic_dist] <- counter_aic_best[best_aic_dist] + 1
    counter_bic_best[best_bic_dist] <- counter_bic_best[best_bic_dist] + 1
    counter_ks_best[best_ks_dist] <- counter_ks_best[best_ks_dist] + 1
  }
  
  ## After going through the 111 ROIs, we print the result for which distribution was best how many times
  table_of_results <- data.frame(Distribution = c("Gamma", "Exponential", "Weibull", "Log-Normal",
                                                  "Log-Logistic", "Bimodal Mixture"),
                                 amount_logl = counter_logl_best, amount_aic = counter_aic_best,
                                 amount_bic = counter_bic_best, amount_ks = counter_ks_best)
  
  print(table_of_results)
}

# Determine the best models for the ROIs of the data blocks
investigate_ROIs_for_data(data_reho_combined)
investigate_ROIs_for_data(data_lfcd_combined)
investigate_ROIs_for_data(data_alff_combined)



##### FUNCTIONS FOR SIMULATING DATA

## Functions used for simulating data
fitted_marginals_for_data <- function(data, distribution_type = "gamma") {
  # Create the list to save the marginal distributions of the data
  result <- list()
  
  # Iterate through the ROIs for the data
  for (i in 1:ncol(data)) {
    if (distribution_type == "bimodal") {
      # For Bimodal, first ensure the data is centered around 0
      data_mean <- mean(data[, i])
      centered_data <- data[, i] - data_mean
      
      # Fit the ROI data to the Bimodal distribution
      ROI_fit <- normalmixEM(centered_data, k = 2, maxit = 5000, epsilon= 1e-8)
      
      result[[i]] <- list(fit = ROI_fit, mean = data_mean)
    } else {
      # Fit the ROI data to the given distribution
      ROI_fit <- fitdist(data[, i], distribution_type, method = "mle")
      
      result[[i]] <- ROI_fit
    }
  }
  
  # Return the list of marginals
  return(result)
}


## Function that turns the simulated uniform pseudo data back to the original marginals
transform_pseudo_to_original_marginals <- function(simulated_data, marginals, distribution_type = "gamma") {
  # Create an empty matrix we will return the transformed data in
  result <- matrix(NA, nrow(simulated_data), length(marginals))
  
  # Loop through the ROIs we will transform
  for (i in 1:ncol(simulated_data)) {
    # For Gamma distribution:
    if (distribution_type == "gamma") {
      param_shape <- marginals[[i]]$estimate[1] # Shape parameter for the Gamma function
      param_rate <- marginals[[i]]$estimate[2] # Rate parameter for the Gamma function
      
      # Transformation using the Gamma marginals
      result[, i] <- qgamma(simulated_data[, i], shape = param_shape, rate = param_rate)
    }
    
    # For Exponential distribution:
    if (distribution_type == "exp") {
      param_rate <- marginals[[i]]$estimate[1] # Rate parameter for Exponential distribution

      # Transformation using the Exponential marginals
      result[, i] <- qexp(simulated_data[, i], rate = param_rate)
    }
    
    # For Weibull distribution:
    if (distribution_type == "weibull") {
      param_shape <- marginals[[i]]$estimate[1] # Shape parameter for the Weibull function
      param_scale <- marginals[[i]]$estimate[2] # Scale parameter for the Weibull function
      
      # Transformation using the Weibull marginals
      result[, i] <- qweibull(simulated_data[, i], shape = param_shape, scale = param_scale)
    }
    
    # For Log-normal distribution:
    if (distribution_type == "lnorm") {
      param_meanlog <- marginals[[i]]$estimate[1]
      param_sdlog <- marginals[[i]]$estimate[2]
      
      # Transformation using the Log-Normal marginals
      result[, i] <- qlnorm(simulated_data[, i], meanlog = param_meanlog, sdlog = param_sdlog)
    }
    
    # For Log-logistic distribution:
    if (distribution_type == "llogis") {
      param_shape <- marginals[[i]]$estimate[1]  # Shape parameter for the Log-Logistic function
      param_scale <- marginals[[i]]$estimate[2]  # Scale parameter for the Log-Logistic function
      
      # Transformation using the Log-Logistic marginals
      result[, i] <- qllogis(simulated_data[, i], shape = param_shape, scale = param_scale)
    }
    
    # For Bimodal Mixture:
    if (distribution_type == "bimodal") {
      ROI_fit <- marginals[[i]]$fit
      ROI_data_mean <- marginals[[i]]$mean
      
      for (j in 1:nrow(simulated_data)) {
        simulated_value <- simulated_data[j, i]
        
        # This function serves as the CDF of our Bimodal distribution
        bimodal_mixture_cdf <- function(x) {
            ROI_fit$lambda[1] * pnorm(x, mean = ROI_fit$mu[1], sd = ROI_fit$sigma[1]) +
            ROI_fit$lambda[2] * pnorm(x, mean = ROI_fit$mu[2], sd = ROI_fit$sigma[2])
        }
        
        # Set an interval in which we try to find quantiles for this ROI
        min_interval <- min(ROI_fit$mu) - 10 * max(ROI_fit$sigma)
        max_interval <- max(ROI_fit$mu) + 10 * max(ROI_fit$sigma)
        
        # Try to find quantiles for the Bimodal distribution appropritate for this ROI data
        drawn_value <- tryCatch({
          uniroot(
            f = function(x) bimodal_mixture_cdf(x) - simulated_value,
            interval = c(min_interval, max_interval),
            tol = 1e-8
          )$root
        }, error = function(e) {
          # If necessary, retry with an even wider interval
          tryCatch({
            uniroot(
              f = function(x) bimodal_mixture_cdf(x) - simulated_value,
              interval = c(min_interval - 10, max_interval + 10),
              tol = 1e-8
            )$root
          })
        })
        
        # Add mean back in
        drawn_value <- drawn_value + ROI_data_mean
        
        # Truncate if negative
        if (drawn_value < 0) {
          result[j, i] <- 0
        } else {
          result[j, i] <- drawn_value
        }
      }
    }
  }
  
  # Return the transformed data
  return(result)
}

# Function to check correlation for original and simulated data for a certain ROI
check_correlation_ROI <- function(original_data, simulated_data) {
  correlation_original <- cor(original_data)
  correlation_simulated <- cor(simulated_data)

  print("Correlation for this ROI for Original Data:")
  print(correlation_original)
  
  print("Correlation for this ROI for Simulated Data:")
  print(correlation_simulated)
}

# Function to compute average Tucker congruence coefficient between two DISCO-SCA models
compute_tucker_coef_disco <- function(model_1, model_2) {
  # Extract the score matrices from both DISCO-SCA models
  scores_model_1 <- model_1$scores
  scores_model_2 <- model_2$scores
  
  # Initialize a matrix to store the Tucker coefficients
  n_comp <- ncol(scores_model_1)
  tucker_scores_matrix <- matrix(0, n_comp, n_comp)
  
  # Compute Tucker congruence coefficient for all possible component pairs
  for (i in 1:n_comp) {
    for (j in 1:n_comp) {
      # Extract the scores of the current components
      component_model_1 <- scores_model_1[, i]
      component_model_2 <- scores_model_2[, j]
      
      # Calculate the dot product of the current components
      dot_product <- sum(component_model_1 * component_model_2)
      
      # Calculate the norms of the current components
      norm_1 <- sqrt(sum(component_model_1^2))
      norm_2 <- sqrt(sum(component_model_2^2))
      
      # Compute the Tucker Congruence Coefficient for the current component pair
      tucker_scores_matrix[i, j] <- abs(dot_product) / (norm_1 * norm_2)
    }
  }
  
  # Find the best component pairings
  best_pairings <- solve_LSAP(tucker_scores_matrix, maximum = TRUE)
  best_component_scores <- tucker_scores_matrix[cbind(1:n_comp, best_pairings)]
  
  # Return the mean of the best pairings
  return(mean(best_component_scores))
}

# Function to compute Pearson coefficients between two DISCO-SCA models
compute_pearson_score_disco <- function(model_1, model_2) {
  # Extract loading matrices from both DISCO-SCA models
  loadings_model_1 <- model_1$loadings
  loadings_model_2 <- model_2$loadings
  
  # Initialize a matrix to store the Pearson coefficients
  n_comp <- ncol(loadings_model_1)
  pearson_scores_matrix <- matrix(0, n_comp, n_comp)
  
  # Compute Pearson correlation coefficient for all possible component pairs
  for (i in 1:n_comp) {
    for (j in 1:n_comp) {
      comp_1 <- loadings_model_1[, i]
      comp_2 <- loadings_model_2[, j]
      pearson_scores_matrix[i, j] <- abs(cor(comp_1, comp_2))
    }
  }
  
  # Find the best component pairings
  best_pairings <- solve_LSAP(pearson_scores_matrix, maximum = TRUE)
  best_component_scores <- pearson_scores_matrix[cbind(1:n_comp, best_pairings)]
  
  # Return the mean of the best component pairings
  return(mean(best_component_scores))
}

# Function to compute average Tucker congruence coefficient between two JIVE models
compute_tucker_coef_jive <- function(model_1, model_2) {
  # Extract joint components from both JIVE models
  joint_model_1 <- model_1$joint
  joint_model_2 <- model_2$joint
  
  # Ensure that both are matrices
  joint_model_1 <- lapply(joint_model_1, as.matrix)
  joint_model_2 <- lapply(joint_model_2, as.matrix)
  
  
  # Initialize matrix to store Tucker congruence between each pair of components
  n_comp <- length(joint_model_1)
  tucker_scores_matrix <- matrix(0, nrow = n_comp, ncol = n_comp)
  
  # Compute Tucker congruence coefficient for all possible component pairs
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


# Function to compute average Pearson score of components between two JIVE models
compute_pearson_score_jive <- function(model_1, model_2) {
  # Extract joint components from both JIVE models
  joint_model_1 <- model_1$joint
  joint_model_2 <- model_2$joint
  
  # Ensure that both are matrices
  joint_model_1 <- lapply(joint_model_1, as.matrix)
  joint_model_2 <- lapply(joint_model_2, as.matrix)
  
  # Initialize a matrix to store the Pearson correlations
  n_comp <- length(joint_model_1)
  pearson_scores_matrix <- matrix(0, n_comp, n_comp)
  
  # Compute Pearson correlation coefficient for all possible component pairs
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



##### SIMULATING DATA

# Select the amount of samples to generate:
set.seed(2025)
n_samples <- 884
n_simulation_datasets <- 10


# Set the distribution types to base simulation on
dist_type_reho <- "bimodal" 
dist_type_lfcd <- "llogis" 
dist_type_alff <- "bimodal" 


# Copy the data so the original data is not altered
data_reho <- data_reho_combined
data_lfcd <- data_lfcd_combined
data_alff <- data_alff_combined


# Create DISCO-SCA and JIVE models on the real data to later use for similarity comparison
data_blocks_real_jive <- list(reho = t(data_reho), lfcd = t(data_lfcd), alff = t(data_alff))

data_list_real_disco <- list(scale(as.matrix(data_reho)), scale(as.matrix(data_lfcd)),
                             scale(as.matrix(data_alff)))

disco_sca_model_real <- disco(data_list_real_disco, ncomp = 3)
jive_model_real <- jive(data_blocks_real_jive, method = "given", rankJ = 2, rankA = c(1, 1, 1))


# Fit the data to the corresponding distribution type:
marginals_reho_combined <- fitted_marginals_for_data(data_reho, dist_type_reho)
marginals_lfcd_combined <- fitted_marginals_for_data(data_lfcd, dist_type_lfcd)
marginals_alff_combined <- fitted_marginals_for_data(data_alff, dist_type_alff)


# Create pseudo-observations of our data:
pseudo_reho <- pobs(data_reho)
pseudo_lfcd <- pobs(data_lfcd)
pseudo_alff <- pobs(data_alff)

# Use the pseudo-observations to find underlying structure of ROIs of the data blocks
vine_copula_model <- RVineStructureSelect(cbind(pseudo_reho, pseudo_lfcd, pseudo_alff),
                                          familyset = c(1, 2), 
                                          selectioncrit = "logLik", indeptest = TRUE)


# Create multiple realizations of the simulation data sets
for (i in 1:n_simulation_datasets) {
  print(paste0("Starting Process to obtain dataset # ", i))
  
  # Simulate 5 datasets of which the best one will be selected
  simulated_datasets <- lapply(1:5, function(j) {
    # Use the vine copulas to obtain simulated pseudo observations
    simulated_pseudo_data <- RVineSim(n_samples, vine_copula_model)
    
    # Use the simulated pseudo observations to obtain simulated REHO data
    simulated_reho_combined <- transform_pseudo_to_original_marginals(
      simulated_pseudo_data[, 1:ncol(data_reho)], marginals_reho_combined, dist_type_reho)
    
    # Use the simulated pseudo observations to obtain simulated LFCD data
    simulated_lfcd_combined <- transform_pseudo_to_original_marginals(
      simulated_pseudo_data[, (ncol(data_reho)+1):(ncol(data_reho) + ncol(data_lfcd))],
      marginals_lfcd_combined, dist_type_lfcd)
    
    # Use the simulated pseudo observations to obtain simulated ALFF data
    simulated_alff_combined <- transform_pseudo_to_original_marginals(
      simulated_pseudo_data[, (ncol(data_reho) + ncol(data_lfcd) + 1):(ncol(data_reho) + ncol(data_lfcd) + ncol(data_alff))],
      marginals_alff_combined, dist_type_alff)
    
    combined_simulated_combined <- cbind(simulated_reho_combined, simulated_lfcd_combined, simulated_alff_combined)
    
    list(
      reho = simulated_reho_combined,
      lfcd = simulated_lfcd_combined,
      alff = simulated_alff_combined,
      combined = combined_simulated_combined
    )
  })
  
  print(paste0("Done creating 5 candidate simulated datasets for # ", i))
  
  # Settings for the parallel package to make us run 5 models at once
  n_cores <- min(5, detectCores() - 1)
  clusters_simulation <- makeCluster(n_cores)
  clusterEvalQ(clusters_simulation, {
    library(multiblock)
    library(r.jive)
    library(scales)
    library(clue)
  })

  
  # Fit the DISCO-SCA models on the 5 datasets simultaneously
  disco_sca_models <- parLapply(clusters_simulation, simulated_datasets, function(dataset) {
    data_list_simulated_disco <- list(
      scale(as.matrix(dataset$reho)),
      scale(as.matrix(dataset$lfcd)),
      scale(as.matrix(dataset$alff))
    )
    disco(data_list_simulated_disco, ncomp = 3)
  })
  
  print(paste0("Done creating DISCO-SCA models for iteration # ", i))
  
  # Fit the JIVE models on the 5 datasets simultaneously
  jive_models <- parLapply(clusters_simulation, simulated_datasets, function(dataset) {
    data_blocks_simulated_jive <- list(
      reho = t(dataset$reho),
      lfcd = t(dataset$lfcd),
      alff = t(dataset$alff)
    )
    jive(data_blocks_simulated_jive, method = "given", rankJ = 2, rankA = c(1, 1, 1))
  })
  
  print(paste0("Done creating JIVE models for iteration # ", i))
  
  # End the clusters of the parallel process
  stopCluster(clusters_simulation)
  
  
  # Set variables which will be used to select best simulated set
  best_score <- 0
  best_simulated_combined <- NULL
  
  # Select the best simulated dataset
  for (j in 1:5) {
    # Compute the Tucker's congruence coefficients for the DISCO-SCA and JIVE models on this realisation
    tucker_disco_sca <- compute_tucker_coef_disco(disco_sca_model_real, disco_sca_models[[j]])
    tucker_jive <- compute_tucker_coef_jive(jive_model_real, jive_models[[j]])
    
    # Compute the Pearson correlation coefficients for the DISCO-SCA and JIVE models on this realisation
    pearson_disco_sca <- compute_pearson_score_disco(disco_sca_model_real, disco_sca_models[[j]])
    pearson_jive <- compute_pearson_score_jive(jive_model_real, jive_models[[j]])
    
    # Compute the averages of the coefficient scores
    average_tucker <- mean(c(tucker_disco_sca, tucker_jive))
    average_pearson <- mean(c(pearson_disco_sca, pearson_jive))
    average_score <- mean(c(average_tucker, average_pearson))
    
    # Check if this is the best dataset so far
    if (average_score > best_score) {
      best_score <- average_score
      best_simulated_combined <- simulated_datasets[[j]]$combined
    }
  }
  
  # Directory to save the file in:
  save_directory <- paste0(project_directory, "Simulated_data/Uncontaminated")
  filename <- file.path(save_directory,  paste0("simulated_data_", i, ".csv"))
  
  # Write the  data to the CSV file
  write.csv(best_simulated_combined, file = filename, row.names = FALSE)
  
  print(paste0("Selected dataset ", i, " with score ", best_score, " at ", Sys.time()))
}
