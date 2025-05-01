### Code For Bootstrapped Analysis of Robustness
# Author: Timo Waling (581706tw)

#### Functions used
# Function to add noise based on the width of data
apply_noise_to_data <- function(data, width, noise_level) {
  # Create a copy in order to preserve the original input data 
  data_copy <- data
  # Simulate random noise based on width and noise level
  noise <- matrix(rnorm(length(data_copy), mean = 0, sd = width * noise_level),
                  nrow = nrow(data_copy), ncol = ncol(data_copy))
  # Add the noise to the data
  result <- data_copy + noise
  return(result)
}

# Function to add outliers
apply_outlier_to_data <- function(data, fraction = 0.05, outlier_strength = 10) {
  # Create a copy in order to preserve the original input data and ensure data is a matrix
  data_copy <- as.matrix(data)
  # Calculate statistics from data
  n_rows <- nrow(data_copy)
  n_cols <- ncol(data_copy)
  sd_from_data <- mean(apply(data_copy, 2, sd))
  # Pick a fraction of rows to contaminate input data
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
# Function that fits the input data into a DISCO-SCA model
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

# Create noise levels based on values of real data
width_reho_combined = (max(data_reho_combined) - min(data_reho_combined))
width_lfcd_combined = (max(data_lfcd_combined) - min(data_lfcd_combined))
width_alff_combined = (max(data_alff_combined) - min(data_alff_combined))


#### Bootstrapping Reho data for noise:
# Reho 1% Noise
reho_noise_01_boot_1 <- apply_noise_to_data(simulated_reho_combined, width_reho_combined, 0.01)
reho_noise_01_boot_2 <- apply_noise_to_data(simulated_reho_combined, width_reho_combined, 0.01)
  
# Reho 5% Noise
reho_noise_05_boot_1 <- apply_noise_to_data(simulated_reho_combined, width_reho_combined, 0.05)
reho_noise_05_boot_2 <- apply_noise_to_data(simulated_reho_combined, width_reho_combined, 0.05)
  
# Reho 10% Noise
reho_noise_10_boot_1 <- apply_noise_to_data(simulated_reho_combined, width_reho_combined, 0.1)
reho_noise_10_boot_2 <- apply_noise_to_data(simulated_reho_combined, width_reho_combined, 0.1)

# Reho 25% Noise
reho_noise_25_boot_1 <- apply_noise_to_data(simulated_reho_combined, width_reho_combined, 0.25)
reho_noise_25_boot_2 <- apply_noise_to_data(simulated_reho_combined, width_reho_combined, 0.25)
  
# Reho 50% Noise
reho_noise_50_boot_1 <- apply_noise_to_data(simulated_reho_combined, width_reho_combined, 0.5)
reho_noise_50_boot_2 <- apply_noise_to_data(simulated_reho_combined, width_reho_combined, 0.5)
  
#### Bootstrapping LFCD data for noise:
# lfcd 1% Noise
lfcd_noise_01_boot_1 <- apply_noise_to_data(simulated_lfcd_combined, width_lfcd_combined, 0.01)
lfcd_noise_01_boot_2 <- apply_noise_to_data(simulated_lfcd_combined, width_lfcd_combined, 0.01)
  
# lfcd 5% Noise
lfcd_noise_05_boot_1 <- apply_noise_to_data(simulated_lfcd_combined, width_lfcd_combined, 0.05)
lfcd_noise_05_boot_2 <- apply_noise_to_data(simulated_lfcd_combined, width_lfcd_combined, 0.05)
  
# lfcd 10% Noise
lfcd_noise_10_boot_1 <- apply_noise_to_data(simulated_lfcd_combined, width_lfcd_combined, 0.1)
lfcd_noise_10_boot_2 <- apply_noise_to_data(simulated_lfcd_combined, width_lfcd_combined, 0.1)
  
# lfcd 25% Noise
lfcd_noise_25_boot_1 <- apply_noise_to_data(simulated_lfcd_combined, width_lfcd_combined, 0.25)
lfcd_noise_25_boot_2 <- apply_noise_to_data(simulated_lfcd_combined, width_lfcd_combined, 0.25)
  
# lfcd 50% Noise
lfcd_noise_50_boot_1 <- apply_noise_to_data(simulated_lfcd_combined, width_lfcd_combined, 0.5)
lfcd_noise_50_boot_2 <- apply_noise_to_data(simulated_lfcd_combined, width_lfcd_combined, 0.5)
    
### Bootstrapping alff data for noise:
# alff 1% Noise
alff_noise_01_boot_1 <- apply_noise_to_data(simulated_alff_combined, width_alff_combined, 0.01)
alff_noise_01_boot_2 <- apply_noise_to_data(simulated_alff_combined, width_alff_combined, 0.01)

# alff 5% Noise
alff_noise_05_boot_1 <- apply_noise_to_data(simulated_alff_combined, width_alff_combined, 0.05)
alff_noise_05_boot_2 <- apply_noise_to_data(simulated_alff_combined, width_alff_combined, 0.05)
  
# alff 10% Noise
alff_noise_10_boot_1 <- apply_noise_to_data(simulated_alff_combined, width_alff_combined, 0.1)
alff_noise_10_boot_2 <- apply_noise_to_data(simulated_alff_combined, width_alff_combined, 0.1)
  
# alff 25% Noise
alff_noise_25_boot_1 <- apply_noise_to_data(simulated_alff_combined, width_alff_combined, 0.25)
alff_noise_25_boot_2 <- apply_noise_to_data(simulated_alff_combined, width_alff_combined, 0.25)

# alff 50% Noise
alff_noise_50_boot_1 <- apply_noise_to_data(simulated_alff_combined, width_alff_combined, 0.5)
alff_noise_50_boot_2 <- apply_noise_to_data(simulated_alff_combined, width_alff_combined, 0.5)
  
### Bootstrapping Reho data for outlier:
# Reho 1% outlier
reho_outlier_01_boot_1 <- apply_outlier_to_data(simulated_reho_combined, 0.01, 7)
reho_outlier_01_boot_2 <- apply_outlier_to_data(simulated_reho_combined, 0.01, 7)
  
# Reho 5% outlier
reho_outlier_05_boot_1 <- apply_outlier_to_data(simulated_reho_combined, 0.05, 7)
reho_outlier_05_boot_2 <- apply_outlier_to_data(simulated_reho_combined, 0.05, 7)
  
# Reho 10% outlier
reho_outlier_10_boot_1 <- apply_outlier_to_data(simulated_reho_combined, 0.1, 7)
reho_outlier_10_boot_2 <- apply_outlier_to_data(simulated_reho_combined, 0.1, 7)

# Reho 25% outlier
reho_outlier_25_boot_1 <- apply_outlier_to_data(simulated_reho_combined, 0.25, 7)
reho_outlier_25_boot_2 <- apply_outlier_to_data(simulated_reho_combined, 0.25, 7)

# Reho 50% outlier
reho_outlier_50_boot_1 <- apply_outlier_to_data(simulated_reho_combined, 0.5, 7)
reho_outlier_50_boot_2 <- apply_outlier_to_data(simulated_reho_combined, 0.5, 7)

### Bootstrapping LFCD data for outlier:
# lfcd 1% outlier
lfcd_outlier_01_boot_1 <- apply_outlier_to_data(simulated_lfcd_combined, 0.01, 7)
lfcd_outlier_01_boot_2 <- apply_outlier_to_data(simulated_lfcd_combined, 0.01, 7)

# lfcd 5% outlier
lfcd_outlier_05_boot_1 <- apply_outlier_to_data(simulated_lfcd_combined, 0.05, 7)
lfcd_outlier_05_boot_2 <- apply_outlier_to_data(simulated_lfcd_combined, 0.05, 7)

# lfcd 10% outlier
lfcd_outlier_10_boot_1 <- apply_outlier_to_data(simulated_lfcd_combined, 0.1, 7)
lfcd_outlier_10_boot_2 <- apply_outlier_to_data(simulated_lfcd_combined, 0.1, 7)

# lfcd 25% outlier
lfcd_outlier_25_boot_1 <- apply_outlier_to_data(simulated_lfcd_combined, 0.25, 7)
lfcd_outlier_25_boot_2 <- apply_outlier_to_data(simulated_lfcd_combined, 0.25, 7)

# lfcd 50% outlier
lfcd_outlier_50_boot_1 <- apply_outlier_to_data(simulated_lfcd_combined, 0.5, 7)
lfcd_outlier_50_boot_2 <- apply_outlier_to_data(simulated_lfcd_combined, 0.5, 7)

### Bootstrapping alff data for outlier:
# alff 1% outlier
alff_outlier_01_boot_1 <- apply_outlier_to_data(simulated_alff_combined, 0.01, 7)
alff_outlier_01_boot_2 <- apply_outlier_to_data(simulated_alff_combined, 0.01, 7)

# alff 5% outlier
alff_outlier_05_boot_1 <- apply_outlier_to_data(simulated_alff_combined, 0.05, 7)
alff_outlier_05_boot_2 <- apply_outlier_to_data(simulated_alff_combined, 0.05, 7)

# alff 10% outlier
alff_outlier_10_boot_1 <- apply_outlier_to_data(simulated_alff_combined, 0.1, 7)
alff_outlier_10_boot_2 <- apply_outlier_to_data(simulated_alff_combined, 0.1, 7)

# alff 25% outlier
alff_outlier_25_boot_1 <- apply_outlier_to_data(simulated_alff_combined, 0.25, 7)
alff_outlier_25_boot_2 <- apply_outlier_to_data(simulated_alff_combined, 0.25, 7)

# alff 50% outlier
alff_outlier_50_boot_1 <- apply_outlier_to_data(simulated_alff_combined, 0.5, 7)
alff_outlier_50_boot_2 <- apply_outlier_to_data(simulated_alff_combined, 0.5, 7)


#### Create datablocks that will be used to make JIVE models
# For noise 1%
block_01_noise_boot_1 <- list(reho = t(reho_noise_01_boot_1), lfcd = t(lfcd_noise_01_boot_1), alff = t(alff_noise_01_boot_1))
block_01_noise_boot_2 <- list(reho = t(reho_noise_01_boot_1), lfcd = t(lfcd_noise_01_boot_2), alff = t(alff_noise_01_boot_2))

# For noise 5%
block_05_noise_boot_1 <- list(reho = t(reho_noise_05_boot_1), lfcd = t(lfcd_noise_05_boot_1), alff = t(alff_noise_05_boot_1))
block_05_noise_boot_2 <- list(reho = t(reho_noise_05_boot_1), lfcd = t(lfcd_noise_05_boot_2), alff = t(alff_noise_05_boot_2))

# For noise 10%
block_10_noise_boot_1 <- list(reho = t(reho_noise_10_boot_1), lfcd = t(lfcd_noise_10_boot_1), alff = t(alff_noise_10_boot_1))
block_10_noise_boot_2 <- list(reho = t(reho_noise_10_boot_1), lfcd = t(lfcd_noise_10_boot_2), alff = t(alff_noise_10_boot_2))

# For noise 25%
block_25_noise_boot_1 <- list(reho = t(reho_noise_25_boot_1), lfcd = t(lfcd_noise_25_boot_1), alff = t(alff_noise_25_boot_1))
block_25_noise_boot_2 <- list(reho = t(reho_noise_25_boot_1), lfcd = t(lfcd_noise_25_boot_2), alff = t(alff_noise_25_boot_2))

# For noise 50%
block_50_noise_boot_1 <- list(reho = t(reho_noise_50_boot_1), lfcd = t(lfcd_noise_50_boot_1), alff = t(alff_noise_50_boot_1))
block_50_noise_boot_2 <- list(reho = t(reho_noise_50_boot_1), lfcd = t(lfcd_noise_50_boot_2), alff = t(alff_noise_50_boot_2))


# For outlier 1%
block_01_outlier_boot_1 <- list(reho = t(reho_outlier_01_boot_1), lfcd = t(lfcd_outlier_01_boot_1), alff = t(alff_outlier_01_boot_1))
block_01_outlier_boot_2 <- list(reho = t(reho_outlier_01_boot_1), lfcd = t(lfcd_outlier_01_boot_2), alff = t(alff_outlier_01_boot_2))

# For outlier 5%
block_05_outlier_boot_1 <- list(reho = t(reho_outlier_05_boot_1), lfcd = t(lfcd_outlier_05_boot_1), alff = t(alff_outlier_05_boot_1))
block_05_outlier_boot_2 <- list(reho = t(reho_outlier_05_boot_1), lfcd = t(lfcd_outlier_05_boot_2), alff = t(alff_outlier_05_boot_2))

# For outlier 10%
block_10_outlier_boot_1 <- list(reho = t(reho_outlier_10_boot_1), lfcd = t(lfcd_outlier_10_boot_1), alff = t(alff_outlier_10_boot_1))
block_10_outlier_boot_2 <- list(reho = t(reho_outlier_10_boot_1), lfcd = t(lfcd_outlier_10_boot_2), alff = t(alff_outlier_10_boot_2))

# For outlier 25%
block_25_outlier_boot_1 <- list(reho = t(reho_outlier_25_boot_1), lfcd = t(lfcd_outlier_25_boot_1), alff = t(alff_outlier_25_boot_1))
block_25_outlier_boot_2 <- list(reho = t(reho_outlier_25_boot_1), lfcd = t(lfcd_outlier_25_boot_2), alff = t(alff_outlier_25_boot_2))

# For outlier 50%
block_50_outlier_boot_1 <- list(reho = t(reho_outlier_50_boot_1), lfcd = t(lfcd_outlier_50_boot_1), alff = t(alff_outlier_50_boot_1))
block_50_outlier_boot_2 <- list(reho = t(reho_outlier_50_boot_1), lfcd = t(lfcd_outlier_50_boot_2), alff = t(alff_outlier_50_boot_2))

######## DISCO-SCA
# Set number of components
n_comp = 3

# Create DISCO-SCA models on the noise data:
DISCO_noise_01_boot_1 <- create_DISCOSCA_model(reho_noise_01_boot_1, lfcd_noise_01_boot_1, alff_noise_01_boot_1, n_comp)
DISCO_noise_01_boot_2 <- create_DISCOSCA_model(reho_noise_01_boot_2, lfcd_noise_01_boot_2, alff_noise_01_boot_2, n_comp)
DISCO_noise_05_boot_1 <- create_DISCOSCA_model(reho_noise_05_boot_1, lfcd_noise_05_boot_1, alff_noise_05_boot_1, n_comp)
DISCO_noise_05_boot_2 <- create_DISCOSCA_model(reho_noise_05_boot_2, lfcd_noise_05_boot_2, alff_noise_05_boot_2, n_comp)
DISCO_noise_10_boot_1 <- create_DISCOSCA_model(reho_noise_10_boot_1, lfcd_noise_10_boot_1, alff_noise_10_boot_1, n_comp)
DISCO_noise_10_boot_2 <- create_DISCOSCA_model(reho_noise_10_boot_2, lfcd_noise_10_boot_2, alff_noise_10_boot_2, n_comp)
DISCO_noise_25_boot_1 <- create_DISCOSCA_model(reho_noise_25_boot_1, lfcd_noise_25_boot_1, alff_noise_25_boot_1, n_comp)
DISCO_noise_25_boot_2 <- create_DISCOSCA_model(reho_noise_25_boot_2, lfcd_noise_25_boot_2, alff_noise_25_boot_2, n_comp)
DISCO_noise_50_boot_1 <- create_DISCOSCA_model(reho_noise_50_boot_1, lfcd_noise_50_boot_1, alff_noise_50_boot_1, n_comp)
DISCO_noise_50_boot_2 <- create_DISCOSCA_model(reho_noise_50_boot_2, lfcd_noise_50_boot_2, alff_noise_50_boot_2, n_comp)

# Create DISCO-SCA models on the outlier data:
DISCO_outlier_01_boot_1 <- create_DISCOSCA_model(reho_outlier_01_boot_1, lfcd_outlier_01_boot_1, alff_outlier_01_boot_1, n_comp)
DISCO_outlier_01_boot_2 <- create_DISCOSCA_model(reho_outlier_01_boot_2, lfcd_outlier_01_boot_2, alff_outlier_01_boot_2, n_comp)
DISCO_outlier_05_boot_1 <- create_DISCOSCA_model(reho_outlier_05_boot_1, lfcd_outlier_05_boot_1, alff_outlier_05_boot_1, n_comp)
DISCO_outlier_05_boot_2 <- create_DISCOSCA_model(reho_outlier_05_boot_2, lfcd_outlier_05_boot_2, alff_outlier_05_boot_2, n_comp)
DISCO_outlier_10_boot_1 <- create_DISCOSCA_model(reho_outlier_10_boot_1, lfcd_outlier_10_boot_1, alff_outlier_10_boot_1, n_comp)
DISCO_outlier_10_boot_2 <- create_DISCOSCA_model(reho_outlier_10_boot_2, lfcd_outlier_10_boot_2, alff_outlier_10_boot_2, n_comp)
DISCO_outlier_25_boot_1 <- create_DISCOSCA_model(reho_outlier_25_boot_1, lfcd_outlier_25_boot_1, alff_outlier_25_boot_1, n_comp)
DISCO_outlier_25_boot_2 <- create_DISCOSCA_model(reho_outlier_25_boot_2, lfcd_outlier_25_boot_2, alff_outlier_25_boot_2, n_comp)
DISCO_outlier_50_boot_1 <- create_DISCOSCA_model(reho_outlier_50_boot_1, lfcd_outlier_50_boot_1, alff_outlier_50_boot_1, n_comp)
DISCO_outlier_50_boot_2 <- create_DISCOSCA_model(reho_outlier_50_boot_2, lfcd_outlier_50_boot_2, alff_outlier_50_boot_2, n_comp)

######## JIVE
# Create JIVE models on the noise data:
JIVE_noise_01_boot_1 <- jive(block_01_noise_boot_1, method = "given", rankJ = 3, rankA = c(2,2,2))
JIVE_noise_01_boot_2 <- jive(block_01_noise_boot_2, method = "given", rankJ = 3, rankA = c(2,2,2))
JIVE_noise_05_boot_1 <- jive(block_05_noise_boot_1, method = "given", rankJ = 3, rankA = c(2,2,2))
JIVE_noise_05_boot_2 <- jive(block_05_noise_boot_2, method = "given", rankJ = 3, rankA = c(2,2,2))
JIVE_noise_10_boot_1 <- jive(block_10_noise_boot_1, method = "given", rankJ = 3, rankA = c(2,2,2))
JIVE_noise_10_boot_2 <- jive(block_10_noise_boot_2, method = "given", rankJ = 3, rankA = c(2,2,2))
JIVE_noise_25_boot_1 <- jive(block_25_noise_boot_1, method = "given", rankJ = 3, rankA = c(2,2,2))
JIVE_noise_25_boot_2 <- jive(block_25_noise_boot_2, method = "given", rankJ = 3, rankA = c(2,2,2))
JIVE_noise_50_boot_1 <- jive(block_50_noise_boot_1, method = "given", rankJ = 3, rankA = c(2,2,2))
JIVE_noise_50_boot_2 <- jive(block_50_noise_boot_2, method = "given", rankJ = 3, rankA = c(2,2,2))

# Create JIVE models on the outlier data:
JIVE_outlier_01_boot_1 <- jive(block_01_outlier_boot_1, method = "given", rankJ = 3, rankA = c(2,2,2))
JIVE_outlier_01_boot_2 <- jive(block_01_outlier_boot_2, method = "given", rankJ = 3, rankA = c(2,2,2))
JIVE_outlier_05_boot_1 <- jive(block_05_outlier_boot_1, method = "given", rankJ = 3, rankA = c(2,2,2))
JIVE_outlier_05_boot_2 <- jive(block_05_outlier_boot_2, method = "given", rankJ = 3, rankA = c(2,2,2))
JIVE_outlier_10_boot_1 <- jive(block_10_outlier_boot_1, method = "given", rankJ = 3, rankA = c(2,2,2))
JIVE_outlier_10_boot_2 <- jive(block_10_outlier_boot_2, method = "given", rankJ = 3, rankA = c(2,2,2))
JIVE_outlier_25_boot_1 <- jive(block_25_outlier_boot_1, method = "given", rankJ = 3, rankA = c(2,2,2))
JIVE_outlier_25_boot_2 <- jive(block_25_outlier_boot_2, method = "given", rankJ = 3, rankA = c(2,2,2))
JIVE_outlier_50_boot_1 <- jive(block_50_outlier_boot_1, method = "given", rankJ = 3, rankA = c(2,2,2))
JIVE_outlier_50_boot_2 <- jive(block_50_outlier_boot_2, method = "given", rankJ = 3, rankA = c(2,2,2))