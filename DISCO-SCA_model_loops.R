### Code For Creating all DISCO-SCA models
# Author: Timo Waling (581706tw)


## Install and load packages
install.packages("multiblock")
install.packages("neuroim")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("parallel")

library(multiblock)
library(neuroim)
library(dplyr)
library(ggplot2)
library(reshape2)
library(parallel)


## Use a seed for consistency in results, and set project directory
set.seed(2025)
project_directory <- "F:/Studie/Thesis/Thesis R Project/"

## Load in the original data
# Read the REHO data
reho_data_combined <- read.csv(file.path(project_directory, "data/reho_data/reho_atlas_combined.csv"), header = FALSE)
reho_data_control <- read.csv(file.path(project_directory, "data/reho_data/reho_atlas_Control_group.csv"), header = FALSE)
reho_data_asd <- read.csv(file.path(project_directory, "data/reho_data/reho_atlas_ASD_group.csv"), header = FALSE)
reho_data_combined <- reho_data_combined[-1, ]
reho_data_control <- reho_data_control[-1, ]
reho_data_asd <- reho_data_asd[-1, ]

# Read the LFCD data
lfcd_data_control <- read.csv(file.path(project_directory, "data/lfcd_data/lfcd_atlas_Control_group.csv"), header = FALSE)
lfcd_data_asd <- read.csv(file.path(project_directory, "data/lfcd_data/lfcd_atlas_ASD_group.csv"), header = FALSE)
lfcd_data_control <- lfcd_data_control[-1, ]
lfcd_data_asd <- lfcd_data_asd[-1, ]
lfcd_data_combined <- read.csv(file.path(project_directory, "data/lfcd_data/lfcd_atlas_combined.csv"), header = FALSE)
lfcd_data_combined <- lfcd_data_combined[-1, ]

# Read the ALFF data
alff_data_control <- read.csv(file.path(project_directory, "data/alff_data/alff_atlas_Control_group.csv"), header = FALSE)
alff_data_asd <- read.csv(file.path(project_directory, "data/alff_data/alff_atlas_ASD_group.csv"), header = FALSE)
alff_data_control <- alff_data_control[-1, ]
alff_data_asd <- alff_data_asd[-1, ]
alff_data_combined <- read.csv(file.path(project_directory, "data/alff_data/alff_atlas_combined.csv"), header = FALSE)
alff_data_combined <- alff_data_combined[-1, ]



##### Obtain DISCO-SCA models for the real data #####
data_list_real_control <- list(scale(as.matrix(reho_data_control)), scale(as.matrix(lfcd_data_control)),
                                  scale(as.matrix(alff_data_control)))

data_list_real_asd <- list(scale(as.matrix(reho_data_asd)), scale(as.matrix(lfcd_data_asd)),
                                  scale(as.matrix(alff_data_asd)))

data_list_real_combined <- list(scale(as.matrix(reho_data_combined)), scale(as.matrix(lfcd_data_combined)),
                                  scale(as.matrix(alff_data_combined)))


# Create DISCO-SCA models for the data
disco_sca_real_control <- disco(data_list_real_control, ncomp = 4)
disco_sca_real_asd <- disco(data_list_real_asd, ncomp = 4)
disco_sca_real_combined <- disco(data_list_real_combined, ncomp = 4)


# Save the DISCO-SCA models for later use
save_directory <- file.path(project_directory, "Results/Real_Data")

saveRDS(disco_sca_real_control, file = file.path(save_directory, "disco_sca_real_control.rds"))
saveRDS(disco_sca_real_asd, file = file.path(save_directory, "disco_sca_real_asd.rds"))
saveRDS(disco_sca_real_combined, file = file.path(save_directory, "disco_sca_real_combined.rds"))


##### Function for processing DISCO-SCA for given data #####
process_disco_sca_function <- function(i, data_directory, data_type, save_directory, save_type) {
  cat("Processing dataset", i, "from", data_type, "...\n")
  
  # Read the data and ensure it is numeric
  filename <- file.path(data_directory, paste0(data_type, i, ".csv"))
  full_data <- read.csv(filename, header = FALSE)[-1, ]
  full_data[] <- lapply(full_data, function(x) as.numeric(as.character(x)))
  
  data_reho <- full_data[, 1:111]
  data_lfcd <- full_data[, 112:222]
  data_alff <- full_data[, 223:333]
  
  data_list_full <- list(
    reho = scale(as.matrix(data_reho)),
    lfcd = scale(as.matrix(data_lfcd)),
    alff = scale(as.matrix(data_alff))
  )
  
  # Obtain the DISCO-SCA model for the data and save it
  disco_sca_model <- disco(data_list_full, ncomp = 4)
  saveRDS(disco_sca_model, file = file.path(save_directory, paste0(save_type, i, ".rds")))
  
  cat("Finished dataset", i, "\n")
}


##### Obtain DISCO-SCA Models for the Simulated Data without contamination #####
set.seed(2025)
n_simulation_datasets <- 10

# Settings for reading the uncontaminated data and saving the contaminated data
data_directory <- file.path(project_directory, "Simulated_data/Uncontaminated")
data_type <- "simulated_data_"
save_directory <- file.path(project_directory, "Results/Simulated_Data")
save_type <- "disco_sca_simulated_"


# Here we use clusters to simultaneously compute multiple models
n_cores <- min(10, detectCores() - 1)
clusters <- makeCluster(n_cores)

clusterExport(clusters, c("data_directory", "data_type", "save_directory", "save_type", "process_disco_sca_function"))
clusterEvalQ(clusters, { library(multiblock) })

parLapply(clusters, 1:n_simulation_datasets, function(i) {
  process_disco_sca_function(i, data_directory, data_type, save_directory, save_type)
})

stopCluster(clusters)



##### Obtain DISCO-SCA Models for the Simulated Data Noise Contamination #####
set.seed(2025)
n_simulation_datasets <- 10
noise_levels <- c("01", "05", "10", "25", "50")

for (noise_level in noise_levels) {
  cat("Processing noise level:", noise_level, "\n")
  
  # Here we use clusters to simultaneously compute multiple models
  n_cores <- min(10, detectCores() - 1)
  clusters <- makeCluster(n_cores)
  
  # Settings for reading the uncontaminated data and saving the contaminated data
  data_directory <- file.path(project_directory, "Simulated_data/Noise", paste0("Level_", noise_level))
  data_type <- paste0("noise_data_", noise_level, "_simulation_")
  save_directory <- file.path(project_directory, "Results/Noise_Results", paste0("Level_", noise_level))
  save_type <- paste0("disco_sca_noise_", noise_level, "_simulation_")
  
  clusterExport(clusters, c("data_directory", "data_type", "save_directory", "save_type", "process_disco_sca_function"))
  clusterEvalQ(clusters, { library(multiblock) })
  
  parLapply(clusters, 1:n_simulation_datasets, function(i) {
    process_disco_sca_function(i, data_directory, data_type, save_directory, save_type)
  })
  
  stopCluster(clusters)
}

##### Obtain DISCO-SCA Models for the Simulated Data Outlier Contamination #####
set.seed(2025)
n_simulation_datasets <- 10
outlier_levels <- c("01", "05", "10", "25", "50")

for (outlier_level in outlier_levels) {
  cat("Processing outlier level:", outlier_level, "\n")
  
  # Here we use clusters to simultaneously compute multiple models
  n_cores <- min(10, detectCores() - 1)
  clusters <- makeCluster(n_cores)
  
  # Settings for reading the uncontaminated data and saving the contaminated data
  data_directory <- file.path(project_directory, "Simulated_data/Outlier", paste0("Level_", outlier_level))
  data_type <- paste0("outlier_data_", outlier_level, "_simulation_")
  save_directory <- file.path(project_directory, "Results/Outlier_Results", paste0("Level_", outlier_level))
  save_type <- paste0("disco_sca_outlier_", outlier_level, "_simulation_")
  
  clusterExport(clusters, c("data_directory", "data_type", "save_directory", "save_type", "process_disco_sca_function"))
  clusterEvalQ(clusters, { library(multiblock) })
  
  parLapply(clusters, 1:n_simulation_datasets, function(i) {
    process_disco_sca_function(i, data_directory, data_type, save_directory, save_type)
  })
  
  stopCluster(clusters)
}
