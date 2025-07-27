### Code For Creating all JIVE models
# Author: Timo Waling (581706tw)


## Install packages and load them:
install.packages("multiblock")
install.packages("neuroim")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("r.jive")
install.packages("reshape2")
install.packages("parallel")

library(multiblock)
library(neuroim)
library(dplyr)
library(ggplot2)
library(r.jive)
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


##### Obtain JIVE models for the real data #####
data_blocks_real_control <- list(reho = t(reho_data_control), lfcd = t(lfcd_data_control), alff = t(alff_data_control))
data_blocks_real_asd <- list(reho = t(reho_data_asd), lfcd = t(lfcd_data_asd), alff = t(alff_data_asd))
data_blocks_real_combined <- list(reho = t(reho_data_combined), lfcd = t(lfcd_data_combined), alff = t(alff_data_combined))

### JIVE for the data
jive_real_control <- jive(data_blocks_real_control, method = "given", rankJ = 3, rankA = c(2, 2, 2))
jive_real_asd <- jive(data_blocks_real_asd, method = "given", rankJ = 3, rankA = c(2, 2, 2))
jive_real_combined <- jive(data_blocks_real_combined, method = "given", rankJ = 3, rankA = c(2, 2, 2))

# Save the JIVE models for later use
save_directory <- file.path(project_directory, "Results/Real_Data")
saveRDS(jive_real_control, file = file.path(save_directory, "jive_real_control.rds"))
saveRDS(jive_real_asd, file = file.path(save_directory, "jive_real_asd.rds"))
saveRDS(jive_real_combined, file = file.path(save_directory, "jive_real_combined.rds"))


##### Function for processing JIVE for given data #####
process_jive_function <- function(i, data_directory, data_type, save_directory, save_type) {
  cat("Processing dataset", i, "from", data_type, "...\n")
  
  # Read the data and ensure it is numeric
  filename <- file.path(data_directory, paste0(data_type, i, ".csv"))
  full_data <- read.csv(filename, header = FALSE)[-1, ]
  full_data[] <- lapply(full_data, function(x) as.numeric(as.character(x)))
  
  data_reho <- full_data[, 1:111]
  data_lfcd <- full_data[, 112:222]
  data_alff <- full_data[, 223:333]
  
  data_blocks_full <- list(
    reho = t(data_reho),
    lfcd = t(data_lfcd),
    alff = t(data_alff)
  )
  
  # Obtain the JIVE model for the data and save it
  jive_model <- jive(data_blocks_full, method = "given", rankJ = 3, rankA = c(2, 2, 2))
  saveRDS(jive_model, file = file.path(save_directory, paste0(save_type, i, ".rds")))
  
  cat("Finished dataset", i, "\n")
}


##### Obtain JIVE Models for the Simulated Data without contamination #####
set.seed(2025)
n_simulation_datasets <- 10

# Settings for reading the uncontaminated data and saving the contaminated data
data_directory <- file.path(project_directory, "Simulated_data/Uncontaminated")
data_type <- "simulated_data_"
save_directory <- file.path(project_directory, "Results/Simulated_Data")
save_type <- "jive_simulated_"

# Here we use clusters to simultaneously compute multiple models
n_cores <- min(10, detectCores() - 1)
clusters <- makeCluster(n_cores)

clusterExport(clusters, c("data_directory", "data_type", "save_directory",
                          "save_type", "process_jive_function"))
clusterEvalQ(clusters, { library(r.jive) })

parLapply(clusters, 1:n_simulation_datasets, function(i) {
  process_jive_function(i, data_directory, data_type, save_directory, save_type)
})

stopCluster(clusters)



##### Obtain JIVE Models for the Simulated Data Noise Contamination #####
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
  save_type <- paste0("jive_noise_", noise_level, "_simulation_")
  
  clusterExport(clusters, c("data_directory", "data_type", "save_directory", "save_type", "process_jive_function"))
  clusterEvalQ(clusters, { library(r.jive) })
  
  parLapply(clusters, 1:n_simulation_datasets, function(i) {
    process_jive_function(i, data_directory, data_type, save_directory, save_type)
  })
  
  stopCluster(clusters)
}

##### Obtain JIVE Models for the Simulated Data Outlier Contamination #####
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
  save_type <- paste0("jive_outlier_", outlier_level, "_simulation_")
  
  clusterExport(clusters, c("data_directory", "data_type", "save_directory", "save_type", "process_jive_function"))
  clusterEvalQ(clusters, { library(r.jive) })
  
  parLapply(clusters, 1:n_simulation_datasets, function(i) {
    process_jive_function(i, data_directory, data_type, save_directory, save_type)
  })
  
  stopCluster(clusters)
}
