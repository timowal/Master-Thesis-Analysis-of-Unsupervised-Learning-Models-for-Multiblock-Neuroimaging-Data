'''
Author: Timo Waling (581706tw)
Python Code For Joint ICA Method
'''

# Import the packages and functions that are used in this python script
import pandas as pd
import numpy as np
import random
from sklearn.decomposition import FastICA
from sklearn.preprocessing import StandardScaler
from jointICA_functions import (load_data, explained_variation_for_modality, joint_ica_explained_variation,
                                run_joint_ica_on_simulated_data)

# Set seed for consistent reproduction of results
random.seed(2025)

# Number of realizations of the simulated data
n_datasets = 10

# Set number of components to be used for jointICA
n_comp = 9

# Settings for which data to run joint ICA for
real_data_combined = False
real_data_control = False
real_data_asd = False
simulated_data_clean = False
noise_contamination_01 = False
noise_contamination_05 = False
noise_contamination_10 = False
noise_contamination_25 = False
noise_contamination_50 = False
outlier_contamination_01 = False
outlier_contamination_05 = False
outlier_contamination_10 = False
outlier_contamination_25 = False
outlier_contamination_50 = True

# Run joint ICA on the full real data if set to TRUE
if real_data_combined:
    # Load each dataset
    reho_data = load_data('data/reho_data/reho_atlas_combined.csv')
    lfcd_data = load_data('data/lfcd_data/lfcd_atlas_combined.csv')
    alff_data = load_data('data/alff_data/alff_atlas_combined.csv')

    explained_total, explained_modalities = joint_ica_explained_variation(reho_data, lfcd_data, alff_data, n_comp)

    print("Results for Real Data Combined:")
    print(f"Total explained variation: {explained_total:.3f}")
    print("Explained variation per modality:")
    for modality, value in explained_modalities.items():
        print(f"{modality}: {value:.3f}")


# Run joint ICA for the real Control group data if set to TRUE
if real_data_control:
    # Load each dataset
    reho_data = load_data('data/reho_data/reho_atlas_Control_group.csv')
    lfcd_data = load_data('data/lfcd_data/lfcd_atlas_Control_group.csv')
    alff_data = load_data('data/alff_data/alff_atlas_Control_group.csv')

    explained_total, explained_modalities = joint_ica_explained_variation(reho_data, lfcd_data, alff_data, n_comp)

    print("Results for Real Data Control Group:")
    print(f"Total explained variation: {explained_total:.3f}")
    print("Explained variation per modality:")
    for modality, value in explained_modalities.items():
        print(f"{modality}: {value:.3f}")


# Run joint ICA for the real ASD data if set to TRUE
if real_data_asd:
    # Load each dataset
    reho_data = load_data('data/reho_data/reho_atlas_ASD_group.csv')
    lfcd_data = load_data('data/lfcd_data/lfcd_atlas_ASD_group.csv')
    alff_data = load_data('data/alff_data/alff_atlas_ASD_group.csv')

    explained_total, explained_modalities = joint_ica_explained_variation(reho_data, lfcd_data, alff_data, n_comp)

    print("Results for Real Data ASD Group:")
    print(f"Total explained variation: {explained_total:.3f}")
    print("Explained variation per modality:")
    for modality, value in explained_modalities.items():
        print(f"{modality}: {value:.3f}")



# Run joint ICA for the uncontaminated simulated data if set to TRUE
if simulated_data_clean:
    data_directory = 'Simulated_data/Uncontaminated/'
    file_name = 'simulated_data_{}.csv'

    average_total, average_per_modality = run_joint_ica_on_simulated_data(data_directory, file_name, n_datasets, n_comp)
    print(f"Average total explained variation over {n_datasets} files: {average_total:.3f}")
    print("Average explained variation by modality:")
    for modality, value in average_per_modality.items():
        print(f"{modality}: {value:.3f}")


# Run joint ICA for the simulated data at 1% noise contamination if set to TRUE
if noise_contamination_01:
    data_directory = 'Simulated_data/Noise/Level_01/'
    file_name = 'noise_data_01_simulation_{}.csv'

    average_total, average_per_modality = run_joint_ica_on_simulated_data(data_directory, file_name, n_datasets, n_comp)
    print(f"Average total explained variation over {n_datasets} files: {average_total:.3f}")
    print("Average explained variation by modality:")
    for modality, value in average_per_modality.items():
        print(f"{modality}: {value:.3f}")


# Run joint ICA for the simulated data at 5% noise contamination if set to TRUE
if noise_contamination_05:
    data_directory = 'Simulated_data/Noise/Level_05/'
    file_name = 'noise_data_05_simulation_{}.csv'

    average_total, average_per_modality = run_joint_ica_on_simulated_data(data_directory, file_name, n_datasets, n_comp)
    print(f"Average total explained variation over {n_datasets} files: {average_total:.3f}")
    print("Average explained variation by modality:")
    for modality, value in average_per_modality.items():
        print(f"{modality}: {value:.3f}")


# Run joint ICA for the simulated data at 10% noise contamination if set to TRUE
if noise_contamination_10:
    data_directory = 'Simulated_data/Noise/Level_10/'
    file_name = 'noise_data_10_simulation_{}.csv'

    average_total, average_per_modality = run_joint_ica_on_simulated_data(data_directory, file_name, n_datasets, n_comp)
    print(f"Average total explained variation over {n_datasets} files: {average_total:.3f}")
    print("Average explained variation by modality:")
    for modality, value in average_per_modality.items():
        print(f"{modality}: {value:.3f}")


# Run joint ICA for the simulated data at 25% noise contamination if set to TRUE
if noise_contamination_25:
    data_directory = 'Simulated_data/Noise/Level_25/'
    file_name = 'noise_data_25_simulation_{}.csv'

    average_total, average_per_modality = run_joint_ica_on_simulated_data(data_directory, file_name, n_datasets, n_comp)
    print(f"Average total explained variation over {n_datasets} files: {average_total:.3f}")
    print("Average explained variation by modality:")
    for modality, value in average_per_modality.items():
        print(f"{modality}: {value:.3f}")


# Run joint ICA for the simulated data at 50% noise contamination if set to TRUE
if noise_contamination_50:
    data_directory = 'Simulated_data/Noise/Level_50/'
    file_name = 'noise_data_50_simulation_{}.csv'

    average_total, average_per_modality = run_joint_ica_on_simulated_data(data_directory, file_name, n_datasets, n_comp)
    print(f"Average total explained variation over {n_datasets} files: {average_total:.3f}")
    print("Average explained variation by modality:")
    for modality, value in average_per_modality.items():
        print(f"{modality}: {value:.3f}")


# Run joint ICA for the simulated data at 1% outlier contamination if set to TRUE
if outlier_contamination_01:
    data_directory = 'Simulated_data/Outlier/Level_01/'
    file_name = 'outlier_data_01_simulation_{}.csv'

    average_total, average_per_modality = run_joint_ica_on_simulated_data(data_directory, file_name, n_datasets, n_comp)
    print(f"Average total explained variation over {n_datasets} files: {average_total:.3f}")
    print("Average explained variation by modality:")
    for modality, value in average_per_modality.items():
        print(f"{modality}: {value:.3f}")


# Run joint ICA for the simulated data at 5% outlier contamination if set to TRUE
if outlier_contamination_05:
    data_directory = 'Simulated_data/Outlier/Level_05/'
    file_name = 'outlier_data_05_simulation_{}.csv'

    average_total, average_per_modality = run_joint_ica_on_simulated_data(data_directory, file_name, n_datasets, n_comp)
    print(f"Average total explained variation over {n_datasets} files: {average_total:.3f}")
    print("Average explained variation by modality:")
    for modality, value in average_per_modality.items():
        print(f"{modality}: {value:.3f}")


# Run joint ICA for the simulated data at 10% outlier contamination if set to TRUE
if outlier_contamination_10:
    data_directory = 'Simulated_data/Outlier/Level_10/'
    file_name = 'outlier_data_10_simulation_{}.csv'

    average_total, average_per_modality = run_joint_ica_on_simulated_data(data_directory, file_name, n_datasets, n_comp)
    print(f"Average total explained variation over {n_datasets} files: {average_total:.3f}")
    print("Average explained variation by modality:")
    for modality, value in average_per_modality.items():
        print(f"{modality}: {value:.3f}")


# Run joint ICA for the simulated data at 25% outlier contamination if set to TRUE
if outlier_contamination_25:
    data_directory = 'Simulated_data/Outlier/Level_25/'
    file_name = 'outlier_data_25_simulation_{}.csv'

    average_total, average_per_modality = run_joint_ica_on_simulated_data(data_directory, file_name, n_datasets, n_comp)
    print(f"Average total explained variation over {n_datasets} files: {average_total:.3f}")
    print("Average explained variation by modality:")
    for modality, value in average_per_modality.items():
        print(f"{modality}: {value:.3f}")

# Run joint ICA for the simulated data at 50% outlier contamination if set to TRUE
if outlier_contamination_50:
    data_directory = 'Simulated_data/Outlier/Level_50/'
    file_name = 'outlier_data_50_simulation_{}.csv'

    average_total, average_per_modality = run_joint_ica_on_simulated_data(data_directory, file_name, n_datasets, n_comp)
    print(f"Average total explained variation over {n_datasets} files: {average_total:.3f}")
    print("Average explained variation by modality:")
    for modality, value in average_per_modality.items():
        print(f"{modality}: {value:.3f}")
