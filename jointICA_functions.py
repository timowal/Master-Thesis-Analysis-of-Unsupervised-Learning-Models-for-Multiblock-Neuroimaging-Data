'''
Author: Timo Waling (581706tw)
Python Code Containing the functions used for jointICA.py
'''

# Import the packages and functions that are used in this python script
import pandas as pd
import numpy as np
from sklearn.decomposition import FastICA
from sklearn.preprocessing import StandardScaler

# Define a function to load and clean each dataset
def load_data(path):
    """
    Read the .csv data files and remove the first row of the data (which contains a row of zeros)\

    Arguments:
        - param path: The path to the .csv file that contains the data to be loaded

    Returns:
        - The data as a np array
    """
    data = pd.read_csv(path, header=None)
    data = data.iloc[1:, :]  # remove first row (all zeros)
    data = data.to_numpy(dtype=float)

    return data

def explained_variation_for_modality(original_data, reconstructed_data):
    """
    Compute and return the ratio of explained variation compared to the original data

    Arguments:
        - param original_data: The original multiblock data
        - param reconstructed_data: Data obtained through the joint ICA method

    Returns:
        - The explained variation compared to the original data
    """
    total_var_mod = np.sum(np.var(original_data, axis=0))
    recon_var_mod = np.sum(np.var(reconstructed_data, axis=0))

    return recon_var_mod / total_var_mod

def joint_ica_explained_variation(data_reho, data_lfcd, data_alff, n_comp):
    """
    Perform joint ICA on the REHO, LFCD and ALFF data blocks. Afterward compute explained
    variation for the data overall and per data block.

    Arguments:
        - param data_reho: Data for the REHO modality
        - param data_lfcd: Data for the LFCD modality
        - param data_alff: Data for the ALFF modality
        - param n_comp: amount of components for the joint ICA method

    Return:
        - explained_var_ratio: total explained variation of the multiblock data
        - explained_var_block: explained variation per data block
    """

    # Load the scaler and standardize the data in the data blocks separately
    scaler = StandardScaler()
    reho_standardized = scaler.fit_transform(data_reho)
    lfcd_standardized = scaler.fit_transform(data_lfcd)
    alff_standardized = scaler.fit_transform(data_alff)

    # Concatenate the standardized data along the features
    joint_data = np.hstack([reho_standardized, lfcd_standardized, alff_standardized])

    # Fit the joint ICA model and extract important characteristics
    joint_ica = FastICA(n_components=n_comp, random_state=0, max_iter=5000)

    # Reconstruct the data
    sources = joint_ica.fit_transform(joint_data)
    mixing = joint_ica.mixing_
    reconstructed_data = np.dot(sources, mixing.T)

    # Explained variation total
    total_original_var = np.sum(np.var(joint_data, axis=0))
    total_reconstructed_var = np.sum(np.var(reconstructed_data, axis=0))
    explained_var_ratio = total_reconstructed_var / total_original_var

    # Index slices for each data block
    n_reho = reho_standardized.shape[1]
    n_lfcd = lfcd_standardized.shape[1]
    n_alff = alff_standardized.shape[1]

    index_reho = slice(0, n_reho)
    index_lfcd = slice(n_reho, n_reho + n_lfcd)
    index_alff = slice(n_reho + n_lfcd, n_reho + n_lfcd + n_alff)

    # Compute the explained variation per data block
    explained_var_block = {
        'REHO': explained_variation_for_modality(joint_data[:, index_reho], reconstructed_data[:, index_reho]),
        'LFCD': explained_variation_for_modality(joint_data[:, index_lfcd], reconstructed_data[:, index_lfcd]),
        'ALFF': explained_variation_for_modality(joint_data[:, index_alff], reconstructed_data[:, index_alff]),
    }

    # Return the total explained variation and the explained variation per data block
    return explained_var_ratio, explained_var_block


def run_joint_ica_on_simulated_data(data_directory, file_name, n_datasets, n_comp):
    """
    Function that performs joint ICA for multiple realisations of simulated data sets.

    Arguments:
        - param data_directory: The path to the folder containing the datasets to cycle through
        - param file_name: The name of the files containing the simulated/contaminated data
        - param n_datasets: The amount of realisations of simulated/contaminated data to cycle through
        - param n_comp: The amount of components for the joint ICA method

    Returns:
        - average_total_explained_var: Average ratio of explained variation over the datasets
        - average_explained_var_modality: average explained variation per type of data block over the datasets
    """

    total_explained_var_list = []
    explained_var_modality_list = []

    for i in range(1, n_datasets + 1):
        # Create the path to the file of this realisation of the simulated/contaminated data, and load the data
        file_path = f"{data_directory}{file_name.format(i)}"
        data = load_data(file_path)

        # Extract the REHO, LFCD and ALFF data from the dataset
        reho_data = data[:, 0:111]
        lfcd_data = data[:, 111:222]
        alff_data = data[:, 222:333]

        # Perform joint ICA for this realisation of the simulated/contaminated data
        total_explained_var, explained_var_modality = joint_ica_explained_variation(reho_data, lfcd_data, alff_data, n_comp)

        # Print the individual results of joint ICA for this realisation of the data
        print(f"File: {file_name.format(i)}")
        print(f"Total explained variation: {total_explained_var:.3f}")
        print("Explained variation by modality:")
        for modality, explained_var in explained_var_modality.items():
            print(f"{modality}: {explained_var:.3f}")

        # Add the results of this realisation to the lists
        total_explained_var_list.append(total_explained_var)
        explained_var_modality_list.append(explained_var_modality)

    # Compute averages
    average_total_explained_var = sum(total_explained_var_list) / n_datasets

    # Initialize sums for each modality
    sum_modality = {modality: 0 for modality in explained_var_modality_list[0].keys()}
    for explained_var_dict in explained_var_modality_list:
        for modality, value in explained_var_dict.items():
            sum_modality[modality] += value
    average_explained_var_modality = {modality: sum_values / n_datasets for modality, sum_values in sum_modality.items()}

    return average_total_explained_var, average_explained_var_modality
