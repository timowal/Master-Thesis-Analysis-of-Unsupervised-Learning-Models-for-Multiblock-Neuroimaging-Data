'''
Author: Timo Waling (581706tw)
Converts 3D REHO,LFCD,and ALFF files into useable ROI data for processing in R
'''

import nibabel as nib
import os
import pandas as pd
from nilearn.input_data import NiftiLabelsMasker

# List with the patient groups
groups = ["ASD_group","Control_group","Combined"]
# List with all brain data types
data_types = ['reho','lfcd','alff']
# Import atlas file for converting 3D files
atlas_dir = "data/atlas/ho_roi_atlas.nii.gz"
atlas_img = nib.load(atlas_dir)
# Create a masker using the atlas
masker = NiftiLabelsMasker(labels_img=atlas_img, standardize=True)
# Create a csv file for each group and each data type
for i in groups:
    for j in data_types:
        # Get directory for all patient files
        cur_dir = f"data/{i}/{j}/"
        # Create a list of files to iterate through
        file_list = [os.path.join(cur_dir, f) for f in os.listdir(cur_dir) if f.endswith(".nii.gz")]
        # Save converted data to be put into csv
        all_cur_values = []
        for k in file_list:
            # Transform file into ROI data
            cur_img = nib.load(k)
            cur_values = masker.fit_transform(cur_img)
            # Check if masked values are the correct shape
            if cur_values.shape[0] > 0 and cur_values.shape[1] > 0:
                print(f"Extraction successful for {os.path.basename(cur_values)}! Shape: {cur_values.shape}")
                all_cur_values.append(cur_values.flatten())
            else:
                print(f"Warning: No data extracted from {cur_values}")
        # Convert the list to a DataFrame for saving to csv
        cur_df = pd.DataFrame(all_cur_values)
        # Show shape and preview
        print(f"\nFinal DataFrame shape: {cur_df.shape}")
        print(cur_df.head())
        # Save current results dataframe into csv
        cur_df.to_csv(f"data/Outputs/{j}_atlas_{i}.csv", index=False)
