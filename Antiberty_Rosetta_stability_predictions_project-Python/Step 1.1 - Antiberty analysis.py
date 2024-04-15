#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import pandas as pd
import numpy as np
import ast
import re
import csv


# ## Load antiberty files, filter them

# In[ ]:


# Define the directory containing your pickle files
directory_h = "Antiberty-h"

# List all files in the directory
files = os.listdir(directory_h)

def read_files_return_dataframes(files):
    antiberty_h_dataframes = {}
    
    """
    take the files from the directory and first read and extract their custom names
    then we drop the first 5 rows since they are useless 
    we reset the index (we need this later when merging with the MUT_AA dataframe
    """
    
    for file in files:
        file_path = os.path.join(directory_h, file)  # adds files to the directory wanted
        
        if os.path.isfile(file_path):
            custom_name = file[0:7]  # extracts the custom name e.g.7n4l_HL
            with open(file_path, 'r') as file:
                loaded_data = pd.read_csv(file_path)

            loaded_data = loaded_data.rename(columns = {"wt_aa": "WT_AA", "scaled_pseudolog_likelihood": "Scaled Antiberty Score", "mut_aa": "MUT_AA", "PDB_numbering": "position"})
            antiberty_h_dataframes[custom_name] = loaded_data

    return antiberty_h_dataframes

antiberty_h_dataframes = read_files_return_dataframes(files)


# In[ ]:


# Define the directory containing your pickle files
directory_l = "Antiberty-l"

# List all files in the directory
files = os.listdir(directory_l)

def read_files_return_dataframes(files):
    antiberty_l_dataframes = {}

    for file in files:
        file_path = os.path.join(directory_l, file)
        if os.path.isfile(file_path):
            custom_name = file[0:7]
            with open(file_path, 'rb') as f:
                loaded_data = pd.read_csv(f)
                
            loaded_data = loaded_data.rename(columns = {"wt_aa": "WT_AA", "scaled_pseudolog_likelihood": "Scaled Antiberty Score", "mut_aa": "MUT_AA", "PDB_numbering": "position"})
            antiberty_l_dataframes[custom_name] = loaded_data

    return antiberty_l_dataframes

antiberty_l_dataframes = read_files_return_dataframes(files)


# ## Save the antiberty files

# In[ ]:


output_folder_H = "AntiBERTy filtered files/H-AntiBERTy"
output_folder_L = "AntiBERTy filtered files/L-AntiBERTy"

def antiberty_to_csv(dictionary, output_folder):
    for custom_name, df in dictionary.items():
        file_path = os.path.join(output_folder, f"{custom_name}.csv")
        df.to_csv(file_path, index = False)

antiberty_to_csv(antiberty_h_dataframes, output_folder_H)
antiberty_to_csv(antiberty_l_dataframes, output_folder_L)

