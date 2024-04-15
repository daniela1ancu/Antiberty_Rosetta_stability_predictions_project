#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import shutil
import pandas as pd
import pickle as pkl
import numpy as np


# ## Sort files
# The original data comes in the folder "pseudo_log_likelihood" which contain both the H-chain and the L-chain pseudologs. We first need to sort this into only the H-chain and L-chain folders to make it easier to work with.

# In[ ]:


input_folder = "Initial data/pseudo_log_likelihood_csv_melted_scaled"
output_folder_h = "Antiberty-h"
output_folder_l = "Antiberty-l"

files = os.listdir(input_folder)  # retrieves the filenames in the folder

for file in files:
    if "_h_" in file:
        source_path = os.path.join(input_folder,file)  # creates the source path
        destination_path = os.path.join(output_folder_h,file)  # creates the destination path
        shutil.move(source_path, destination_path)  # moves the file from source to destination
    else:
        source_path = os.path.join(input_folder,file)
        destination_path = os.path.join(output_folder_l,file)
        shutil.move(source_path, destination_path)

