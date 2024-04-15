#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import shutil
import pandas as pd
import numpy as np


# ## Read the files and sort them
# the files are names something like 1rzi_NM_N or 1rzi_NM_M.
# 
# The first letter is always the heavy chain whilst the second is always the ligth chain.

# In[2]:


os.makedirs("Surface information/H chain", exist_ok=True)
os.makedirs("Surface information/L chain", exist_ok=True)

input_folder = "Initial data/iso_rpopsResidue"
surface_info_dataframes = {}

def sort_information(input_folder):
    """
    for a file named 1rzi_NM_N or 1rzi_NM_M:
    file[8:9] = NM
    file[5:6] = N (heavy chain)
    file[6:7] = M (light chain)
    this function is used to copy files into a new folder and rename them to contain the first 6 characters in the file: 1rzi_NM
    """
    
    for file in os.listdir(input_folder):
        file_path = os.path.join(input_folder, file)
        # print(file[5:6], file[6:7], file[:7])
        
        if file[8:9] == file[5:6]:
            # Copy file to H chain directory
            custom_name = file[:7]
            destination_path = os.path.join("Surface information/H chain", custom_name)
            shutil.copyfile(file_path, destination_path)

        elif file[8:9] == file[6:7]:
            # Copy file to L chain directory
            custom_name = file[:7]
            destination_path = os.path.join("Surface information/L chain", custom_name)
            shutil.copyfile(file_path, destination_path)
            
        else:
            print(file_path)

sort_information(input_folder)


# ## Convert to dataframes
# 
# the text files will be converted to dataframes, which will then be concatenated so each file will contain both the heavy and the light chains. 

# In[5]:


heavy_chain_input_folder = "Surface information/H chain"
light_chain_input_folder = "Surface information/L chain"
HC_dataframe = {}
LC_dataframe = {}
conc_dataframes = {}

def read_files(input_folder):
    """
    this function takes the files converts them to dataframes and then formats the columns 
    NOTE - had to do it this way since the format of the file was a bit off.
    """

    
    files = os.listdir(input_folder)
    for file in files:
        file_path = os.path.join(input_folder, file)
        custom_name = file[:7]

        column_names = ["ResidNe", "Chain", "ResidNr", "iCode", "Phob/A^2", "Phil/A^2", "SASA/A^2", "Q(SASA)", "N(overl)", "Surf/A^2"]
        
        with open(file_path, "r") as f:
            # Read dataframe from file with specified column names and delimiter
            df = pd.read_csv(f, delimiter="\t", names=column_names, skiprows=1)
        
        if input_folder == heavy_chain_input_folder:
            # df = df.rename(columns={"Unnamed: 5": "Phil/A^2"})
            df["Chain"] = "H"
            HC_dataframe[custom_name] = df
        elif input_folder == light_chain_input_folder:
            df["Chain"] = "L"
            LC_dataframe[custom_name] = df

read_files(heavy_chain_input_folder)
read_files(light_chain_input_folder)


# In[6]:


for custom_name, df in HC_dataframe.items():
    if custom_name in LC_dataframe:
        L_chain = LC_dataframe[custom_name]
        conc_dataframes[custom_name] = pd.concat([df, L_chain])
    else:
        print(custom_name)


"""
i had a look and it appears that these names are not in the "clean files".
"""


# ## Add a column with the one letter code
# 
# This will help us successfuly merge these files with the "clean files"

# In[ ]:


"""
started off with only one dataframe
"""
amino_acids_mapping = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V"
}

dataframe = conc_dataframes["7mn8_DC"].copy()

dataframe["WT_AA"] = dataframe["ResidNe"].str.strip()
dataframe["WT_AA"] = dataframe["WT_AA"].map(amino_acids_mapping)


# In[8]:


# Define the directory path where CSV files will be saved
file_path = "Surface information/surface_access_clean_files"
os.makedirs(file_path, exist_ok=True)

# Mapping of three-letter amino acid codes to one-letter amino acid codes
amino_acids_mapping = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V"
}

# Dictionary to store modified DataFrames (not used in the current code)
solvent_access_info = {}

def add_one_letter_code():
    """
    Function to add a new column containing one-letter amino acid codes to each DataFrame.
    
    Creates a copy of the DataFrame to avoid modifying the original.
    Strips any whitespace from the "ResidNe" column values and adds them to a new column "WT_AA".
    Maps the three-letter amino acid codes to one-letter codes using the defined mapping.
    Constructs the file path for saving the DataFrame to a CSV file.
    Saves the modified DataFrame to a CSV file without the index.
    Optionally, stores the modified DataFrame in the dictionary (currently not used).
    """

    
    for custom_name, df in conc_dataframes.items():
        dataframe = df.copy()
        dataframe["WT_AA"] = dataframe["ResidNe"].str.strip()
        dataframe["WT_AA"] = dataframe["WT_AA"].map(amino_acids_mapping)
        dataframe_file_path = os.path.join(file_path, f"{custom_name}.csv")
        dataframe.to_csv(dataframe_file_path, index=False)
        # solvent_access_info[custom_name] = dataframe
        

add_one_letter_code()

