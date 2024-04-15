#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import numpy as np
import math


# ## Merge dataframes
# 
# Read the files which contain the solvent accessible information and then also the files which contain the Rosetta and AntiBERTy scores.

# In[12]:


solvent_access_folder = "Surface information/surface_access_clean_files"
dataframe_folder = "Clean files"

solvent_access_info = {}
dataframes = {}

def read_files(input_folder, output_dictionary):
    for file in os.listdir(input_folder):
        file_path = os.path.join(input_folder, file)
        
        if file.endswith(".csv"):
            custom_name = file[:7]
            with open(file_path, "r") as f:
                df = pd.read_csv(f)
                output_dictionary[custom_name] = df

read_files(solvent_access_folder, solvent_access_info)
read_files(dataframe_folder, dataframes)


# ## Create a column named PDB_num which contains the necessary info for merging

# In[13]:


"""
for the solvent accessible information the icode contains a dash line "-" which would need to be replaced to an empty space, 
since the "clean files" contain "NaN" (empty) values.
then we create a column named PDB_num which would contain the PDB numbering and if applicable, add the insertion code (e.g. 108A or 108B)
"""
for custom_name, df in solvent_access_info.items():
    df = df.copy()
    df["iCode"] = df["iCode"].replace({"-": ""})
    df["PDB_num"] = df["ResidNr"].astype(str) + df["iCode"]
    df["PDB_num"] = df["PDB_num"].str.strip()
    solvent_access_info[custom_name] = df


# In[14]:


"""
for this I just conver the position to a string (since now it is an int)
"""

for custom_name, df in dataframes.items():
    df = df.copy()
    # df["insertion_code"] = df["insertion_code"].fillna("-")
    df["position"] = df["position"].astype(str)
    dataframes[custom_name] = df


# ## Merge and save the files 

# In[65]:


"""
lastly combine the two dataframes. The PDB_num and position would correcly match the two dataframes based on the amino acid type,
the position and the chain type.
We sort them (so they are alphabetically ordered) and then save the files in the desired folder.
"""

output_folder = "Surface information"
os.makedirs(output_folder, exist_ok=True)

merged_dataframes = {}

for file_name, df in dataframes.items():
    if file_name in solvent_access_info:
        custom_name = file_name[:7]
        solvent_access_df = solvent_access_info[custom_name]
        merged_dataframe = df.merge(solvent_access_df, 
                                        left_on=["WT_AA", "position", "chain"], 
                                        right_on=["WT_AA", "PDB_num", "Chain"], 
                                        how="left")
        merged_dataframe.sort_values(by="IMGT numbering", ascending=True)

        file_path = f"{output_folder}/clean_files_solvent_access_info/{custom_name}.csv"
        merged_dataframe.to_csv(file_path, index=False)
        
        merged_dataframes[custom_name] = merged_dataframe
        


# ## Split rows into dataframes - core and surface (one example)
# 
# take one dataframe (Trastuzumab) and split the dataframes into core and surface rows, depending on their QSASA score.
# 
# > <0.15 - they are considered core
# 
# > larger 0.15 - considered surface

# In[71]:


dataframe = merged_dataframes["7mn8_DC"].copy()

core_residues = []
surface_residues = []
custom_name = "7mn8_DC"
dataframe["Q(SASA)"] = dataframe["Q(SASA)"].astype(float)

for index, row in dataframe.iterrows():
    if not math.isnan(row["Q(SASA)"]):
        if row["Q(SASA)"] < 0.15:
            core_residues.append({
                "chain": row["chain"], 
                "position": row["position"],
                "IMGT numbering": row["IMGT numbering"],
                "WT_AA": row["WT_AA"],
                "MUT_AA": row["MUT_AA"],
                "region": row["region"],
                "Rosetta score": row["Rosetta normalised"],
                "AntiBERTy scores": row["AntiBERTy normalised"],  
                "QSASA": row["Q(SASA)"]
            })
        elif row["Q(SASA)"] >= 0.15:
            surface_residues.append({
                "chain": row["chain"], 
                "position": row["position"],
                "IMGT numbering": row["IMGT numbering"],
                "WT_AA": row["WT_AA"],
                "MUT_AA": row["MUT_AA"],
                "region": row["region"],
                "Rosetta score": row["Rosetta normalised"],
                "AntiBERTy scores": row["AntiBERTy normalised"],  
                "QSASA": row["Q(SASA)"]
            })
        else:
            print(f"{row} has no QSASA")
                
    dataframe_core = pd.DataFrame(core_residues)
    dataframe_core.to_csv(f"{custom_name}_Core.csv", index=False) 
            
    dataframe_surface = pd.DataFrame(surface_residues)
    dataframe_surface.to_csv(f"{custom_name}_surface.csv", index=False)


# ## Split all dataframes depeding on their QSASA score

# In[3]:


input_folder = "Surface information/clean_files_solvent_access_info"
output_folder = "Surface information"
os.makedirs(f"{output_folder}/core_residues_dataframes", exist_ok=True)
os.makedirs(f"{output_folder}/surface_residues_dataframes", exist_ok=True)

for file_name in os.listdir(input_folder):
    file_path = os.path.join(input_folder, file_name)
    if file_name.endswith(".csv"):
        df = pd.read_csv(file_path)
        custom_name = file_name[:7]
        
        core_residues = []
        surface_residues = []
        df["Q(SASA)"] = df["Q(SASA)"].astype(float)
    
    
        for index, row in df.iterrows():
            if not math.isnan(row["Q(SASA)"]):
                Qsasa = row["Q(SASA)"]
                
                if Qsasa < 0.15:
                    core_residues.append({
                        "chain": row["chain"], 
                        "position": row["position"],
                        "IMGT numbering": row["IMGT numbering"],
                        "WT_AA": row["WT_AA"],
                        "MUT_AA": row["MUT_AA"],
                        "region": row["region"],
                        "Rosetta score": row["Rosetta normalised"],
                        "AntiBERTy scores": row["AntiBERTy normalised"],  
                        "QSASA": Qsasa
                    })
                elif Qsasa >= 0.15:
                    surface_residues.append({
                        "chain": row["chain"], 
                        "position": row["position"],
                        "IMGT numbering": row["IMGT numbering"],
                        "WT_AA": row["WT_AA"],
                        "MUT_AA": row["MUT_AA"],
                        "region": row["region"],
                        "Rosetta score": row["Rosetta normalised"],
                        "AntiBERTy scores": row["AntiBERTy normalised"],  
                        "QSASA": Qsasa
                    })
                    
        dataframe_core = pd.DataFrame(core_residues)
        dataframe_core.to_csv(f"{output_folder}/core_residues_dataframes/{custom_name}.csv", index=False) 
                
        dataframe_surface = pd.DataFrame(surface_residues)
        dataframe_surface.to_csv(f"{output_folder}/surface_residues_dataframes/{custom_name}.csv", index=False)


# In[4]:


def process_file(input_file, output_folder):
    df = pd.read_csv(input_file)
    custom_name = os.path.basename(input_file)[:7]
    df["Q(SASA)"] = df["Q(SASA)"].astype(float)

    core_residues = df[df["Q(SASA)"] < 0.15]
    surface_residues = df[df["Q(SASA)"] >= 0.15]

    core_residues.to_csv(f"{output_folder}/core_residues_dataframes/{custom_name}.csv", index=False)
    surface_residues.to_csv(f"{output_folder}/surface_residues_dataframes/{custom_name}.csv", index=False)

def process_files(input_folder, output_folder):
    os.makedirs(f"{output_folder}/core_residues_dataframes", exist_ok=True)
    os.makedirs(f"{output_folder}/surface_residues_dataframes", exist_ok=True)

    for file_name in os.listdir(input_folder):
        if file_name.endswith(".csv"):
            file_path = os.path.join(input_folder, file_name)
            process_file(file_path, output_folder)

input_folder = "Surface information/clean_files_solvent_access_info"
output_folder = "Surface information"
process_files(input_folder, output_folder)

