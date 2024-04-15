#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import pandas as pd


# ## open antiberty filtered and rosetta files

# In[ ]:


output_folder_H = "AntiBERTy filtered files/H-AntiBERTy"
output_folder_L = "AntiBERTy filtered files/L-AntiBERTy"

# Create dictionaries to store DataFrames
antiberty_h_dataframes = {}
antiberty_l_dataframes = {}

# Read CSV files for H-chain
for custom_name in os.listdir(output_folder_H):
    file_path = os.path.join(output_folder_H, custom_name)
    if os.path.isfile(file_path) and file_path.endswith(".csv"):
        antiberty_h_dataframes[custom_name[:-4]] = pd.read_csv(file_path)

# Read CSV files for L-chain
for custom_name in os.listdir(output_folder_L):
    file_path = os.path.join(output_folder_L, custom_name)
    if os.path.isfile(file_path) and file_path.endswith(".csv"):
        antiberty_l_dataframes[custom_name[:-4]] = pd.read_csv(file_path)


# In[ ]:


output_folder_H = "Rosetta filtered files/H-rosetta"
output_folder_L = "Rosetta filtered files/L-rosetta"

# Create dictionaries to store DataFrames
rosetta_h_dataframes = {}
rosetta_l_dataframes = {}

# Read CSV files for H-chain
for custom_name in os.listdir(output_folder_H):  # lists all the files in the folder
    file_path = os.path.join(output_folder_H, custom_name) # os.path.join takes the path and joins them
    if os.path.isfile(file_path) and file_path.endswith(".csv"):
        rosetta_h_dataframes[custom_name[:-4]] = pd.read_csv(file_path)  # saves the custom names without the ".csv"

# Read CSV files for L-chain
for custom_name in os.listdir(output_folder_L):
    file_path = os.path.join(output_folder_L, custom_name)  # in brakets - output folder - is were we save it, custom_name is the file (we iterate)
    if os.path.isfile(file_path) and file_path.endswith(".csv"):  # ensures the file_path points to a file and that the file ends with csv
        rosetta_l_dataframes[custom_name[:-4]] = pd.read_csv(file_path)


# ## Merge the rosetta and antiberty dataframe
# polish the dataframes, rename desired columns

# In[ ]:


antibody_h_chain_scores ={}
antibody_l_chain_scores ={}

for custom_name, df in antiberty_h_dataframes.items():  # iterate trough the dictionary
    if custom_name in rosetta_h_dataframes:  # makes sure the names match

        """
        because the positions in both dataframes do not match (be a string) we need to make them a string
        then we merge the dataframes based on the position (foremost)
        we drop unecessary columns and rename the columns for ease
        laslty we reorder the columns in the desired order (to make it easier to read)
        """
        df["position"] = df["position"].astype(str)
        rosetta_h_dataframe = rosetta_h_dataframes[custom_name].astype({"position": str})
        df = df.merge(rosetta_h_dataframe, on=["position", "WT_AA", "MUT_AA"])     
        df = df.drop(columns = ["res_code", "pdb_numbering","VorC"])
        df = df.rename(columns = {"mean_ddG": "Rosetta score", "imgt_numbering": "IMGT numbering"})        
        df = df.reindex(columns = ["chain", "position", "IMGT numbering", "WT_AA", "MUT_AA", "Scaled Antiberty Score", "Rosetta score"])
        antibody_h_chain_scores[custom_name] = df


for custom_name, df in antiberty_l_dataframes.items():
    if custom_name in rosetta_l_dataframes:
        df["position"] = df["position"].astype(str)
        rosetta_l_dataframe = rosetta_l_dataframes[custom_name].astype({"position": str})
        df = df.merge(rosetta_l_dataframe, on=["position", "WT_AA", "MUT_AA"])        
        df = df.drop(columns = ["res_code", "pdb_numbering","VorC"])
        df = df.rename(columns = {"mean_ddG": "Rosetta score", "imgt_numbering": "IMGT numbering"})
        df = df.reindex(columns = ["chain", "position","IMGT numbering", "WT_AA", "MUT_AA", "Scaled Antiberty Score", "Rosetta score"])
        antibody_l_chain_scores[custom_name] = df

# for custom_name, df in antibody_h_chain_scores.items():
#     print(f"Custom Name: {custom_name}")
#     display(df)
# for custom_name, df in antibody_l_chain_scores.items():
#     print(f"Custom Name: {custom_name}")
#     display(df)


# # Merge rosetta and antiberty H- and L-chain dataframes
# calculate the arcsin (to normalise data)
# 
# ## Remove empty dataframes
# After I performed the correlation calculations, I found that there are some empty dataframes. I checked to make sure there were no mistakes, and it appears that the raw files (original -Rosetta) have always been empty. Since this is not a mistake that happened during the table merging, I will remove them.

# In[ ]:


clean_dataframes = {}

# for custom_name in set(antibody_l_chain_scores.keys()) & set(antibody_h_chain_scores.keys()):
#     concatenated_df = pd.concat([antibody_l_chain_scores[custom_name], antibody_h_chain_scores[custom_name]])
#     clean_dataframes[custom_name] = concatenated_df

for custom_name, df_l in antibody_l_chain_scores.items():
    if custom_name in antibody_h_chain_scores:
        df_h = antibody_h_chain_scores[custom_name]  # saves dataframe under value_H
        concatenated_df = pd.concat([df_l, df_h])
        concatenated_df["IMGT numbering"] = concatenated_df["IMGT numbering"].astype(str)
        concatenated_df["pos"] = concatenated_df["IMGT numbering"].apply(lambda x: ''.join(filter(str.isdigit, x)))
        concatenated_df["insertion_code"] = concatenated_df["IMGT numbering"].apply(lambda x: ''.join(filter(str.isalpha, x)))
        concatenated_df["AntiBERTy normalised"] = np.arcsinh(concatenated_df["Scaled Antiberty Score"])
        concatenated_df["Rosetta normalised"] = np.arcsinh(concatenated_df["Rosetta score"])
        if not concatenated_df.empty:
            clean_dataframes[custom_name] = concatenated_df

# for custom_name, df in clean_dataframes.items():
#     print(f"Custom Name: {custom_name}")
#     display(df)


# In[ ]:


def assign_region(row):
    pos = row['pos']
    if pos <= 15:
        return 'FR1-A'
    elif 16 <= pos <= 26:
        return 'FR1-B'
    elif 27 <= pos <= 38:
        return 'CDR1'
    elif 39 <= pos <= 46:
        return 'FR2-C'
    elif 47 <= pos <= 55:
        return "FR2-C'"
    elif 56 <= pos <= 65:
        return 'CDR2'
    elif 66 <= pos <= 74:
        return 'FR3-C"'
    elif 75 <= pos <= 84:
        return 'FR3-D'
    elif 85 <= pos <= 96:
        return 'FR3-E'
    elif 97 <= pos <= 104:
        return 'FR3-F'
    elif 105 <= pos <= 117:
        return 'CDR3'
    else:
        return 'FR4'

for custom_name, df in clean_dataframes.items():
    if 'pos' in df.columns:
        # Convert 'pos' column to int
        df['pos'] = df['pos'].astype(int)

        # Create a new column "region" and apply the logic
        df['region'] = df.apply(assign_region, axis=1)
    else:
        print(f'Warning: DataFrame {custom_name} does not have a "pos" column.')


# ## save the files as csv files

# In[ ]:


output_folder = "Clean files/"

for custom_name, df in clean_dataframes.items():
    file_path = os.path.join(output_folder, f"{custom_name}.csv")  # joins files to the output folder
    df.to_csv(file_path, index=False)  # index = False removes the index in the csv files

