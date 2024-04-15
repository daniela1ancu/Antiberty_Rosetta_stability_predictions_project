#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import pandas as pd


# In[ ]:


Hrosetta_directory = "Initial data/cleaned_result/h_result"
Hrosetta_files = os.listdir(Hrosetta_directory)  # lists all the files in the dictionary
# display(Hrosetta_files)

def read_h_rosetta_files(files):

    """
    from the list of files from the directory, iterate through each file to
    - extract the custom name e.g. 3gje_BA (BA needs to be there since there may be files like 3gje_HL)
    - read the csv file (convert it to dataframe)

    we check if VorC is there (this should be in here but without it it comes up with an error)
    we filter the dataframe to only have the rows from the V region only (rosetta has both V and C in the same table)

    once again we filter the dataframe to contain only the necessary columns 
    technically we do not need the res_code and pdb_numbering but I included it so that we can compare it (see last function)
    """
    
    Hrosetta_dataframes ={}
    
    for file in files:
        file_path = os.path.join(Hrosetta_directory, file)  # joins file to the directory
        
        if os.path.isfile(file_path):
            custom_name = file[:7]  # custom name = "3gje_BA" not "3gje_BA_likelihood.csv"
            Hrosetta_dataframe = pd.read_csv(file_path)
                        
            if 'VorC' in Hrosetta_dataframe.columns:  # Check if 'VorC' is in the columns
                # Filter rows where 'VorC' is 'V'
                Hrosetta_dataframe_filtered = Hrosetta_dataframe[Hrosetta_dataframe["VorC"] == "V"]

                # Check if other desired columns are present
                desired_columns = ["pos", "mean_ddG", "position", "res_code", "pdb_numbering", "VorC", "imgt_numbering"]
                if all(col in Hrosetta_dataframe_filtered.columns for col in desired_columns):
                    Hrosetta_dataframe_filtered = Hrosetta_dataframe_filtered[desired_columns]
                    Hrosetta_dataframe_filtered.insert(0, "chain", "H")
                    Hrosetta_dataframes[custom_name] = Hrosetta_dataframe_filtered
                    
    return Hrosetta_dataframes


Hrosetta_dataframes = read_h_rosetta_files(Hrosetta_files)


Hrosetta_dataframes_new = {}

for custom_name,df in Hrosetta_dataframes.items():

    """
    iterate through each file in the dictionary to:
    the format of the "pos" column looks something like this: H-E23AY
    
   - WT_AA : extract the third character (E)
   - position : 23A (everything between E and Y)
   - MUT_AA : Y (very last character)
    """
    df["WT_AA"] = df["pos"].str[2]
    df["position"] = df["pos"].str[3:-1]
    df["MUT_AA"] = df["pos"].str[-1]
    df = df.drop(columns = ["pos"])  # drops the column   
    df = df.drop_duplicates(subset=["position", "WT_AA", "MUT_AA"])
    Hrosetta_dataframes_new[custom_name] = df
    

# for custom_name, df in Hrosetta_dataframes_new.items():
#     print(f"Custom Name: {custom_name}")
#     display(df)


# In[ ]:


Lrosetta_directory = "Initial data/cleaned_result/l_result"
Lrosetta_files = os.listdir(Lrosetta_directory)
# display(Lrosetta_files)

def read_l_rosetta_files(files):
    Lrosetta_dataframes ={}
    for file in files:
        file_path = os.path.join(Lrosetta_directory, file)
        
        if os.path.isfile(file_path):
            custom_name = file[:7]
            Lrosetta_dataframe = pd.read_csv(file_path)
            
            # Check if 'VorC' is in the columns
            if 'VorC' in Lrosetta_dataframe.columns:
                # Filter rows where 'VorC' is 'V'
                Lrosetta_dataframe_filtered = Lrosetta_dataframe[Lrosetta_dataframe["VorC"] == "V"]

                # Check if other desired columns are present
                desired_columns = ["pos", "mean_ddG", "position", "res_code", "pdb_numbering", "VorC", "imgt_numbering"]
                if all(col in Lrosetta_dataframe_filtered.columns for col in desired_columns):
                    Lrosetta_dataframe_filtered = Lrosetta_dataframe_filtered[desired_columns]
                    Lrosetta_dataframe_filtered.insert(0, "chain", "L")
                    Lrosetta_dataframes[custom_name] = Lrosetta_dataframe_filtered
                    
    return Lrosetta_dataframes


Lrosetta_dataframes = read_l_rosetta_files(Lrosetta_files)


Lrosetta_dataframes_new = {}

for custom_name,df in Lrosetta_dataframes.items():
    df["WT_AA"] = df["pos"].str[2]
    df["position"] = df["pos"].str[3:-1]
    df["MUT_AA"] = df["pos"].str[-1]
    df = df.drop(columns = ["pos"])
    df = df.drop_duplicates(subset=["position", "WT_AA", "MUT_AA"])
    Lrosetta_dataframes_new[custom_name] = df
    

# for custom_name, df in Lrosetta_dataframes_new.items():
#     print(f"Custom Name: {custom_name}")
#     display(df)


# In[ ]:


output_folder_H = "Rosetta filtered files/H-rosetta"
output_folder_L = "Rosetta filtered files/L-rosetta"

def Rosetta_to_csv(dataframes, output_folder):
    for custom_name, df in dataframes.items():
        file_path = os.path.join(output_folder, f"{custom_name}.csv")
        df.to_csv(file_path, index = False)

Rosetta_to_csv(Hrosetta_dataframes_new, output_folder_H)
Rosetta_to_csv(Lrosetta_dataframes_new, output_folder_L)

