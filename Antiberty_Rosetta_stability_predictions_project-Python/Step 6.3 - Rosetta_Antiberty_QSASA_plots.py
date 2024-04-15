#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt


# ## Read files which are either core or surface

# In[2]:


core_input_folder = "Surface information/core_residues_dataframes"
surface_input_folder = "Surface information/surface_residues_dataframes"
core_dataframes = {}
surface_dataframes = {}

def read_files(input_folder, dataframes):
    for file in os.listdir(input_folder):
        file_path = os.path.join(input_folder, file)
        # print(f"{file_path}")
        
        if file.endswith(".csv"):
            custom_name = file[:7]
            
            with open(file_path, "r") as f:
                df = pd.read_csv(f)
                dataframes[custom_name] = df


read_files(core_input_folder, core_dataframes)
read_files(surface_input_folder, surface_dataframes)


# In[3]:


base_folder = "Final work/Validation against clinically approved drugs"
targets = {
    "3eoa": "3eoa_BA.csv",
    "1bey": "1bey_HL.csv",
    "1ce1": "1ce1_HL.csv",
    "1n8z": "1n8z_BA.csv",
    "7mn8": "7mn8_DC.csv",
    "1l7i": "1l7i_HL.csv",
    
    # mutant variants
    "2fjg": "2fjg_BA.csv", 
    "5jw3_HL": "5jw3_HL.csv",
    "5jw4_MN": "5jw4_MN.csv",
    "5jw4_QR": "5jw4_QR.csv",
    "5jw4_ST": "5jw4_ST.csv",
    "5jw5_AB": "5jw5_AB.csv",
    "5jw5_HL": "5jw5_HL.csv",
    "5fha_HL": "5fha_HL.csv",
    "5fhc_HL": "5fhc_HL.csv",
    "6ws6_AB": "6ws6_AB.csv",
    "6ws6_CD": "6ws6_CD.csv",
    "6ws6_EF": "6ws6_EF.csv",
    "6xdg_CA": "6xdg_CA.csv",
    "6xdg_BD": "6xdg_BD.csv",
}

for target in targets:
    os.makedirs(f"{base_folder}/{target}/{target}_pictures_and_pdb_files", exist_ok =True)


# ## Create the plots showing the QSASA scores across the AntiBERTy and Rosetta scores
# 

# In[6]:


base_folder = "Final work/Validation against clinically approved drugs"
targets = {
    "3eoa": "3eoa_BA",
    "1bey": "1bey_HL",
    "1ce1": "1ce1_HL",
    "1n8z": "1n8z_BA",
    "7mn8": "7mn8_DC",
    "1l7i": "1l7i_HL",
    
    # mutant variants
    "2fjg": "2fjg_BA", 
    "5jw3_HL": "5jw3_HL",
    "5jw4_MN": "5jw4_MN",
    "5jw4_QR": "5jw4_QR",
    "5jw4_ST": "5jw4_ST",
    "5jw5_AB": "5jw5_AB",
    "5jw5_HL": "5jw5_HL",
    "5fha_HL": "5fha_HL",
    "5fhc_HL": "5fhc_HL",
    "6ws6_AB": "6ws6_AB",
    "6ws6_CD": "6ws6_CD",
    "6ws6_EF": "6ws6_EF",
    "6xdg_CA": "6xdg_CA",
    "6xdg_BD": "6xdg_BD",
}


def QSASA_regression_plots(df, target, solvent_access):
    """
    Plot Antiberty and Rosetta scores against QSASA for core or surface residues.

    Parameters:
    - df: DataFrame containing the data.
    - target: Target PDB identifier.
    - solvent_access: String indicating whether the residues are core or surface.

    The plot will calculate the means for Antiberty, Rosetta, and QSASA for each position in the variable region.
    """
    # df["IMGT numbering"] = df["IMGT numbering"].astype(int)
    # calculate the means for the scores (Antiberty, rosetta and QSASA) after grouping them by the wild-type amino acid, position and chain
    position_means = df.groupby(["WT_AA", "IMGT numbering", "chain"])[["AntiBERTy normalised", "Rosetta normalised", "Q(SASA)"]].mean()
    mean_dataframes = pd.DataFrame(position_means)
    mean_filtered = mean_dataframes[mean_dataframes.index.get_level_values("WT_AA") != "C"]  # Remove C values since these are always 0 (in Rosetta)

    plt.clf()
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))
    sns.set_theme()

    title = "Core" if solvent_access == "core" else "Surface"
    plt.suptitle(f"{title} Residues: Antiberty and Rosetta Scores vs. QSASA\n"
                 f"(PDB: {target})")

    # Plot for AntiBERTy scores
    sns.scatterplot(data=mean_filtered, x="Q(SASA)", y="AntiBERTy normalised", hue="chain", ax=axes[0])
    for point, row in mean_filtered.iterrows():
        # for each row in the dataframe it combines the amino acid type (row.name[0]) and position (row.name[1]) - these will be the labels
        # (row["Q(SASA)"], row["AntiBERTy normalised"]) access the values for the columns (so they can match the points to the labels)
        axes[0].annotate(f"{row.name[0]} {row.name[1]}", (row["Q(SASA)"], row["AntiBERTy normalised"]), textcoords="offset points", fontsize=7, xytext=(0,7), ha='center')
    sns.regplot(data=mean_filtered, x="Q(SASA)", y="AntiBERTy normalised", scatter=False, color='red', line_kws={'label':"Correlation line"}, ax=axes[0])

    axes[0].set_xlabel("Q(SASA)")
    axes[0].set_ylabel("AntiBERTy scores")
    
    # Plot for Rosetta scores
    sns.scatterplot(data=mean_filtered, x="Q(SASA)", y="Rosetta normalised", hue="chain", ax=axes[1])
    for point, row in mean_filtered.iterrows():
        axes[1].annotate(f"{row.name[0]} {row.name[1]}", (row["Q(SASA)"], row["Rosetta normalised"]), textcoords="offset points", fontsize=7, xytext=(0,7), ha='center')
    sns.regplot(data=mean_filtered, x="Q(SASA)", y="Rosetta normalised", scatter=False, color='red', line_kws={'label':"Correlation line"}, ax=axes[1])

    axes[1].set_xlabel("Q(SASA)")
    axes[1].set_ylabel("Rosetta scores")

    # Define the output folder and file path for saving the plot
    output_folder = f"{target}/{target}_pictures_and_pdb_files"
    file_path = os.path.join(base_folder, output_folder, f"{solvent_access}_position_scores_over_QSASA.png")

    plt.savefig(file_path, bbox_inches = "tight")
    plt.close()

def process_targets(targets, solvent_access):
    for target, file_name in targets.items():
        df = eval(f"{solvent_access}_dataframes[file_name]")  # eval - strings the variables together (so they can be accessed easier)
        QSASA_regression_plots(df, target, solvent_access)

process_targets(targets, "core")
process_targets(targets, "surface")

