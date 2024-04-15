#!/usr/bin/env python
# coding: utf-8

# In[3]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# ## Open the 7mn8_DC file
# we are working with Trastuzumab (7mn8)

# In[4]:


Trastuzumab = pd.read_csv("Clean files/7mn8_DC.csv")


# ## Heatmaps

# In[5]:


"""
First we look at the entire chain (this would reflect the matrices produced by AntiBERTy)
IMGT numbering scheme is used
"""

L_chain = Trastuzumab[Trastuzumab["chain"] == "L"]
H_chain = Trastuzumab[Trastuzumab["chain"] == "H"]

def get_heatmaps(df, chain_type):
    heatmap_data = df.pivot_table(index=['MUT_AA'], columns=['IMGT numbering', 'WT_AA'], values='AntiBERTy normalised')
    plt.clf()
    
    plt.figure(figsize=(30, 7))
    sns.heatmap(heatmap_data, cmap='Oranges', linewidths=0.25)
    plt.title(f'Mutational Impact Heatmap (AntiBERTy) of {chain_type} chain variable region')
    plt.xlabel('Position')
    plt.ylabel('Mutant Amino Acid', size=25)
    plt.savefig(f"Final work/Validation against clinically approved drugs/7mn8/7mn8_pictures_and_pdb_files/{chain_type}V_mutational_profile_heatmap_antiberty.png", bbox_inches="tight")
    plt.show()

get_heatmaps(L_chain, "L")
get_heatmaps(H_chain, "H")


# In[29]:


"""
First we look at the entire chain (this would reflect the matrices produced by Rosetta)
IMGT numbering scheme is used
"""

L_chain = Trastuzumab[Trastuzumab["chain"] == "L"]
H_chain = Trastuzumab[Trastuzumab["chain"] == "H"]

def get_heatmaps(df, chain_type):
    heatmap_data = df.pivot_table(index=['MUT_AA'], columns=['IMGT numbering', 'WT_AA'], values='Rosetta normalised')
    plt.clf()
    
    plt.figure(figsize=(30, 7))
    sns.heatmap(heatmap_data, cmap='Oranges', linewidths=0.25)
    #  cbar_kws={'orientation': 'horizontal', 'shrink': 0.30} - for the colourmap
    plt.title(f'Mutational Impact Heatmap (Rosetta) of {chain_type} chain variable region')
    plt.xlabel('Position')
    plt.ylabel('Mutant Amino Acid', size=25)
    plt.savefig(f"Final work/Validation against clinically approved drugs/7mn8/7mn8_pictures_and_pdb_files/{chain_type}V_mutational_profile_heatmap_rosetta.png", bbox_inches="tight")
    plt.show()

get_heatmaps(L_chain, "L")
get_heatmaps(H_chain, "H")


# In[5]:


"""
Next we check the first two framework regions.
In these regions we noticed that there is a high level of disagreement and we wanted to check why that is.
PDB numbering scheme is used
"""

L_chain = Trastuzumab[Trastuzumab["chain"] == "L"]
L_chain_FR1 = L_chain[L_chain["region"].isin(["FR1-A", "FR1-B"])]
H_chain = Trastuzumab[Trastuzumab["chain"] == "H"]
H_chain_FR1 = H_chain[H_chain["region"].isin(["FR1-A", "FR1-B"])]


def get_heatmaps(df, chain_type):
    heatmap_data = df.pivot_table(index=['MUT_AA'], columns=['position', 'WT_AA'], values='AntiBERTy normalised')
    plt.clf()
    
    plt.figure(figsize=(15, 7))
    sns.heatmap(heatmap_data, cmap='Oranges', linewidths=0.5)
    plt.title(f'Mutational Impact Heatmap (AntiBERTy) of {chain_type} chain variable region')
    plt.xlabel('Position')
    plt.ylabel('MUT_AA')
    plt.savefig(f"Final work/Validation against clinically approved drugs/Trastuzumab/Trastuzumab pictures and pdb files/{chain_type}V_FR1_mutational_profile_heatmap_antiberty.png", bbox_inches="tight")
    plt.show()

get_heatmaps(L_chain_FR1, "L")
get_heatmaps(H_chain_FR1, "H")


# In[4]:


"""
This are the Rosetta heatmaps.
In these regions we noticed that there is a high level of disagreement and we wanted to check why that is.
PDB numbering scheme is used
"""

L_chain = Trastuzumab[Trastuzumab["chain"] == "L"]
L_chain_FR1 = L_chain[L_chain["region"].isin(["FR1-A", "FR1-B"])]
H_chain = Trastuzumab[Trastuzumab["chain"] == "H"]
H_chain_FR1 = H_chain[H_chain["region"].isin(["FR1-A", "FR1-B"])]


def get_heatmaps(df, chain_type):
    heatmap_data = df.pivot_table(index=['MUT_AA'], columns=['position', 'WT_AA'], values='Rosetta normalised')
    plt.clf()
    
    plt.figure(figsize=(15, 7))
    sns.heatmap(heatmap_data, cmap='bwr', linewidths=0.5)
    plt.title(f'Mutational Impact Heatmap (Rosetta) of {chain_type} chain variable region')
    plt.xlabel('Position')
    plt.ylabel('MUT_AA')
    plt.savefig(f"Final work/Validation against clinically approved drugs/Trastuzumab/Trastuzumab pictures and pdb files/{chain_type}V_FR1_mutational_profile_heatmap_rosetta.png", bbox_inches="tight")
    plt.show()


get_heatmaps(L_chain_FR1, "L")
get_heatmaps(H_chain_FR1, "H")


# ## Confirm results across other strctures 
# We observed that Rosetta marks certain residues seemed to mark certain residues (mainly proline, some glycine residues and some bulkly residues mutations) as destabilising. so we wanted to confirm this against other structures besides 7mn8.   

# In[10]:


Trastuzumab_3eoa = pd.read_csv("Clean files/3eoa_BA.csv")

L_chain = Trastuzumab_3eoa[Trastuzumab_3eoa["chain"] == "L"]
H_chain = Trastuzumab_3eoa[Trastuzumab_3eoa["chain"] == "H"]

def get_heatmaps(df, chain_type):
    heatmap_data = df.pivot_table(index=['MUT_AA'], columns=['IMGT numbering', 'WT_AA'], values='AntiBERTy normalised')
    plt.clf()
    
    plt.figure(figsize=(30, 7))
    sns.heatmap(heatmap_data, cmap='Oranges', linewidths=0.25)
    plt.title(f'Mutational Impact Heatmap (AntiBERTy) of {chain_type} chain variable region')
    plt.xlabel('Position')
    plt.ylabel('Mutant Amino Acid', size=25)
    plt.savefig(f"Final work/Validation against clinically approved drugs/3eoa/3eoa_pictures_and_pdb_files/{chain_type}V_mutational_profile_heatmap_antiberty.png", bbox_inches="tight")
    plt.show()

get_heatmaps(L_chain, "L")
get_heatmaps(H_chain, "H")


# In[11]:


L_chain = Trastuzumab_3eoa[Trastuzumab_3eoa["chain"] == "L"]
H_chain = Trastuzumab_3eoa[Trastuzumab_3eoa["chain"] == "H"]

def get_heatmaps(df, chain_type):
    heatmap_data = df.pivot_table(index=['MUT_AA'], columns=['IMGT numbering', 'WT_AA'], values='Rosetta normalised')
    plt.clf()
    
    plt.figure(figsize=(30, 7))
    sns.heatmap(heatmap_data, cmap='Oranges', linewidths=0.25)
    plt.title(f'Mutational Impact Heatmap (AntiBERTy) of {chain_type} chain variable region')
    plt.xlabel('Position')
    plt.ylabel('Mutant Amino Acid', size=25)
    plt.savefig(f"Final work/Validation against clinically approved drugs/3eoa/3eoa_pictures_and_pdb_files/{chain_type}V_mutational_profile_heatmap_rosetta.png", bbox_inches="tight")
    plt.show()

get_heatmaps(L_chain, "L")
get_heatmaps(H_chain, "H")


# In[18]:


Trastuzumab_1ce1 = pd.read_csv("Clean files/1ce1_HL.csv")

L_chain = Trastuzumab_1ce1[Trastuzumab_1ce1["chain"] == "L"]
H_chain = Trastuzumab_1ce1[Trastuzumab_1ce1["chain"] == "H"]

def get_heatmaps(df, chain_type):
    heatmap_data = df.pivot_table(index=['MUT_AA'], columns=['IMGT numbering', 'WT_AA'], values='AntiBERTy normalised')
    plt.clf()
    
    plt.figure(figsize=(30, 7))
    sns.heatmap(heatmap_data, cmap='Oranges', linewidths=0.25)
    plt.title(f'Mutational Impact Heatmap (AntiBERTy) of {chain_type} chain variable region')
    plt.xlabel('Position')
    plt.ylabel('Mutant Amino Acid', size=25)
    plt.savefig(f"Final work/Validation against clinically approved drugs/1ce1/1ce1_pictures_and_pdb_files/{chain_type}V_mutational_profile_heatmap_antiberty.png", bbox_inches="tight")
    plt.show()

get_heatmaps(L_chain, "L")
get_heatmaps(H_chain, "H")


# In[19]:


L_chain = Trastuzumab_1ce1[Trastuzumab_1ce1["chain"] == "L"]
H_chain = Trastuzumab_1ce1[Trastuzumab_1ce1["chain"] == "H"]

def get_heatmaps(df, chain_type):
    heatmap_data = df.pivot_table(index=['MUT_AA'], columns=['IMGT numbering', 'WT_AA'], values='Rosetta normalised')
    plt.clf()
    
    plt.figure(figsize=(30, 7))
    sns.heatmap(heatmap_data, cmap='Oranges', linewidths=0.25)
    plt.title(f'Mutational Impact Heatmap (AntiBERTy) of {chain_type} chain variable region')
    plt.xlabel('Position')
    plt.ylabel('Mutant Amino Acid', size=25)
    plt.savefig(f"Final work/Validation against clinically approved drugs/1ce1/1ce1_pictures_and_pdb_files/{chain_type}V_mutational_profile_heatmap_rosetta.png", bbox_inches="tight")
    plt.show()

get_heatmaps(L_chain, "L")
get_heatmaps(H_chain, "H")

