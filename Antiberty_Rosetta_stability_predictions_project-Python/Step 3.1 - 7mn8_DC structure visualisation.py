#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from Bio import PDB
from Bio.PDB import PDBIO, Selection
import urllib.request
import nglview as nv
from matplotlib import colormaps as cmap


# # Import the pdb file for Trastuzumab
# This is because we do not have it downloaded (we can have downloaded and then read)

# In[2]:


"""
for chain.id we can also use the residue.id to find the residues in the structure
    for model in Trastuzumab:
        for chain in model:
            residues = [residue.id for residue in chain]

this creates a list for all the residues in the chain (they are together)
if you want them seperately you need to write something like this:

all_residue_ids = []
for model in Trastuzumab:
    for chain in model:
        residues = [residue.id for residue in chain]
        all_residues_ids.append(residues)

Selection.unfold_entities(selection, target_level) 
- selection can be both the name of the structure or just called "structure"
    - others: models (may be more structures), chain, residue, atom
- target_level - how detailed we want this to be
    - "structure", "model", "chain", "residue", or "atom".
"""


def download_and_parse_pdb(pdb_code):
    pdb_filename = f"{pdb_code.lower()}.pdb"
    pdb_url = f"https://files.rcsb.org/download/{pdb_code.upper()}.pdb"  # link needed to obtain the pdb file

    try:
        with open(pdb_filename, "r"):
            pass
    except FileNotFoundError:
        urllib.request.urlretrieve(trastuzumab_url, trastuzumab_filename)  # download the file and save it locally

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_code, pdb_filename)

    return structure

Trastuzumab = download_and_parse_pdb("7MN8")
light_chain = Selection.unfold_entities(Trastuzumab, "C")[3]


# ## Display the light chain

# In[3]:


print(nv.color.COLOR_SCHEMES)


# In[3]:


chains = [chain.id for chain in Trastuzumab[0]]  
# print(chains)

# select the chains needed 
light_chain = Selection.unfold_entities(Trastuzumab, "C")[3]  # in this case LV was chain 3
heavy_chain = Selection.unfold_entities(Trastuzumab, "C")[4]  # "C" stands for chain (can also be "R" - residues)
structure_light_chain = nv.show_biopython(light_chain)
structure_heavy_chain = nv.show_biopython(heavy_chain)
structure_light_chain.clear_representations()
structure_light_chain.add_representation(repr_type = "cartoon", colorScheme = "bfactor")

structure_light_chain
# structure_heavy_chain


# # Colour light chain (Variable region) of Trastuzumab
# Replace the B-factor values for correlations calculated in the 78mn_DC analysis file (chain_L_correlation_df.csv file)
# Do this for both the Pearson and Spearman correlations (produce 2 different structures)

# In[4]:


light_chain = Selection.unfold_entities(Trastuzumab, "C")[3]


LC_dataframe = pd.read_csv("Final work/Validation against clinically approved drugs/7mn8/7mn8_correlations-DF/Chains/chain_L_correlations.csv")
LC_dataframe.set_index("Position", inplace=True)
cmap = cmap["viridis"]
residue_colours = {}

for residue in light_chain:
    residue_id = residue.get_id()[1]
    if residue_id in LC_dataframe.index:
        correlation = LC_dataframe.loc[residue_id, "Pearson correlation"]
        colour = cmap(correlation)
        residue_colours[residue] = colour

        for atom in residue:
            atom.set_bfactor(correlation)
    else:
        for atom in residue:
            atom.set_bfactor(-1)

io = PDBIO()
io.set_structure(light_chain)
io.save("Final work/Validation against clinically approved drugs/7mn8/7mn8_pictures_and_pdb_files/coloured_light_chain_pearson.pdb")

view = nv.show_structure_file("Final work/Validation against clinically approved drugs/7mn8/7mn8_pictures_and_pdb_files/coloured_light_chain_pearson.pdb")
view.clear_representations()
view.add_representation(repr_type = "surface", colorScheme = "bfactor")
view


# In[5]:


residue_colours = {}

for residue in light_chain:
    residue_id = residue.get_id()[1]
    if residue_id in LC_dataframe.index:
        correlation = LC_dataframe.loc[residue_id, "Spearman correlation"]
        colour = cmap(correlation)
        residue_colours[residue] = colour

        for atom in residue:
            atom.set_bfactor(correlation)
    else:
        for atom in residue:
            atom.set_bfactor(-1)

io = PDBIO()
io.set_structure(light_chain)
io.save("Final work/Validation against clinically approved drugs/7mn8/7mn8_pictures_and_pdb_files/coloured_light_chain_spearman.pdb")

view = nv.show_structure_file("Final work/Validation against clinically approved drugs/7mn8/7mn8_pictures_and_pdb_files/coloured_light_chain_spearman.pdb")
view.clear_representations()
view.add_representation(repr_type = "surface", colorScheme = "bfactor")
view


# # Colour heavy chain (Variable region) of Trastuzumab
# 
# Replace the B-factor values for correlations calculated in the 78mn_DC analysis file (chain_H_correlation_df.csv file) Do this for both the Pearson and Spearman correlations (produce 2 different structures)

# In[6]:


heavy_chain = Selection.unfold_entities(Trastuzumab, "C")[4]


HC_dataframe = pd.read_csv("Final work/Validation against clinically approved drugs/7mn8/7mn8_correlations-DF/Chains/chain_H_correlations.csv")
HC_dataframe.set_index("Position", inplace=True)

residue_colours_h = {}

for residue in heavy_chain:
    residue_id = residue.get_id()[1]
    if residue_id in HC_dataframe.index:
        correlation = HC_dataframe.loc[residue_id, "Pearson correlation"]
        colour = cmap(correlation)
        residue_colours[residue] = colour

        for atom in residue:
            atom.set_bfactor(correlation)
    else:
        for atom in residue:
            atom.set_bfactor(-1)

io = PDBIO()
io.set_structure(heavy_chain)
io.save("Final work/Validation against clinically approved drugs/7mn8/7mn8_pictures_and_pdb_files/coloured_heavy_chain_pearson.pdb")

view = nv.show_structure_file("Final work/Validation against clinically approved drugs/7mn8/7mn8_pictures_and_pdb_files/coloured_heavy_chain_pearson.pdb")
view.clear_representations()
view.add_representation(repr_type = "surface", colorScheme = "bfactor")
view


# In[7]:


residue_colours = {}

for residue in heavy_chain:
    residue_id = residue.get_id()[1]
    if residue_id in HC_dataframe.index:
        correlation = HC_dataframe.loc[residue_id, "Spearman correlation"]
        colour = cmap(correlation)
        residue_colours[residue] = colour

        for atom in residue:
            atom.set_bfactor(correlation)
    else:
        for atom in residue:
            atom.set_bfactor(-1)

io = PDBIO()
io.set_structure(heavy_chain)
io.save("Final work/Validation against clinically approved drugs/7mn8/7mn8_pictures_and_pdb_files/coloured_heavy_chain_spearman.pdb")

view = nv.show_structure_file("Final work/Validation against clinically approved drugs/7mn8/7mn8_pictures_and_pdb_files/coloured_heavy_chain_spearman.pdb")
view.clear_representations()
view.add_representation(repr_type = "surface", colorScheme = "bfactor")
view

