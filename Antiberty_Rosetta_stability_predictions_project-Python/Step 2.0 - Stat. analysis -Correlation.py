#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import pandas as pd
import numpy as np
import scipy.stats


# In[ ]:


input_folder = "Clean files/"
dataframes = {}

def open_dataframes(input_folder):
    for file in os.listdir(input_folder):
        file_path = os.path.join(input_folder, file)

        if file_path.endswith(".csv"):
            custom_name = file[:-4]
            df = pd.read_csv(file_path)
            dataframes[custom_name] = df

open_dataframes(input_folder)


# ## Calculate the correlation for the entire variable region
# .corr to calculate the Pearsons correlation (default) and then also the Spearman (to see if there are some non-linear trends.

# In[ ]:


positive_corr_df_pearson = {}
positive_corr_df_spearman = {}
negative_corr_df = {}
unknown_dataframes = {}


for custom_name, df in dataframes.items():
    pearson_correlation = df["AntiBERTy normalised"].corr(df["Rosetta normalised"], method = "pearson")
    spearman_correlation = df["AntiBERTy normalised"].corr(df["Rosetta normalised"], method = "spearman")
    
    # print(f"AntiBERTy and Rosetta scores correlations {custom_name}: Pearson:{pearson_correlation}, Spearman:{spearman_correlation}\n")

    if pearson_correlation > 0:
        positive_corr_df_pearson[custom_name] = df
    elif spearman_correlation > 0:
        positive_corr_df_spearman[custom_name] = df
    elif pearson_correlation < 0 and spearman_correlation < 0:
        negative_corr_df[custom_name] = df
    else:
        unknown_dataframes[custom_name] = df


# # Calculate correlation for H and L chain:
# Using the filtered dataframes we can calculate the correlation between the scores and see which ones are the highest, lowest as well as which chain has the stronger correlation.

# In[ ]:


"""
To do this we can create a dictionary where all the correlations, variance and standard deviation
THIS DOES NOT SAVE THE DATAFRAMES, only the values. The dataframes can be accessed like normal by calling "dataframes"
there are four main categories two for each chain - positive and negative 
this is further divided based on which correlation type is positive/negative (both, only spearman or only pearson)
[f"positive_{chain_category.lower()}"] this would result in positive_l or positive_h based on the chain_category
"""
correlation_data = {
    "positive_h": {"pearson": {}, "spearman": {}, "both": {}},
    "negative_h": {"pearson": {}, "spearman": {}, "both": {}},
    "positive_l": {"pearson": {}, "spearman": {}, "both": {}},
    "negative_l": {"pearson": {}, "spearman": {}, "both": {}},
    "unknown": {}
}

def calculate_correlation(chain_value, dictionary, custom_name):
    pearson_correlation = chain_value["AntiBERTy normalised"].corr(chain_value["Rosetta normalised"])
    spearman_correlation = chain_value["AntiBERTy normalised"].corr(chain_value["Rosetta normalised"], method = "spearman")
    variance = np.var(chain_value["AntiBERTy normalised"] - chain_value["Rosetta normalised"])
    std_dev = np.sqrt(variance)

    chain_category = chain_value.iloc[0]["chain"]  # this locates the chain type
    
    if pearson_correlation > 0 and spearman_correlation > 0:
        dictionary[f"positive_{chain_category.lower()}"]["both"][custom_name] = (pearson_correlation, spearman_correlation, variance, std_dev)
    elif pearson_correlation < 0 and spearman_correlation < 0:
        dictionary[f"negative_{chain_category.lower()}"]["both"][custom_name] = (pearson_correlation, spearman_correlation, variance, std_dev)
    elif pearson_correlation < 0 and spearman_correlation > 0:
        dictionary[f"positive_{chain_category.lower()}"]["spearman"][custom_name] = (pearson_correlation, spearman_correlation, variance, std_dev)
    elif pearson_correlation > 0 and spearman_correlation < 0:
        dictionary[f"positive_{chain_category.lower()}"]["pearson"][custom_name] = (pearson_correlation, spearman_correlation, variance, std_dev)
    else:
        dictionary["unknown"][custom_name] = (pearson_correlation, spearman_correlation, variance, std_dev)
    
    return pearson_correlation, spearman_correlation, variance, std_dev
    
for custom_name, df in dataframes.items():
    chain_data_h = df[df["chain"] == "H"]  # this would be the chain_value
    chain_data_l = df[df["chain"] == "L"]
    
    pearson_correlation_h, spearman_correlation_h, variability_h, std_dev_h = calculate_correlation(
        chain_data_h, correlation_data, custom_name
    )
    pearson_correlation_l, spearman_correlation_l, variability_l, std_dev_l = calculate_correlation(
        chain_data_l, correlation_data, custom_name
    )


# In[ ]:


print(len(correlation_data["positive_h"]["pearson"]))
print(len(correlation_data["positive_h"]["spearman"]))
print(len(correlation_data["positive_h"]["both"]))
print(len(correlation_data["negative_h"]["both"]))

print(len(correlation_data["positive_l"]["pearson"]))
print(len(correlation_data["positive_l"]["spearman"]))
print(len(correlation_data["positive_l"]["both"]))
print(len(correlation_data["negative_l"]["both"]))


# # Calculate correlations for CDRs only
# Take the amino acids (for each file) and calculate the correlations for this region

# In[ ]:


CDRs_correlation_data = {
    "positive_h": {"pearson": {}, "spearman": {}, "both": {}},
    "negative_h": {"pearson": {}, "spearman": {}, "both": {}},
    "positive_l": {"pearson": {}, "spearman": {}, "both": {}},
    "negative_l": {"pearson": {}, "spearman": {}, "both": {}},
}

def calculate_correlation_regions(chain_value, dictionary, custom_name):
    region_value = chain_value[chain_value["region"].isin(["CDR1", "CDR2", "CDR3"])]
    pearson_correlation = region_value["AntiBERTy normalised"].corr(region_value["Rosetta normalised"])
    spearman_correlation = region_value["AntiBERTy normalised"].corr(region_value["Rosetta normalised"], method="spearman")
    variance = np.var(region_value["AntiBERTy normalised"] - region_value["Rosetta normalised"])
    std_dev = np.sqrt(variance)
    
    chain_category = region_value.iloc[0]["chain"]

    if pearson_correlation > 0 and spearman_correlation > 0:
        dictionary[f"positive_{chain_category.lower()}"]["both"][custom_name] = (pearson_correlation, spearman_correlation, variance, std_dev)
    elif pearson_correlation < 0 and spearman_correlation < 0:
        dictionary[f"negative_{chain_category.lower()}"]["both"][custom_name] = (pearson_correlation, spearman_correlation, variance, std_dev)
    elif pearson_correlation < 0 and spearman_correlation > 0:
        dictionary[f"positive_{chain_category.lower()}"]["spearman"][custom_name] = (pearson_correlation, spearman_correlation, variance, std_dev)
    elif pearson_correlation > 0 and spearman_correlation < 0:
        dictionary[f"positive_{chain_category.lower()}"]["pearson"][custom_name] = (pearson_correlation, spearman_correlation, variance, std_dev)
    else:
        dictionary["unknown_dataframes"][custom_name] = (pearson_correlation, spearman_correlation, variance, std_dev)
    
    return pearson_correlation, spearman_correlation, variance, std_dev

for custom_name, df in dataframes.items():
    chain_data_h = df[df["chain"] == "H"]
    chain_data_l = df[df["chain"] == "L"]

    pearson_correlation_h, spearman_correlation_h, variability_h, std_dev_h = calculate_correlation_regions(
        chain_data_h, CDRs_correlation_data, custom_name
    )
    pearson_correlation_l, spearman_correlation_l, variability_l, std_dev_l = calculate_correlation_regions(
        chain_data_l, CDRs_correlation_data, custom_name
    )


# In[ ]:


print(len(CDRs_correlation_data["positive_h"]["pearson"]))
print(len(CDRs_correlation_data["positive_h"]["spearman"]))
print(len(CDRs_correlation_data["positive_h"]["both"]))
print(len(CDRs_correlation_data["negative_h"]["both"]))

print(len(CDRs_correlation_data["positive_l"]["pearson"]))
print(len(CDRs_correlation_data["positive_l"]["spearman"]))
print(len(CDRs_correlation_data["positive_l"]["both"]))
print(len(CDRs_correlation_data["negative_l"]["both"]))


# # Calculate correlations for each region
# Take the amino acids (for each file) and calculate the correlations for every region in the variable region

# In[ ]:


result_dict = {}

regions = ["FR1-A", "FR1-B", "CDR1", "FR2-C", "FR2-C'", "CDR2", 'FR3-C"', "FR3-D", "FR3-E", "FR3-F", "CDR3", "FR4"]

for region in regions:
    result_dict[region] = {
        "positive_corr_df": {"H": {"both": {}, "pearson": {}, "spearman": {}}, "L": {"both": {}, "pearson": {}, "spearman": {}}},
        "negative_corr_df": {"H": {"both": {}, "pearson": {}, "spearman": {}}, "L": {"both": {}, "pearson": {}, "spearman": {}}},
    }

unknown_dataframes = {}

min_corr_H = float('inf')
max_corr_H = float('-inf')

min_corr_L = float('inf')
max_corr_L = float('-inf')

def calculate_correlation_regions(df, chain, result_dict, custom_name):
    global min_corr_H, max_corr_H, min_corr_L, max_corr_L
    # print(f"Processing {custom_name} for region {region} and chain {chain}")
    chain_value = df[df["chain"] == chain]
    region_value = chain_value[chain_value["region"] == region]
    pearson_correlation = region_value["AntiBERTy normalised"].corr(region_value["Rosetta normalised"])
    spearman_correlation = region_value["AntiBERTy normalised"].corr(region_value["Rosetta normalised"], method = "spearman")
    variance = np.var(region_value["AntiBERTy normalised"] - region_value["Rosetta normalised"])
    std_dev = np.sqrt(variance)

    if chain == "H":
        if spearman_correlation > max_corr_H:
            max_corr_H = spearman_correlation
        elif spearman_correlation < min_corr_H:
            min_corr_H = spearman_correlation
    elif chain == "L":
        if spearman_correlation > max_corr_L:
            max_corr_L = spearman_correlation
        elif spearman_correlation < min_corr_L:
            min_corr_L = spearman_correlation

    
    if pearson_correlation > 0 and spearman_correlation > 0:
        result_dict[region]["positive_corr_df"][chain]["both"][custom_name] = region_value
    elif pearson_correlation < 0 and spearman_correlation < 0:
        result_dict[region]["negative_corr_df"][chain]["both"][custom_name] = region_value
    elif pearson_correlation < 0 and spearman_correlation > 0:
        result_dict[region]["positive_corr_df"][chain]["spearman"][custom_name] = region_value
    elif pearson_correlation > 0 and spearman_correlation < 0:
        result_dict[region]["positive_corr_df"][chain]["pearson"][custom_name] = region_value
    else:
        unknown_dataframes[custom_name] = region_value
    
    return pearson_correlation, spearman_correlation, variance, std_dev, min_corr_H, max_corr_H, min_corr_L, max_corr_L

for custom_name, df in dataframes.items():
    for region in regions:
        for chain in ["H", "L"]:
            calculate_correlation_regions(df, chain, result_dict, custom_name)
print(f"Overall Minimum Spearman Correlation for Chain H: {min_corr_H}")
print(f"Overall Maximum Spearman Correlation for Chain H: {max_corr_H}")
print(f"Overall Minimum Spearman Correlation for Chain L: {min_corr_L}")
print(f"Overall Maximum Spearman Correlation for Chain L: {max_corr_L}")


# ## Print length of lists 
# Printing counts of positive and negative correlations (structures) for each region and antibody types (H and L) using various methods (both, Pearson, Spearman)

# In[ ]:


regions = ["FR1-A", "FR1-B", "CDR1", "FR2-C", "FR2-C'", "CDR2", 'FR3-C"', "FR3-D", "FR3-E", "FR3-F", "CDR3", "FR4"]

for region in regions:
    positive_corr_H_both = len(result_dict[region]["positive_corr_df"]["H"]["both"])
    positive_corr_L_both = len(result_dict[region]["positive_corr_df"]["L"]["both"])
    negative_corr_H_both = len(result_dict[region]["negative_corr_df"]["H"]["both"])
    negative_corr_L_both = len(result_dict[region]["negative_corr_df"]["L"]["both"])

    positive_corr_H_pearson = len(result_dict[region]["positive_corr_df"]["H"]["pearson"])
    positive_corr_L_pearson = len(result_dict[region]["positive_corr_df"]["L"]["pearson"])
    negative_corr_H_pearson = len(result_dict[region]["negative_corr_df"]["H"]["pearson"])
    negative_corr_L_pearson = len(result_dict[region]["negative_corr_df"]["L"]["pearson"])

    positive_corr_H_spearman = len(result_dict[region]["positive_corr_df"]["H"]["spearman"])
    positive_corr_L_spearman= len(result_dict[region]["positive_corr_df"]["L"]["spearman"])
    negative_corr_H_spearman = len(result_dict[region]["negative_corr_df"]["H"]["spearman"])
    negative_corr_L_spearman = len(result_dict[region]["negative_corr_df"]["L"]["spearman"])

    print(f"Region: {region}")
    print(f"positive_corr_H_both: {positive_corr_H_both}")
    print(f"negative_corr_H_both: {negative_corr_H_both}")
    print(f"positive_corr_H_pearson: {positive_corr_H_pearson}")
    print(f"negative_corr_H_pearson: {negative_corr_H_pearson}")
    print(f"positive_corr_H_spearman: {positive_corr_H_spearman}")
    print(f"negative_corr_H_spearman: {negative_corr_H_spearman}")
    print("")
    
    print(f"positive_corr_L_both: {positive_corr_L_both}")
    print(f"positive_corr_L_pearson: {positive_corr_L_pearson}")
    print(f"positive_corr_L_spearman: {positive_corr_L_spearman}")
    print(f"negative_corr_L_both: {negative_corr_L_both}")
    print(f"negative_corr_L_pearson: {negative_corr_L_pearson}")
    print(f"negative_corr_L_spearman: {negative_corr_L_spearman}")
    print("")
    print("")


# ## In certain cases, there are structures which do not have certain regions (e.g. CDR2) this is something to look into

# In[ ]:


regions = ["FR1-A", "FR1-B", "CDR1", "FR2-C", "FR2-C'", "CDR2", 'FR3-C"', "FR3-D", "FR3-E", "FR3-F", "CDR3", "FR4"]


for custom_name, df in unknown_dataframes.items():
    for region in regions:
        chain_value = df[df["chain"] == chain]
        region_value = chain_value[chain_value["region"] == region]
        pearson_correlation = region_value["AntiBERTy normalised"].corr(region_value["Rosetta normalised"])
        spearman_correlation = region_value["AntiBERTy normalised"].corr(region_value["Rosetta normalised"], method="spearman")
        variability = np.var(region_value["AntiBERTy normalised"] - region_value["Rosetta normalised"])
        std_dev = np.sqrt(variability)
        # print(custom_name, region, chain,  pearson_correlation, spearman_correlation, variability, std_dev)


# ## Regional data - dataframe
# Taking the correlations calculated back for each region of the antibody, we take this and create a dataframe which will contain the type of correlation  (positive vs negative), by which method (only pearson, only spearman or both), as well as the region and which antibody falls in this category. 

# In[ ]:


data_region_list = []
positive = "positive_corr_df"
negative = "negative_corr_df"

def summary_dataframe(result_dict):
    for region, region_dict in result_dict.items():
        for correlation_type, correlation_type_dict in region_dict.items():
            if correlation_type == positive:
                for chain, chain_dict in correlation_type_dict.items():
                    for correlation, df_dict in chain_dict.items():
                        for antibody_iden, df in df_dict.items():
                            data_region_list.append({
                                "Region": region,
                                "Chain": chain,
                                "Correlation type": "Positive",
                                "Method": correlation,
                                "Antibody IDEN code": antibody_iden,
                            })
            elif correlation_type == negative:
                for chain, chain_dict in correlation_type_dict.items():
                    for correlation, df_dict in chain_dict.items():
                        for antibody_iden, df in df_dict.items():
                            data_region_list.append({
                                "Region": region,
                                "Chain": chain,
                                "Correlation type": "Negative",
                                "Method": correlation,
                                "Antibody IDEN code": antibody_iden,
                            })


summary_dataframe(result_dict)

regional_data = pd.DataFrame(data_region_list)
# regional_data.to_csv("test.csv", index = False)

display(regional_data)


# ## Split regional Dataframe 
# we can then take this dataframe and split it by the region (e.g. CDR2)

# In[ ]:


directory_path = "Region Correlations DF"

if not os.path.exists(directory_path):
    os.makedirs(directory_path)

else:
    pass

region_dataframes = {}
dataframes_path = "Region Correlations DF/"

for region, group in regional_data.groupby("Region"):
    region_dataframes[region + '_dataframe'] = group
    file_path = os.path.join(dataframes_path, f"{region}_dataframe.csv")
    group.to_csv(file_path, index = False)

display(region_dataframes["CDR2_dataframe"])

