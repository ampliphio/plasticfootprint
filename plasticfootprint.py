# %% [markdown]
# # Import Packages

# %%
#import packages 

import pandas as pd
from statistics import mean
import numpy as np

# %% [markdown]
# # Step 1: Files import

# %%
#import input_file 1/3

df_ClientName_Survey = pd.read_csv('/TestInputFile_Survey.csv', sep=';')
df_ClientName_Survey

# %%
#import input_file 2/3

df_ProxyValues = pd.read_csv('/InputFile_ProxyValues_AllCountries.csv', sep=',')
df_ProxyValues

# %%
print(df_ProxyValues.loc[[501]])

# %%
#import input_file 3/3

df_ClientName_SKUs = pd.read_csv('//TestInputFile_SKU_Input.csv', delimiter=',', names=list(range(23)))
new_header = df_ClientName_SKUs.iloc[0] #grab the first row for the header
df_ClientName_SKUs = df_ClientName_SKUs[1:] #take the data less the header row
df_ClientName_SKUs.columns = new_header #set the header row as the df header
df_ClientName_SKUs


# %% [markdown]
# # Step 2: Updates

# %%
#Replace apostroph in string (Add further columns where apostrophs are used as thousands seperator)

df_ClientName_SKUs['Input_SKU_2aPlasticWasteGeneration'].replace('\’','', regex=True, inplace=True)
df_ClientName_SKUs['Input_SKU_BaselinePeriodProductPurchasesProcurement'].replace('\’','', regex=True, inplace=True)

# %%
# change to float
df_ClientName_SKUs['Input_SKU_2aPlasticWasteGeneration'] = df_ClientName_SKUs['Input_SKU_2aPlasticWasteGeneration'].astype(float)
df_ClientName_SKUs['Input_SKU_BaselinePeriodProductPurchasesProcurement'] = df_ClientName_SKUs['Input_SKU_BaselinePeriodProductPurchasesProcurement'].astype(float)
df_ClientName_SKUs['Input_SKU_PlasticWeightperUnit'] = df_ClientName_SKUs['Input_SKU_PlasticWeightperUnit'].astype(float)


# %%
#check data type of a given column

df_ClientName_SKUs.Input_SKU_2aPlasticWasteGeneration.dtype

# %% [markdown]
# # Quality Assurance Checkpoint (QAC) 1: Check File Upload and Contents

# %%
#Generate Country List to be considered in analysis
salesRegions = df_ClientName_SKUs['Input_SKU_SalesRegion'].unique()
print(salesRegions)
len(salesRegions)


# %% [markdown]
# # Step 3: Calcul. 2a Normalization factor

# %%
#Calculate 2a normalization Factor for Approximations
print(df_ClientName_SKUs['Input_SKU_2aPlasticWasteGeneration'].sum())
PWG = df_ClientName_SKUs['Input_SKU_2aPlasticWasteGeneration'].sum() / (df_ProxyValues.iloc[2]['Input_External_Value'] * 100)
PWG

# %% [markdown]
# # Step 4 Calculate Total PWG for Supply Chain Activity Overview

# %%
#Step 4: Calculate Total PWG for Supply Chain Activity Overview

#1a. Operations & Warehouse
OpeWar = PWG * (df_ProxyValues.iloc[0]['Input_External_Value'] * 100)

#1b. Office & Workplace
OffWor = PWG * (df_ProxyValues.iloc[1]['Input_External_Value'] * 100)

#2a. Product & Primary Packaging
ProPri = df_ClientName_SKUs['Input_SKU_2aPlasticWasteGeneration'].sum()

#2b. Shipping & Distribution
ShiDis = PWG * (df_ProxyValues.iloc[3]['Input_External_Value'] * 100)

#3a. Suppliers
Sup = PWG * (df_ProxyValues.iloc[4]['Input_External_Value'] * 100)

#3b. Retailers & Distributors
RetDis = PWG * (df_ProxyValues.iloc[5]['Input_External_Value'] * 100)

#3c. Primary Microplastics
PriMic = PWG * (df_ProxyValues.iloc[6]['Input_External_Value'] * 100)

#Calculate Total PWG
TotPWG = (OpeWar+OffWor+ProPri+ShiDis+Sup+RetDis+PriMic)


# %% [markdown]
# # Step 5: Prepare Greenhouse Gas Emission Calculator

# %%
print(df_ProxyValues.iloc[84])

# %%
# Step 5: Prepare Greenhouse Gas Emission Calculator
# Calculate Emission Factors
# Currently we have only one value to load to calculate average, might be useful in the future

# PP (CO2e/kg)
EmiFac_PP_Avg = np.mean(df_ProxyValues.iloc[940]['Input_External_Value'])

# MW (CO2e/kg)
EmiFac_MW_Avg = np.mean(df_ProxyValues.iloc[941]['Input_External_Value'])

# W2E (CO2e/kg)
EmiFac_W2E_Avg = np.mean(df_ProxyValues.iloc[942]['Input_External_Value'])

# LDF (CO2e/kg)
EmiFac_LDF_Avg = np.mean(df_ProxyValues.iloc[943]['Input_External_Value'])

# REC (CO2e/kg)
EmiFac_REC_Avg = np.mean(df_ProxyValues.iloc[944]['Input_External_Value'])

# Calculate Waste Management Treatment Rate Averages (Activity Data)

# Total IPP (kg)
# 3a. Suppliers
ActDat_IPP = (df_ClientName_SKUs['Input_SKU_BaselinePeriodProductPurchasesProcurement'] * df_ClientName_SKUs['Input_SKU_PlasticWeightperUnit']).sum() 

# Total PWG (kg)
# Calculated in previous steps, variable: 'TotPwg'

# %% [markdown]
# # Create Dataframe for Activity Data and Emission Data (Excel-Sheet "GHG Emissions"); Zwischenschritt 

# %%
# Create Dataframe for Activity Data

# Get all numbers where MWI is mentioned to get an average_mwi

# Copy df with only needed columns
df_mwi = df_ProxyValues[['Input_External_DataPoint', 'Input_External_Value']]
df_mwi = df_mwi[df_mwi['Input_External_DataPoint'].str.contains('Mismanaged Plastic Weight Index')]
df_mwi = df_mwi[df_mwi['Input_External_DataPoint'].str.contains('|'.join(salesRegions))]


# Calculate average MWI
# Could be calculated in one line
Mwi_Avg = np.mean(df_mwi['Input_External_Value'])

# Calculate average Waste2Energy (W2E) Rate
df_W2E = df_ProxyValues[['Input_External_DataPoint', 'Input_External_Value']]
df_W2E = df_W2E[df_W2E['Input_External_DataPoint'].str.contains('Waste2Energy')]
df_W2E = df_W2E[~df_W2E['Input_External_DataPoint'].str.contains('GHG')]
df_W2E = df_W2E[df_W2E['Input_External_DataPoint'].str.contains('|'.join(salesRegions))]

# Calculate average W2E
W2E_Avg = np.mean(df_W2E['Input_External_Value'])

# Build dataframe for Activity Data
df_ActivityData = pd.DataFrame({'Total IPP': [0, 0, 0, 0, ActDat_IPP, 0, 0],
    'Total PWG': [OpeWar, OffWor, ProPri, ShiDis, Sup, RetDis, PriMic],
    'Total MWI': [OpeWar*Mwi_Avg, OffWor*Mwi_Avg, ProPri*Mwi_Avg, ShiDis*Mwi_Avg, Sup*Mwi_Avg, RetDis*Mwi_Avg, PriMic*Mwi_Avg]})

# Add W2E column
df_ActivityData['Total W2E'] = df_ActivityData['Total PWG'] * W2E_Avg

# Calculate average Landfill (LDF) Rate
df_LDF = df_ProxyValues[['Input_External_DataPoint', 'Input_External_Value']]   
df_LDF = df_LDF[df_LDF['Input_External_DataPoint'].str.contains('Landfill')]
df_LDF = df_LDF[~df_LDF['Input_External_DataPoint'].str.contains('GHG')]
df_LDF = df_LDF[df_LDF['Input_External_DataPoint'].str.contains('|'.join(salesRegions))]


# Calculate average LDF
LDF_Avg = np.mean(df_LDF['Input_External_Value'])

# Add LDF column
df_ActivityData['Total LDF'] = df_ActivityData['Total PWG'] * LDF_Avg


# Calculate average Recycling (REC) Rate
df_REC = df_ProxyValues[['Input_External_DataPoint', 'Input_External_Value']]
df_REC = df_REC[df_REC['Input_External_DataPoint'].str.contains('Recycling')]
df_REC = df_REC[~df_REC['Input_External_DataPoint'].str.contains('GHG')]
df_REC = df_REC[~df_REC['Input_External_DataPoint'].str.contains('Leakage')]
df_REC = df_REC[df_REC['Input_External_DataPoint'].str.contains('|'.join(salesRegions))]

# Calculate average REC
REC_Avg = np.mean(df_REC['Input_External_Value'])

# Add REC column
df_ActivityData['Total REC'] = df_ActivityData['Total PWG'] * REC_Avg

# Calculate GHG Emissions

df_GHG_Emissions = pd.DataFrame()

df_GHG_Emissions['Total IPP'] = df_ActivityData['Total IPP'] * EmiFac_PP_Avg
df_GHG_Emissions['Total MWI'] = df_ActivityData['Total MWI'] * EmiFac_MW_Avg
df_GHG_Emissions['Total W2E'] = df_ActivityData['Total W2E'] * EmiFac_W2E_Avg
df_GHG_Emissions['Total LDF'] = df_ActivityData['Total LDF'] * EmiFac_LDF_Avg
df_GHG_Emissions['Total REC'] = df_ActivityData['Total REC'] * EmiFac_REC_Avg

#Total (in kg CO2e)
Tot_GHG_Emi_KG_CO2e = df_GHG_Emissions['Total IPP'].sum() + df_GHG_Emissions['Total MWI'].sum() + df_GHG_Emissions['Total W2E'].sum() + df_GHG_Emissions['Total LDF'].sum() + df_GHG_Emissions['Total REC'].sum()
 
#Co2e/kg of PWG
CO2E_of_PWG = Tot_GHG_Emi_KG_CO2e / TotPWG 

#Total (in mt CO2e)
Tot_GHG_Emi_MT_CO2e = Tot_GHG_Emi_KG_CO2e / 1000

# %% [markdown]
# # QAC 2: Check considered Sales Regions and GHG Factors

# %%
#Print Random GHG Test Values
print(EmiFac_PP_Avg, EmiFac_MW_Avg, EmiFac_W2E_Avg, ActDat_IPP)

# %%
#Print Overview of All Factors
df_EoLratetest = df_ProxyValues[['Input_External_DataPoint', 'Input_External_Value']]
df_EoLratetest = df_EoLratetest[df_EoLratetest['Input_External_DataPoint'].str.contains('|'.join(salesRegions))]
df_EoLratetest = df_EoLratetest[df_EoLratetest['Input_External_DataPoint'].str.contains('W2E')]


print(df_EoLratetest)

#Check average
print(np.mean(df_EoLratetest['Input_External_Value']))

# %%
#Print Average End of Life Factors
EoLFactors = {'EoL':['MWI', 'W2E', 'LDF', 'REC'],'Avg Rates':[Mwi_Avg, W2E_Avg, LDF_Avg, REC_Avg]}
df_EoLFactors = pd.DataFrame(EoLFactors)
print(df_EoLFactors)

#Check if sum of % is 100%
np.sum(df_EoLFactors['Avg Rates'])

# %%
df_mwi

# %%
df_ActivityData

# %%
df_GHG_Emissions


# %% [markdown]
# # Step 6: Calculate BO Analytics Level

# %%
#Calculate Relative Sales Power by Region
df_salesRegions = df_ClientName_SKUs[['Input_SKU_SalesRegion', 'Input_SKU_2aPlasticWasteGeneration']]
df_salesRegions = df_salesRegions.groupby ('Input_SKU_SalesRegion').agg({sum})

#Calculate total 2a PWG
Total_2a_PWG = np.sum(df_salesRegions['Input_SKU_2aPlasticWasteGeneration'])
df_salesRegions["PWG_Prop"] = df_salesRegions['Input_SKU_2aPlasticWasteGeneration'] / Total_2a_PWG

# Hide 2a PWG
df_PWGProp = df_salesRegions[['PWG_Prop']]

df_PWGProp

# len(df_salesRegions)

# %%
# Reset index
df_PWGProp = df_PWGProp.reset_index()
df_mwi = df_mwi.reset_index()

# Rename string in 'Input_External_DataPoint' column
df_mwi['Input_External_DataPoint'] = df_mwi['Input_External_DataPoint'].str.replace('Mismanaged Plastic Weight Index', '')

# Select only desired columns
df_mwi = df_mwi[['Input_External_DataPoint', 'Input_External_Value']]

# Group by 'Input_External_DataPoint' and calculate the mean value of 'Input_External_Value'
df_mwi = df_mwi.groupby('Input_External_DataPoint').agg({'Input_External_Value': 'mean'})

# Reset the index and column names
df_mwi = df_mwi.reset_index()
df_mwi.columns = ['Input_SKU_SalesRegion', 'Mean_Mismanaged_Plastic_Weight_Index']

# Remove empty spaces from the beginning of the values in the 'Input_SKU_SalesRegion' column
df_mwi['Input_SKU_SalesRegion'] = df_mwi['Input_SKU_SalesRegion'].str.strip()

# %%
# Merge df_mwi and df_PWGProp on 'Input_SKU_SalesRegion'
df3 = df_mwi.merge(df_PWGProp, on='Input_SKU_SalesRegion', how='left')

df3 = df3.rename(columns={('PWG_Prop', ''): 'PWG_Prop'})

df3

# %%
# # Here's how it would work (once it works)
df3['avgLR'] = df3['Mean_Mismanaged_Plastic_Weight_Index'] * df3['PWG_Prop'] #to create a new column with the product of rows A and B
avg_LR = df3['avgLR'].sum()


# %%
#Quick Fix Step to move on << not dynamic atm because the join function above doesn't work
df_SupplyChain_LR = pd.DataFrame({'Leakage Rate': [avg_LR, avg_LR, avg_LR, avg_LR, avg_LR, avg_LR, df_ProxyValues.iloc[71]['Input_External_Value']]})

df_SupplyChain_LR

# %%
# prework for Analytics Level df
# Sub-Step 1: Isolate the relevant rows from Proxy Value df to obtain the right External Values
# Only rows with string of interest
LeakageRates = df_ProxyValues[df_ProxyValues['Input_External_DataPoint'].str.contains('Leakage Rate:')]

# Sub-Step 2: Isolate the two relevant columns
df_LeakageRates = LeakageRates[['Input_External_DataPoint',"Input_External_Value"]]

# Sub-Step 3: Reset index
df_LR = df_LeakageRates.reset_index(drop=True)

# Show intermediate df_LR
df_LR

# %%
# Create a empty dataframe which will become Analytics Level df and fill with the values

df_pre_Analytics_Level = pd.DataFrame()

df_pre_Analytics_Level["Input_External_DataPoint"] = df_LR["Input_External_DataPoint"]
df_pre_Analytics_Level["GHG"] = df_GHG_Emissions.sum(axis=1) * 0.001       # simple calculation from GHG_Emissions df
df_pre_Analytics_Level["preEI"] = df_ClientName_Survey.iat[0,7]  # extract value from relevant cell -> ecosystem impact
df_pre_Analytics_Level["prePL"] = df_SupplyChain_LR["Leakage Rate"]            # see steps 1-3 prework
#df_pre_Analytics_Level["prePL"] = df_LR["Input_External_Value"]            # see steps 1-3 prework
df_pre_Analytics_Level["PW"] = df_ActivityData["Total PWG"] * 0.001        # simple calculation from Activity df
df_pre_Analytics_Level

# %%
# Create columns EI and PL

PL = df_pre_Analytics_Level["prePL"] * df_pre_Analytics_Level["PW"]
EI = df_ClientName_Survey.iat[0,7] * PL

# %%
# Add EI and PL to final dataframe: df_Analytics_Level

df_Analytics_Level = pd.DataFrame()

df_Analytics_Level["Input_External_DataPoint"] = df_LR["Input_External_DataPoint"]
df_Analytics_Level['GHG'] = df_pre_Analytics_Level["GHG"]
df_Analytics_Level['EI'] = EI
df_Analytics_Level['PL'] = PL
df_Analytics_Level['PW'] = df_pre_Analytics_Level["PW"]

df_Analytics_Level

# %% [markdown]
# # Step 7: Calculate BO CountryData

# %%
# Step 8: Calculate BO CountryData

df_CountryData = df_ClientName_SKUs.groupby("Input_SKU_SalesRegion", as_index=False)["Input_SKU_2aPlasticWasteGeneration"].sum()

df_CountryData['Percent'] = (df_CountryData['Input_SKU_2aPlasticWasteGeneration'] / df_CountryData['Input_SKU_2aPlasticWasteGeneration'].sum())

df_CountryData['GHG'] = df_CountryData['Percent'] * Tot_GHG_Emi_MT_CO2e

df_CountryData['PW'] = df_Analytics_Level['PW'].sum() * df_CountryData['Percent']


# %%
# Setup MicMac df
df_MicMac_LR = pd.DataFrame()
df_MicMac_LR = df_ClientName_SKUs.groupby("Input_SKU_SalesRegion", as_index=False)["Input_SKU_2aPlasticWasteGeneration"].sum()


# Below is how I'd like to do it (but it doesn't work): Calculate Leakage Rates differentiating between Microplastics and other Supply Chain Activities
df_MicMac_LR['Micropl_LR'] = df_CountryData['PW'] * df_ProxyValues.iloc[6]['Input_External_Value'] * df_SupplyChain_LR.iloc[6]['Leakage Rate']

# %%
#Merge df df_MicMac_LR with df_mwi to create alignment
df_MicMac_LR = df_MicMac_LR.merge(df_mwi, on='Input_SKU_SalesRegion', how='left')
df_MicMac_LR

# %%
df_MicMac_LR2 = df_MicMac_LR.merge(df_CountryData, on='Input_SKU_SalesRegion', how='left')
df_MicMac_LR2

# %%
#Combine Micro- and Macroplastic leakage

df_MicMac_LR2['Marcropl_LR'] = df_MicMac_LR2['PW']  * (1 - df_ProxyValues.iloc[6]['Input_External_Value'])  * df_MicMac_LR2['Mean_Mismanaged_Plastic_Weight_Index']

df_MicMac_LR2['Sum'] = df_MicMac_LR2['Micropl_LR'] + df_MicMac_LR2['Marcropl_LR']

# Simplify df viz
df_MicMac_LR2.drop(['Input_SKU_2aPlasticWasteGeneration_x', 'Input_SKU_2aPlasticWasteGeneration_y', 'Percent', 'GHG' ], axis=1)

# %%
# Add PL Column to df

df_CountryData['PL'] = df_MicMac_LR2['Sum']

df_CountryData['EI'] = df_CountryData['PL'] * df_ClientName_Survey.iat[0,7]

df_CountryData

# %% [markdown]
# # Step 8: Calculate BO Footprint

# %%
#Take Pivot Table Format to proceed

preFootPrints_df = df_ClientName_SKUs[["Input_SKU_SalesRegion", "Input_SKU_FormDescription", "Input_SKU_MaterialType", "Input_SKU_2aPlasticWasteGeneration"]]
preFootPrints_df

#Print United States for crosscheck results with Pivot tab in Excel
# preFootPrints_df.loc[preFootPrints_df['Input_SKU_SalesRegion'] == 'United States']

# %%
# Create a new pre Footprint table with one extra row per Country
preFootPrints_df_2 = pd.concat([preFootPrints_df,
                (preFootPrints_df.groupby('Input_SKU_SalesRegion', sort=False, as_index=False).tail(1)
                   .assign(Input_SKU_FormDescription='Other', Input_SKU_MaterialType='#7, Other', 
                           Input_SKU_2aPlasticWasteGeneration=0))]
               ).sort_index(ignore_index=True)

preFootPrints_df_2

#Note: Remember that the Other row contains values from the 2a SKU import AND the non-2a plastics used e.g. in operations and supply chain

# %%
# Add relevant rows to calculate Footprint values, *0.001 to downscale from mt to kg

preFootPrints_df_2['PWG_Rest'] = preFootPrints_df_2['Input_SKU_2aPlasticWasteGeneration'] * 0.001
preFootPrints_df_2['PWG_Sum'] = preFootPrints_df_2['PWG_Rest'].groupby(preFootPrints_df_2['Input_SKU_SalesRegion']).transform('sum')

# preFootPrints_df_2

#Print United States for crosscheck results with Pivot tab in Excel
preFootPrints_df_2.loc[preFootPrints_df_2['Input_SKU_SalesRegion'] == 'Uganda']

# %%
#These are the values we'll need (instead of the percentages from before with sales proportions) << my problem is that I cannot merge the dfs
column_names = list(df_mwi.columns.values)
column_names


df_mwi

# %%
# merge with df_CountryData to obtain PW value from df_CountryData (needed for further calculation)

# isolate relevant columns of df_CountryData

slice_df_CountryData = df_CountryData[["Input_SKU_SalesRegion", "PW"]]

# merge with preFOotPrints_df_2
preFootPrints_df_3 = pd.merge(preFootPrints_df_2, slice_df_CountryData, on='Input_SKU_SalesRegion')
preFootPrints_df_3.loc[preFootPrints_df_3['Input_SKU_SalesRegion'] == 'British Virgin Islands']

# %%
#Check if we can give 0 values for empty cells (i.e., PWG_Other should be 0 for Input_SKU_MaterialType other than 'Other') << this would require an if clause (see note further up in the script)

preFootPrints_df_3['PWG_Other'] = preFootPrints_df_3['PW'] - preFootPrints_df_3['PWG_Sum']

#Test Sum of PWG and make sure you don't get any negative values
preFootPrints_df_3.loc[preFootPrints_df_3['Input_SKU_MaterialType'] == '#7, Other']

# Comment & Interpretation of Values: PWG_Rest in output below would have to be the sum of all non-others. If it's Zero, we'll get negative results... Not sure if an if-statement or so would be required?


# %%
# Take Proxy values df and isolate rows with mismanaged plastic weight

df_ProxyValues_slice_MPW = df_ProxyValues[df_ProxyValues['Input_External_DataPoint'].str.match('Mismanaged Plastic Weight Index')]
df_ProxyValues_slice_MPW

# %%
# Isolate Country Names from column

df_ProxyValues_slice_MPW['Input_SKU_SalesRegion'] = df_ProxyValues_slice_MPW['Input_External_DataPoint'].str[32:]

df_ProxyValues_slice_MPW

# %%
# Take a subset of df_ProxyValues_slice_MPW with relevant columns

df_ProxyValues_slice_MPW_Country = df_ProxyValues_slice_MPW.loc[:,['Input_SKU_SalesRegion','Input_External_Value']]
# df_ProxyValues_slice_MPW_Country.loc[df_ProxyValues_slice_MPW_Country['Input_SKU_SalesRegion'] == 'Venezuela, RB'] #Check individual country for QA

df_ProxyValues_slice_MPW_Country


# %%
# merge with preFootPrints_df_3

df_pre_BO_Footprints = pd.merge(preFootPrints_df_3, df_ProxyValues_slice_MPW_Country, on='Input_SKU_SalesRegion')
# df_pre_BO_Footprints

#Calculate PL Rest
df_pre_BO_Footprints["PL_Rest"] = df_pre_BO_Footprints["Input_External_Value"] * df_pre_BO_Footprints["PWG_Rest"]

df_pre_BO_Footprints

# %%
# merge df_CountryData with df_pre_BO_Footprints

df_merged_Footprints_CountryData = pd.merge(df_pre_BO_Footprints, df_CountryData, on='Input_SKU_SalesRegion')
df_merged_Footprints_CountryData

# %%
#Combine PWG Rest + Other for QAC 
# df_merged_Footprints_CountryData['PWG_All'] = (df_merged_Footprints_CountryData['PWG_Other'] + df_merged_Footprints_CountryData['PWG_Rest']) # Here we have to add the if clause that PWG_All does not sum up PWG_Rest and PWG_Other but just takes one

df_merged_Footprints_CountryData['PWG_All'] = np.where(df_merged_Footprints_CountryData['Input_SKU_MaterialType'] == '#7, Other', df_merged_Footprints_CountryData['PWG_Other'], df_merged_Footprints_CountryData['PWG_Rest'])
df_merged_Footprints_CountryData.head(40)

# Add column PL_Other (Which is the PL value only for material type "Other")

df_merged_Footprints_CountryData["PL_Other"] = df_merged_Footprints_CountryData['PL'] - df_merged_Footprints_CountryData['PL_Rest'].groupby(df_merged_Footprints_CountryData['Input_SKU_SalesRegion']).transform('sum')


df_merged_Footprints_CountryData['PL_All'] = np.where(df_merged_Footprints_CountryData['Input_SKU_MaterialType'] == '#7, Other', df_merged_Footprints_CountryData['PL_Other'], df_merged_Footprints_CountryData['PL_Rest'])

df_merged_Footprints_CountryData_sliced = df_merged_Footprints_CountryData.loc[:,['Input_SKU_SalesRegion','Input_SKU_FormDescription', "Input_SKU_MaterialType", "PWG_Rest",
                            "PWG_Other","PL_Rest","PL_Other", "PL_All", "PWG_All"]]
df_merged_Footprints_CountryData_sliced

# %%
# Add EI and GHG columns to df_merged_Footprints_CountryData_sliced

df_merged_Footprints_CountryData_sliced["EI_Rest"] =  df_merged_Footprints_CountryData_sliced["PL_Rest"] * df_ClientName_Survey.iat[0,7]
df_merged_Footprints_CountryData_sliced["EI_Other"] =  df_merged_Footprints_CountryData_sliced["PL_Other"] * df_ClientName_Survey.iat[0,7]
df_merged_Footprints_CountryData_sliced['EI_All'] = np.where(df_merged_Footprints_CountryData_sliced['Input_SKU_MaterialType'] == '#7, Other', df_merged_Footprints_CountryData_sliced['EI_Other'], df_merged_Footprints_CountryData_sliced['EI_Rest'])

# df_merged_Footprints_CountryData_sliced["EI_All"] =  df_merged_Footprints_CountryData_sliced["EI_Rest"] + df_merged_Footprints_CountryData_sliced["EI_Other"]  #Combine Rest + Other for QAC

df_merged_Footprints_CountryData_sliced["GHG_Rest"] =  df_merged_Footprints_CountryData_sliced["PWG_Rest"] * CO2E_of_PWG
df_merged_Footprints_CountryData_sliced["GHG_Other"] =  df_merged_Footprints_CountryData_sliced["PWG_Other"] * CO2E_of_PWG
df_merged_Footprints_CountryData_sliced['GHG_All'] = np.where(df_merged_Footprints_CountryData_sliced['Input_SKU_MaterialType'] == '#7, Other', df_merged_Footprints_CountryData_sliced['GHG_Other'], df_merged_Footprints_CountryData_sliced['GHG_Rest'])

# df_merged_Footprints_CountryData_sliced["GHG_All"] =  df_merged_Footprints_CountryData_sliced["GHG_Rest"] + df_merged_Footprints_CountryData_sliced["GHG_Other"] #Combine Rest + Other for QAC

df_merged_Footprints_CountryData_sliced

# %%
#Check if we can give 0 values for empty cells (i.e., PWG_Other should be 0 for Input_SKU_MaterialType other than 'Other') << this would require an if clause (see note further up in the script)

# Reformatted df of Footprints
df_Footprints = df_merged_Footprints_CountryData_sliced.groupby(['Input_SKU_SalesRegion', 'Input_SKU_FormDescription', 'Input_SKU_MaterialType']).sum()
df_Footprints.head(9)



# %% [markdown]
# # Step 9: Calculate BO Material Type
# 

# %%
# Load
df_BO_Material_Type = df_ClientName_SKUs.groupby(['Input_SKU_MaterialType'])['Input_SKU_2aPlasticWasteGeneration'].agg('sum').reset_index()

df_BO_Material_Type.rename(columns = {'Input_SKU_MaterialType':'Category', 'Input_SKU_2aPlasticWasteGeneration':'PW'}, inplace = True)

df_BO_Material_Type['PW'] = df_BO_Material_Type['PW'] / 1000 # Downsize from mt to kg
df_BO_Material_Type

# Calculate value for Others

Sum_PW_Rest = df_BO_Material_Type['PW'].sum()
Sum_PW_Total = df_Analytics_Level['PW'].sum()

OtherCalc = Sum_PW_Total - Sum_PW_Rest
OtherCalc

# Add additional row for Others

df_other = {'Category': '#7, Other', 'PW': OtherCalc}

df_BO_Material_Type = df_BO_Material_Type.append(df_other, ignore_index = True)

# # # #Calculate Plastic Leakage

df_BO_Material_Type_PL_calc = df_Footprints.groupby(['Input_SKU_MaterialType'])['PL_Rest'].agg('sum').reset_index()

df_BO_Material_Type = pd.merge(df_BO_Material_Type, df_BO_Material_Type_PL_calc, left_on='Category', right_on='Input_SKU_MaterialType')

#Update df to fit

df_BO_Material_Type.rename(columns = {'PL_Rest':'PL'}, inplace = True)

df_BO_Material_Type = df_BO_Material_Type.drop('Input_SKU_MaterialType', axis=1)

# # # #Calculate EI

df_BO_Material_Type['EI'] = df_BO_Material_Type['PL'] * df_ClientName_Survey.iat[0,7] #SurveyFile 3333

# # # #Calculate GHG

df_BO_Material_Type['GHG'] = df_BO_Material_Type['PW'] * CO2E_of_PWG 

# # #Print Final BO

df_BO_Material_Type

# %%
#Calculate sum of all PL values for Others
df_B0_Material_Type_Other = df_Footprints.loc[:, ["PL_Other"]]
df_B0_Material_Type_Other = df_B0_Material_Type_Other.filter (like = '#7, Other', axis = 0)

# %%
#Insert PL & EI values for Others
Sum_PL_Rest = df_B0_Material_Type_Other["PL_Other"].sum()

df_BO_Material_Type.at [7, 'PL'] = Sum_PL_Rest 
df_BO_Material_Type.at [7, 'EI'] = df_ClientName_Survey.iat[0,7] * Sum_PL_Rest
df_BO_Material_Type

# %% [markdown]
# # QAC 3: Check & Compare Results of the different BOs

# %%
#Visualize Summary Results in DataFrame
BO_Summary = {'BO':['BO Analytics', 'BO CountryData', 'BO Footprints', 'BO MaterialType'],
              'GHG':[df_Analytics_Level['GHG'].sum(), df_CountryData['GHG'].sum(), df_Footprints['GHG_All'].sum(), df_BO_Material_Type['GHG'].sum()],
              'EI': [df_Analytics_Level['EI'].sum(),df_CountryData['EI'].sum(), df_Footprints['EI_All'].sum(), df_BO_Material_Type['EI'].sum()], 
              'PL': [df_Analytics_Level['PL'].sum(), df_CountryData['PL'].sum(),df_Footprints['PL_All'].sum(),df_BO_Material_Type['PL'].sum()], 
              'PW': [df_Analytics_Level['PW'].sum(),df_CountryData['PW'].sum(), df_Footprints['PWG_All'].sum(),df_BO_Material_Type['PW'].sum()]}
df_BO_Summary = pd.DataFrame(BO_Summary)

print(df_BO_Summary)

# %% [markdown]
# # Step 12: Calculate BO SKU Analysis & Compliance

# %%
#Take df_ClientName_SKUs as basic df for df_BO_SKUAnalysis and
#Transform Input_SKU_%ofRecycledPlasticContent column to float

df_pre_BO_SKUAnalysis = df_ClientName_SKUs[["Input_SKU_FormDescription", "Input_SKU_%ofRecycledPlasticContent", "Input_SKU_MaterialType", "Input_SKU_2aPlasticWasteGeneration"]]
df_pre_BO_SKUAnalysis['Input_SKU_ofRecycledPlasticContent'] = df_pre_BO_SKUAnalysis['Input_SKU_%ofRecycledPlasticContent'].str.rstrip('%').astype('float') 


#group by Forms and create new columns with mean and sum of respective columns

df_BO_SKUAnalysis = (df_pre_BO_SKUAnalysis.groupby('Input_SKU_FormDescription', as_index=False)
       .agg({'Input_SKU_2aPlasticWasteGeneration':'sum','Input_SKU_ofRecycledPlasticContent':'mean'})
       .rename(columns={'Input_SKU_2aPlasticWasteGeneration':'Input_SKU_2aPlasticWasteGeneration_Sum', 'Input_SKU_ofRecycledPlasticContent':'Input_SKU_ofRecycledPlasticContent_Mean'}))

#Divide Input_SKU_2aPlasticWasteGeneration_Sum by 1k
df_BO_SKUAnalysis['Input_SKU_2aPlasticWasteGeneration_Sum'] = df_BO_SKUAnalysis['Input_SKU_2aPlasticWasteGeneration_Sum']/1000


df_BO_SKUAnalysis

# %% [markdown]
# # Step 13: Calculate BO Compliance

# %%
# Comment on Script: Compliance calcs are not completely finalised. Integration with Carlotta's script planned for future.

# %% [markdown]
# # Step 14: Merge all BO df and Print Output File

# %%
#Prepare df_Analytics_Level

df_Analytics_Level.rename(columns={"GHG": "Analytics_GHG", "EI": "Analytics_EI", "PL": "Analytics_PL", 
                                               "PW": "Analytics_PW"},inplace=True,)
df_Analytics_Level


# %%
#Prepare df_CountryData

df_CountryData.rename(columns={"Input_SKU_SalesRegion": "Country_Input_SKU_SalesRegion", "Input_SKU_2aPlasticWasteGeneration": "Country_Input_SKU_2aPlasticWasteGeneration", "Percent": "Country_Percent", 
            "GHG": "Country_GHG", "PW": "Country_PW", "PL": "Country_PL", "EI": "Country_EI"},inplace=True,)

df_CountryData


# %%
df_Footprints

# %%
#Prepare df_BO_Material_Type

df_BO_Material_Type.rename(columns={"Category": "Material_Category", "PW": "Material_PW", "PL": "Material_PL", 
                "EI": "Material_EI", "GHG": "Material_GHG"},inplace=True,)
df_BO_Material_Type

# %%
#Prepare df_BO_SKUAnalysis

df_BO_SKUAnalysis.rename(columns={"Input_SKU_FormDescription": "SKU_Input_SKU_FormDescription", "Input_SKU_2aPlasticWasteGeneration_Sum": "SKU_Input_SKU_2aPlasticWasteGeneration_Sum", 
                "Input_SKU_ofRecycledPlasticContent_Mean": "SKU_Input_SKU_ofRecycledPlasticContent_Mean",},inplace=True,)
df_BO_SKUAnalysis

# %%
# Merge all BO dataframes together to obtain a single df with relevant results for dashboard

Output_ClientName_File = pd.concat([df_Analytics_Level, df_CountryData, df_Footprints, df_BO_Material_Type, df_BO_SKUAnalysis], axis=1)
Output_ClientName_File

# %%
Output_ClientName_File.to_csv('Output_ClientName_File.csv')


