# -*- coding: utf-8 -*-
"""
Author: Madeleine NICOLAS       
        madeleine.nicolas@univ-rennes1.fr

Purpose: Calculates PG, Kc and RET for each grid cell for HydMod

Requirements: PET, land use properties per grid cell

Remarks: 

Date:   April 2018
"""

import os
import pandas as pd
import datetime
import itertools
import string

AlphaIndex = list(itertools.chain(string.ascii_uppercase,(''.join(pair) for pair in itertools.product(string.ascii_uppercase, repeat=2))))

os.chdir(r"C:\Users\Madeleine\Desktop\MAHESHWARAM\01 - SOIL MOISTURE MODEL\02 - INPUT DATA") # Sets working directory

########## INPUT DATA ##########

InputClimate = pd.read_csv("INPUT_climate_1995-2015.txt",sep='\t',header=0,index_col=0) # Potential evaporation & rainfall (in mm)
InputCategories = pd.read_csv("INPUT_categories_gridded_3.txt",sep='\t',header=0,index_col=0) # SoilClass and LU for each grid cell
InputCategoriesUnique = InputCategories.drop_duplicates() # SoilClass and LU combinations are considered once each
InputCategoriesUnique.index = AlphaIndex[0:len(InputCategoriesUnique)] # Sets indexes as letters

LU_categories = list(set(InputCategories['LandUse']))
SoilType_categories = list(set(InputCategories['SoilClass']))

InputKc = pd.DataFrame(columns=InputCategoriesUnique.index,index = pd.to_datetime(InputClimate.index))
InputCf = pd.DataFrame(columns=InputCategoriesUnique.index,index = pd.to_datetime(InputClimate.index))
InputPG = pd.DataFrame(columns=InputCategoriesUnique.index,index = pd.to_datetime(InputClimate.index))
InputRET = pd.DataFrame(columns=InputCategoriesUnique.index,index = pd.to_datetime(InputClimate.index)) 

# CROP COEFFICIENT

InputKc[[i for i, x in enumerate(InputCategoriesUnique['LandUse'] == 'Scrub') if x]] = 0.7
InputKc[[i for i, x in enumerate(InputCategoriesUnique['LandUse'] == 'Forest') if x]] = 1.2      
InputKc[[i for i, x in enumerate(InputCategoriesUnique['LandUse'] == 'Fruit') if x]] = 0.8  
InputKc[[i for i, x in enumerate(InputCategoriesUnique['LandUse'] == 'Tank') if x]] = 1  
InputKc[[i for i, x in enumerate(InputCategoriesUnique['LandUse'] == 'Urban') if x]] = 1  
    
for yr in range(min(InputKc.index.year),max(InputKc.index.year)+1):

    InputKc.at[(InputKc.index.date>=datetime.date(yr, 6, 01))*(InputKc.index.date<=datetime.date(yr, 6, 30)),InputCategoriesUnique['LandUse'] == 'Vegetable'] = 0.95
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 7, 01))*(InputKc.index.date<=datetime.date(yr, 7, 31)),InputCategoriesUnique['LandUse'] == 'Vegetable'] = 0.7
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 8, 01))*(InputKc.index.date<=datetime.date(yr, 11, 25)),InputCategoriesUnique['LandUse'] == 'Vegetable'] = 1.05
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 11, 26))*(InputKc.index.date<=datetime.date(yr, 12, 31)),InputCategoriesUnique['LandUse'] == 'Vegetable'] = 0.95
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 01, 01))*(InputKc.index.date<=datetime.date(yr, 1, 10)),InputCategoriesUnique['LandUse'] == 'Vegetable'] = 0.95
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 1, 11))*(InputKc.index.date<=datetime.date(yr, 2, 10)),InputCategoriesUnique['LandUse'] == 'Vegetable'] = 0.7       
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 2, 11))*(InputKc.index.date<=datetime.date(yr, 5, 31)),InputCategoriesUnique['LandUse'] == 'Vegetable'] = 0.85

    InputKc.at[(InputKc.index.date>=datetime.date(yr, 6, 01))*(InputKc.index.date<=datetime.date(yr, 6, 30)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 1.05
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 7, 01))*(InputKc.index.date<=datetime.date(yr, 10, 31)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 1.2
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 11, 01))*(InputKc.index.date<=datetime.date(yr, 11, 15)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 1
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 11, 16))*(InputKc.index.date<=datetime.date(yr, 12, 15)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 1.05
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 12, 16))*(InputKc.index.date<=datetime.date(yr, 12, 31)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 1.2
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 1, 1))*(InputKc.index.date<=datetime.date(yr, 3, 25)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 1.2      
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 3, 26))*(InputKc.index.date<=datetime.date(yr, 5, 31)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 1
    
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 6, 01))*(InputKc.index.date<=datetime.date(yr, 6, 30)),InputCategoriesUnique['LandUse'] == 'PaddyK'] = 1.05
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 7, 01))*(InputKc.index.date<=datetime.date(yr, 10, 31)),InputCategoriesUnique['LandUse'] == 'PaddyK'] = 1.2
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 11, 01))*(InputKc.index.date<=datetime.date(yr, 12, 31)),InputCategoriesUnique['LandUse'] == 'PaddyK'] = 1
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 1, 1))*(InputKc.index.date<=datetime.date(yr, 5, 31)),InputCategoriesUnique['LandUse'] == 'PaddyK'] = 1     

# EVAPORATION

InputRET = InputKc.multiply(InputClimate['PET_mm'],axis="index")*0.9 #0.9 is the regional coefficient Kp
InputRET = InputRET/1000 # Inputs must be in m

# PUMPING

InputPG[[i for i, x in enumerate(InputCategoriesUnique['LandUse'] == 'Scrub') if x]] = 0
InputPG[[i for i, x in enumerate(InputCategoriesUnique['LandUse'] == 'Forest') if x]] = 0
InputPG[[i for i, x in enumerate(InputCategoriesUnique['LandUse'] == 'Tank') if x]] = 0
InputPG[[i for i, x in enumerate(InputCategoriesUnique['LandUse'] == 'Fruit') if x]] = 1.9
InputPG[[i for i, x in enumerate(InputCategoriesUnique['LandUse'] == 'Urban') if x]] = 0
InputPG[[i for i, x in enumerate(InputCategoriesUnique['LandUse'] == 'Vegetable') if x]] = 7.7

for yr in range(min(InputPG.index.year),max(InputPG.index.year)+1):
    
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 6, 1))*(InputPG.index.date<=datetime.date(yr, 6, 30)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 10.1*0.1 # Nursery (Kharif)     
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 7, 1))*(InputPG.index.date<=datetime.date(yr, 9, 30)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 10.1 # Irrigation (Kharif) 
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 10, 1))*(InputPG.index.date<=datetime.date(yr, 10, 15)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 0 # Harvesting (Kharif) 
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 10, 16))*(InputPG.index.date<=datetime.date(yr, 12, 15)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 15.2*0.1 # Maintenance + Nursery
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 12, 15))*(InputPG.index.date<=datetime.date(yr, 12, 31)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 15.2 # Irrigation (Rabi)
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 1, 1))*(InputPG.index.date<=datetime.date(yr, 4, 15)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 15.2 # Irrigation (Rabi)
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 4, 16))*(InputPG.index.date<=datetime.date(yr, 4, 30)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 0 # Harvesting (Rabi)
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 5, 1))*(InputPG.index.date<=datetime.date(yr, 5, 31)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 15.2*0.1 # Maintenance + Nursery (Rabi)
                
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 6, 1))*(InputPG.index.date<=datetime.date(yr, 6, 30)),InputCategoriesUnique['LandUse'] == 'PaddyK'] = 10.1*0.1 # Nursery (Kharif)     
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 7, 1))*(InputPG.index.date<=datetime.date(yr, 9, 30)),InputCategoriesUnique['LandUse'] == 'PaddyK'] = 10.1 # Irrigation (Kharif) 
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 10, 1))*(InputPG.index.date<=datetime.date(yr, 12, 31)),InputCategoriesUnique['LandUse'] == 'PaddyK'] = 0 # Harvesting (Kharif) + No pumping
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 1, 1))*(InputPG.index.date<=datetime.date(yr, 5, 31)),InputCategoriesUnique['LandUse'] == 'PaddyK'] = 0 # No pumping

InputPG = InputPG/1000 # Inputs must be in m

# PERCENTAGE OF RECHARGE EQUAL TO RETURN FLOW

InputCf[[i for i, x in enumerate(InputCategoriesUnique['LandUse'] == 'Scrub') if x]] = 1
InputCf[[i for i, x in enumerate(InputCategoriesUnique['LandUse'] == 'Forest') if x]] = 1      
InputCf[[i for i, x in enumerate(InputCategoriesUnique['LandUse'] == 'Fruit') if x]] = 1  
InputCf[[i for i, x in enumerate(InputCategoriesUnique['LandUse'] == 'Tank') if x]] = 1  
InputCf[[i for i, x in enumerate(InputCategoriesUnique['LandUse'] == 'Urban') if x]] = 1
InputCf[[i for i, x in enumerate(InputCategoriesUnique['LandUse'] == 'Vegetable') if x]] = 0.24

for yr in range(min(InputCf.index.year),max(InputCf.index.year)+1):
    
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 6, 1))*(InputCf.index.date<=datetime.date(yr, 6, 30)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 0.51 # Nursery (Kharif)     
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 7, 1))*(InputCf.index.date<=datetime.date(yr, 9, 30)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 0.51 # Irrigation (Kharif) 
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 10, 1))*(InputCf.index.date<=datetime.date(yr, 10, 15)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 1 # Harvesting (Kharif) 
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 10, 16))*(InputCf.index.date<=datetime.date(yr, 12, 15)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 0.48 # Maintenance + Nursery
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 12, 15))*(InputCf.index.date<=datetime.date(yr, 12, 31)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 0.48 # Irrigation (Rabi)
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 1, 1))*(InputCf.index.date<=datetime.date(yr, 4, 15)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 0.48 # Irrigation (Rabi)
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 4, 16))*(InputCf.index.date<=datetime.date(yr, 4, 30)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 1 # Harvesting (Rabi)
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 5, 1))*(InputCf.index.date<=datetime.date(yr, 5, 31)),InputCategoriesUnique['LandUse'] == 'PaddyR'] = 0.48 # Maintenance + Nursery (Rabi)
     
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 6, 1))*(InputCf.index.date<=datetime.date(yr, 6, 30)),InputCategoriesUnique['LandUse'] == 'PaddyK'] = 0.51 # Nursery (Kharif)     
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 7, 1))*(InputCf.index.date<=datetime.date(yr, 9, 30)),InputCategoriesUnique['LandUse'] == 'PaddyK'] = 0.51 # Irrigation (Kharif) 
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 10, 1))*(InputCf.index.date<=datetime.date(yr, 12, 31)),InputCategoriesUnique['LandUse'] == 'PaddyK'] = 1 # Harvesting (Kharif) + No pumping
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 1, 1))*(InputCf.index.date<=datetime.date(yr, 5, 31)),InputCategoriesUnique['LandUse'] == 'PaddyK'] = 1 # No pumping

    InputCf.at[(InputCf.index.date>=datetime.date(yr, 6, 1))*(InputCf.index.date<=datetime.date(yr, 10, 31)),InputCategoriesUnique['LandUse'] == 'Vegetable'] = 0.26 # Kharif
     