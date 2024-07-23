# -*- coding: utf-8 -*-
"""
Author: Madeleine NICOLAS       
        madeleine.nicolas@univ-rennes1.fr

Purpose: Plot inputs at a fine grid scale used for HydMod (PG, RET, Cf)

Output: Figures

Requirements: 

Remarks: 

Date:   August 2018
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import shapefile

import HydModel_inputs

InputDir = r"C:\Users\Madeleine\Desktop\MAHESHWARAM\01 - SOIL MOISTURE MODEL\02 - INPUT DATA"

#==============================================================================
# INPUT PARAM TABLES 
#==============================================================================

SoilTypesTable = ['Alfisols 1','Alfisols 2','Inceptisols','Entisols','Tank']

## DEFAULT
KsTab = np.multiply([1.87E-5, 2.48E-5, 1.24E-5, 4.23E-6, 2.5E-7],86400) # Saturated Hydraulic conductivity (m/s)
SoilThickTab = [0.20,  0.20, 0.10, 0.20, 0.50] # Soil thickness (m)
SigmaSTab = [0.32, 0.36, 0.35, 0.36, 0.60] # Saturated soil moisture content (m3/m3)
ThreshTab = [0.12, 0.02] # Height difference between top of soil and edge bordering field (m)
LambdaTab = [0.175, 0.193, 0.129, 0.117, 0.033] # Pore-size-distribution parameter
hbcTab = [-0.089, -0.111, -0.060, -0.037, -0.312] # Air-entry pressure head (m)
EtaTab = [14.8, 13.3, 19.3, 20.2, 63.7] # Texture-dependent conductivity shape parameter

## PADDY
KsPaddy = 2.50E-7 * 86400
SigmaSPaddy = 0.36
LambdaPaddy = 0.099
hbcPaddy = -0.05
EtaPaddy = 23.2

## FOREST
KsTabForest = np.multiply([2.62E-5, 3.15E-5, 2.09E-5, 1.35E-5, 2.5E-7],86400)
SigmaSTabForest = [0.25, 0.32, 0.27, 0.36, 0.60]
LambdaTabForest = [0.177, 0.193, 0.129, 0.117, 0.033]
hbcTabForest = [-0.117, -0.132, -0.098, -0.066, -0.312]
EtaTabForest = [14.8, 13.3, 19.3, 20.2, 63.7]

## SCRUB
KsTabScrub = np.multiply([3.27E-5, 3.71E-5, 2.83E-5, 2.195E-5, 2.5E-7],86400)
SigmaSTabScrub = [0.25, 0.28, 0.27, 0.28, 0.60]
LambdaTabScrub = [0.177, 0.193, 0.129, 0.117, 0.033]
hbcTabScrub = [-0.136, -0.142, -0.125, -0.102, -0.312]
EtaTabScrub = [14.8, 13.3, 19.3, 20.2, 63.7]

#==============================================================================
# INPUT DATA 
#==============================================================================

# Coordinates
InputCoord = pd.read_csv(InputDir + "\INPUT_coord_finegrid.txt",sep='\t',header=0,index_col=0)

# Inputs (PG + Rainfall)
InputClimate = HydModel_inputs.InputClimate
InputPG = HydModel_inputs.InputPG 

# Evapotranspiration
InputRET = HydModel_inputs.InputRET # meters

# Partition coef between IRF and natural Recharge
InputCf = HydModel_inputs.InputCf 

# Input categories for discretizing basin
InputCategories = HydModel_inputs.InputCategories
InputCategoriesUnique = HydModel_inputs.InputCategoriesUnique

#==============================================================================
# DATA PER CATEGORY
#==============================================================================

#########################
## INTER ANNUAL VALUES ##
#########################

# Initialize dataframes
InputCfInterYearCat = pd.DataFrame(columns=InputCategoriesUnique.index)
InputPGInterYearCat = pd.DataFrame(columns=InputCategoriesUnique.index)
InputRETInterYearCat = pd.DataFrame(columns=InputCategoriesUnique.index)

for cat in InputCategoriesUnique.index : 
    
    # Input Cf
    InputCfInterYear_Inter = pd.Series(InputCf[cat],InputCf.index,dtype=float)
    InputCfInterYearCat[cat] = InputCfInterYear_Inter.groupby([InputCfInterYear_Inter.index.month,InputCfInterYear_Inter.index.day]).mean() 
    
    # Input PG   
    InputPGInterYear_Inter = pd.Series(InputPG[cat],InputPG.index,dtype=float)
    InputPGInterYearCat[cat] = InputPGInterYear_Inter.groupby([InputPGInterYear_Inter.index.month,InputPGInterYear_Inter.index.day]).mean() 
  
    #InputRET
    InputRETInterYear_Inter = pd.Series(InputRET[cat],InputRET.index,dtype=float)
    InputRETInterYearCat[cat] = InputRETInterYear_Inter.groupby([InputRETInterYear_Inter.index.month,InputRETInterYear_Inter.index.day]).mean() 

InputCfInterYearCat = InputCfInterYearCat.drop((2L, 29L)) # Remove Feb 29
InputPGInterYearCat = InputPGInterYearCat.drop((2L, 29L))
InputRETInterYearCat = InputRETInterYearCat.drop((2L, 29L))

###################
## YEARLY VALUES ##
###################

InputCfYearCat = pd.DataFrame(InputCf,dtype=float).resample("A").mean()
InputPGYearCat = InputPG.resample("A").sum()*1000
InputRETYearCat = InputRET.resample("A").sum()*1000

InputCfYearCat.index = InputCfYearCat.index.year
InputPGYearCat.index = InputPGYearCat.index.year
InputRETYearCat.index = InputRETYearCat.index.year

#==============================================================================
# SPATIALIZE DATA 
# Obtain values at each grid cell
#==============================================================================

#Initialize dataframes
InputCfYearGrid = pd.DataFrame(index=InputCategories.index, columns=InputCfYearCat.index)
InputPGYearGrid = pd.DataFrame(index=InputCategories.index, columns=InputPGYearCat.index)
InputRETYearGrid = pd.DataFrame(index=InputCategories.index, columns=InputRETYearCat.index)

# Obtain appropriate recharge for soil and land use combination at each grid cell

for j in range(len(InputCategories.index)): 

    MatchingIndex = InputCategoriesUnique.index[[i for i, x in enumerate((InputCategoriesUnique['SoilClass']==InputCategories['SoilClass'][j]) & (InputCategoriesUnique['LandUse']==InputCategories['LandUse'][j])) if x]]
    
    MatchingDataCf = InputCfYearCat[MatchingIndex].transpose()
    MatchingDataPG = InputPGYearCat[MatchingIndex].transpose()
    MatchingDataRET = InputRETYearCat[MatchingIndex].transpose()
    
    InputCfYearGrid.loc[j] = np.asarray(MatchingDataCf)
    InputPGYearGrid.loc[j] = np.asarray(MatchingDataPG)
    InputRETYearGrid.loc[j] = np.asarray(MatchingDataRET)
    
InputCfYearGrid = InputCoord.join(pd.DataFrame(InputCfYearGrid.mean(axis=1)))
InputCfYearGrid.columns = ['X','Y','Cf']

InputPGYearGrid  = InputCoord.join(pd.DataFrame(InputPGYearGrid.mean(axis=1)))
InputPGYearGrid.columns = ['X','Y','PG']

InputRETYearGrid  = InputCoord.join(InputRETYearGrid)

##==============================================================================
## PLOT
##==============================================================================

shpFilePath = InputDir + "\GIS\watershed.shp"
os.chdir(r"C:\Users\Madeleine\Desktop\MAHESHWARAM\01 - SOIL MOISTURE MODEL\01 - MODELS\Figures")

## RETURN FLOW COEFFICIENT

plt.figure('Cf')
Val = InputCfYearGrid['Cf']

listx=[]
listy=[]
test = shapefile.Reader(shpFilePath)
for sr in test.shapeRecords():
    for xNew,yNew in sr.shape.points:
        listx.append(xNew)
        listy.append(yNew)
plt.gca().set_aspect('equal', adjustable='box')
plt.plot(listx,listy,color='k')

plt.title("Cf")
plt.xlim((min(listx)-400,max(listx)+400))
plt.ylim((min(listy)-400,max(listy)+400))
plt.scatter(InputCfYearGrid['X'],InputCfYearGrid['Y'],c=Val, cmap='YlGnBu',s=100,vmin=-0,vmax=+1,marker='s',edgecolor='face') # Default is jet
plt.colorbar()

plt.savefig('InputCf.png', dpi = 900)
plt.show()

## PUMPING

plt.figure('Pumping')
Val = InputPGYearGrid['PG']

listx=[]
listy=[]
test = shapefile.Reader(shpFilePath)
for sr in test.shapeRecords():
    for xNew,yNew in sr.shape.points:
        listx.append(xNew)
        listy.append(yNew)
plt.gca().set_aspect('equal', adjustable='box')
plt.plot(listx,listy,color='k')

plt.title("Pumping")
plt.xlim((min(listx)-400,max(listx)+400))
plt.ylim((min(listy)-400,max(listy)+400))
plt.scatter(InputPGYearGrid['X'],InputPGYearGrid['Y'],c=Val, cmap='Blues',s=100,vmin=-0,vmax=+3000,marker='s',edgecolor='face') # Default is jet
plt.colorbar()

plt.savefig('InputPG.png', dpi = 900)
plt.show()

## EVAPORATION

for yr in InputRETYearGrid.columns[2:]:
        
    ## RECHARGE
    plt.figure('RET' + str(yr))
    plt.suptitle(yr, fontsize=18)
    Val = InputRETYearGrid[yr]
    
    listx=[]
    listy=[]
    test = shapefile.Reader(shpFilePath)
    for sr in test.shapeRecords():
        for xNew,yNew in sr.shape.points:
            listx.append(xNew)
            listy.append(yNew)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.plot(listx,listy,color='k')
    
    plt.title("Evaporation")
    plt.xlim((min(listx)-400,max(listx)+400))
    plt.ylim((min(listy)-400,max(listy)+400))
    plt.scatter(InputRETYearGrid['X'],InputRETYearGrid['Y'],c=Val, cmap='Blues',s=100,vmin=-0,vmax=+3000,marker='s',edgecolor='face') # Default is jet
    plt.colorbar()
    
    plt.savefig('InputRET_' + str(yr) +'.png', dpi = 900)
    plt.show()