# -*- coding: utf-8 -*-
"""
Author: Madeleine NICOLAS       
        madeleine.nicolas@univ-rennes1.fr

Purpose: Physically based hydraulic model from Dewandel et al. 2008 applied to a fine grid accounting for Land Use and
         Soil class, parameters are obtaines from De Condappa (2005)

Output: Yearly recharge at each grid cell

Requirements: HydModel_inputs is run to obtain input data

Remarks: 

Date:   April 2018
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from scipy.stats import hmean, gmean
import shapefile   
import time
import winsound 

t = time.time()

os.chdir(r"C:\Users\Madeleine\Desktop\MAHESHWARAM\01 - SOIL MOISTURE MODEL\03 - PYTHON CODES") # Sets working directory
InputDir = r"C:\Users\Madeleine\Desktop\MAHESHWARAM\01 - SOIL MOISTURE MODEL\02 - INPUT DATA"
OutputDir = r"C:\Users\Madeleine\Desktop\MAHESHWARAM\01 - SOIL MOISTURE MODEL\01 - MODELS\02 - HydMod\01 - Outputs"

from HydModel__def import HydMod_fun
import HydModel_FineGrid_inputs as HydModel_inputs

ScenName = '2Lay_4pt5cm'

#==============================================================================
# INPUT PARAM TABLES 
#==============================================================================

SoilTypesTable = ['Alfisols 1','Alfisols 2','Inceptisols','Entisols','Tank']
ThreshTab = [0.12, 0.045] # Height difference between top of soil and edge bordering field (m)

############################## LAYER 1 ########################################
###################### Eluviation horizon hE ##################################

## GENERAL
SoilThickTab1 = [0.20,  0.20, 0.10, 0.20, 0.50] # Soil thickness (hE) (m)
#SoilThickTab1 = np.multiply(SoilThickTab1,2)

## DEFAULT
KsTab1 = np.multiply([1.87E-5, 2.48E-5, 1.24E-5, 4.23E-6, 2.5E-7],86400) # Saturated Hydraulic conductivity (m/day)
ThetaSTab1 = [0.32, 0.36, 0.35, 0.36, 0.60] # Saturated soil moisture content (m3/m3)
LambdaTab1 = [0.175, 0.193, 0.129, 0.117, 0.033] # Pore-size-distribution parameter
hbcTab1 = [-0.089, -0.111, -0.060, -0.037, -0.312] # Air-entry pressure head (m)
EtaTab1 = [14.8, 13.3, 19.3, 20.2, 63.7] # Texture-dependent conductivity shape parameter

## PADDY
KsPaddy = 2.50E-7 * 86400
ThetaSPaddy = 0.36
LambdaPaddy = 0.099
hbcPaddy = -0.05
EtaPaddy = 23.2

## FOREST
KsTabForest1 = np.multiply([2.62E-5, 3.15E-5, 2.09E-5, 1.35E-5, 2.5E-7],86400)
ThetaSTabForest1 = [0.25, 0.32, 0.27, 0.36, 0.60]
LambdaTabForest1 = [0.177, 0.193, 0.129, 0.117, 0.033]
hbcTabForest1 = [-0.117, -0.132, -0.098, -0.066, -0.312]
EtaTabForest1 = [14.8, 13.3, 19.3, 20.2, 63.7]

## SCRUB
KsTabScrub1 = np.multiply([3.27E-5, 3.71E-5, 2.83E-5, 2.195E-5, 2.5E-7],86400)
ThetaSTabScrub1 = [0.25, 0.28, 0.27, 0.28, 0.60]
LambdaTabScrub1 = [0.177, 0.193, 0.129, 0.117, 0.033]
hbcTabScrub1 = [-0.136, -0.142, -0.125, -0.102, -0.312]
EtaTabScrub1 = [14.8, 13.3, 19.3, 20.2, 63.7]

############################## LAYER 2 ########################################
###################### Illuviation horizon hB #################################

## GENERAL
SoilThickTab2 = [1.20, 1.40, 1.13, 2.20, 2.20] # Soil thickness (hB) (m)
#SoilThickTab2 = np.multiply(SoilThickTab2,2)

## DEFAULT
KsTab2 = np.multiply([5.19E-6, 6.06E-6, 1.08E-5, 3.10E-6, 7.50E-7],86400) # Saturated Hydraulic conductivity (m/s)
ThetaSTab2 = [0.37, 0.37, 0.35, 0.35, 0.55] # Saturated soil moisture content (m3/m3)
LambdaTab2 = [0.077, 0.100, 0.096, 0.080, 0.052] # Pore-size-distribution parameter
hbcTab2 = [-0.041, -0.053, -0.123, -0.036, -0.055] # Air-entry pressure head (m)
EtaTab2 = [29.5, 23.0, 24.1, 28.2, 41.7] # Texture-dependent conductivity shape parameter

## FOREST
KsTabForest2 = np.multiply([1.07E-05, 1.14E-05, 1.62E-05, 6.67E-06, 7.50E-07],86400)
ThetaSTabForest2 = [0.26, 0.37, 0.27, 0.35, 0.55]
LambdaTabForest2 = [0.076, 0.100, 0.096, 0.080, 0.052]
hbcTabForest2 = [-0.122, -0.132, -0.201, -0.089, -0.055]
EtaTabForest2 = [29.7, 23.0, 24.1, 28.2, 41.7]

## SCRUB
KsTabScrub2 = np.multiply([1.61E-05, 1.68E-05, 2.16E-05, 1.22E-05, 7.50E-7],86400)
ThetaSTabScrub2 = [0.26, 0.26, 0.27, 0.35, 0.55]
LambdaTabScrub2 = [0.076, 0.100, 0.096, 0.080, 0.052]
hbcTabScrub2 = [-0.201, -0.210, -0.279, -0.144, -0.055]
EtaTabScrub2 = [29.7, 23.0, 24.1, 28.2, 41.7]

#==============================================================================
# INPUT DATA 
#==============================================================================

# Coordinates
InputCoord = pd.read_csv(InputDir + "\INPUT_coord_finegrid.txt",sep='\t',header=0,index_col=0)
InputCoordRef = pd.read_csv(InputDir + "\INPUT_coord.txt",sep='\t',header=0,index_col=2)

# Inputs (PG + Rainfall)
InputClimate = HydModel_inputs.InputClimate
InputPG = HydModel_inputs.InputPG 
InputPG_Yrly = InputPG.resample("A").sum()
InputPG_Yrly.index = InputPG_Yrly.index.year

# Evapotranspiration
InputRET = HydModel_inputs.InputRET # meters
#InputRET = InputRET * Kp

# Partition coef between IRF and natural Recharge
InputCf = HydModel_inputs.InputCf 

# Input categories for discretizing basin
InputCategories = HydModel_inputs.InputCategories
InputCategoriesUnique = HydModel_inputs.InputCategoriesUnique

#==============================================================================
# FUNCTION RUN 
#==============================================================================

## Initialize dataframes
OutputHydMod_Cat_Rech1 = pd.DataFrame(columns=InputCategoriesUnique.index,index = range(len(InputClimate))) # Layer 1
OutputHydMod_Cat_Deficit1 = pd.DataFrame(columns=InputCategoriesUnique.index,index = range(len(InputClimate)))
OutputHydMod_Cat_AET1 = pd.DataFrame(columns=InputCategoriesUnique.index,index = range(len(InputClimate)))
OutputHydMod_Cat_Rnff1 = pd.DataFrame(columns=InputCategoriesUnique.index,index = range(len(InputClimate)))
OutputHydMod_Cat_Rech2 = pd.DataFrame(columns=InputCategoriesUnique.index,index = range(len(InputClimate))) # Layer 2
OutputHydMod_Cat_Deficit2 = pd.DataFrame(columns=InputCategoriesUnique.index,index = range(len(InputClimate)))
OutputHydMod_Cat_AET2 = pd.DataFrame(columns=InputCategoriesUnique.index,index = range(len(InputClimate)))

OutputVar_Cat_Kh1 = pd.DataFrame(columns=InputCategoriesUnique.index,index = range(len(InputClimate))) # Layer 1
OutputVar_Cat_hUnsat1 = pd.DataFrame(columns=InputCategoriesUnique.index,index = range(len(InputClimate)))
OutputVar_Cat_Theta1 = pd.DataFrame(columns=InputCategoriesUnique.index,index = range(len(InputClimate)))
OutputVar_Cat_Kh2 = pd.DataFrame(columns=InputCategoriesUnique.index,index = range(len(InputClimate))) # Layer 2
OutputVar_Cat_hUnsat2 = pd.DataFrame(columns=InputCategoriesUnique.index,index = range(len(InputClimate)))
OutputVar_Cat_Theta2 = pd.DataFrame(columns=InputCategoriesUnique.index,index = range(len(InputClimate)))

for cat in InputCategoriesUnique.index[28]:
    
    ## Select parameters for specific category
    
    SoilIndex = SoilTypesTable.index(InputCategoriesUnique['SoilClass'][cat])
    
    SoilThick1 = SoilThickTab1[SoilIndex]
    SoilThick2 = SoilThickTab2[SoilIndex]  
    
    Thresh2 = 100 # No runoff from Lay2
    
    if (InputCategoriesUnique['LandUse'][cat]=='PaddyK') | (InputCategoriesUnique['LandUse'][cat]=='PaddyR'):
        
        ##### LAYER 1 #####    
        Thresh1 = ThreshTab[0]
        Ks1 = KsPaddy
        ThetaS1 = ThetaSPaddy
        Lambda1 = LambdaPaddy
        hbc1 = hbcPaddy
        Eta1 = EtaPaddy
        
        ##### LAYER 2 #####   
        Ks2 = KsTab2[SoilIndex]
        ThetaS2 = ThetaSTab2[SoilIndex]
        Lambda2 = LambdaTab2[SoilIndex]
        hbc2 = hbcTab2[SoilIndex]
        Eta2 = EtaTab2[SoilIndex]
        
    else:
        
        Thresh1 = ThreshTab[1]
        
        if  (InputCategoriesUnique['LandUse'][cat] == 'Forest'):   
            
            ##### LAYER 1 #####   
            Ks1 = KsTabForest1[SoilIndex]
            ThetaS1 = ThetaSTabForest1[SoilIndex]
            Lambda1 = LambdaTabForest1[SoilIndex]
            hbc1 = hbcTabForest1[SoilIndex]
            Eta1 = EtaTabForest1[SoilIndex]
            
            ##### LAYER 2 #####   
            Ks2 = KsTabForest2[SoilIndex]
            ThetaS2 = ThetaSTabForest2[SoilIndex]
            Lambda2 = LambdaTabForest2[SoilIndex]
            hbc2 = hbcTabForest2[SoilIndex]
            Eta2 = EtaTabForest2[SoilIndex]
            
        elif (InputCategoriesUnique['LandUse'][cat] == 'Scrub'):
            
            ##### LAYER 1 #####  
            Ks1 = KsTabScrub1[SoilIndex]
            ThetaS1 = ThetaSTabScrub1[SoilIndex]
            Lambda1 = LambdaTabScrub1[SoilIndex]
            hbc1 = hbcTabScrub1[SoilIndex]
            Eta1 = EtaTabScrub1[SoilIndex]
            
            ##### LAYER 2 #####  
            Ks2 = KsTabScrub2[SoilIndex]
            ThetaS2 = ThetaSTabScrub2[SoilIndex]
            Lambda2 = LambdaTabScrub2[SoilIndex]
            hbc2 = hbcTabScrub2[SoilIndex]
            Eta2 = EtaTabScrub2[SoilIndex]
            
        else:

            ##### LAYER 1 #####  
            Ks1 = KsTab1[SoilIndex]
            ThetaS1 = ThetaSTab1[SoilIndex]
            Lambda1 = LambdaTab1[SoilIndex]
            hbc1 = hbcTab1[SoilIndex]
            Eta1 = EtaTab1[SoilIndex]
    
            ##### LAYER 2 #####  
            Ks2 = KsTab2[SoilIndex]
            ThetaS2 = ThetaSTab2[SoilIndex]
            Lambda2 = LambdaTab2[SoilIndex]
            hbc2 = hbcTab2[SoilIndex]
            Eta2 = EtaTab2[SoilIndex]
    
    ############################# LAYER 1 #####################################
    
    ## Get input data for specific category
    
    InputData1 = pd.DataFrame(columns=['R','PG','RET'])
    InputData1['R'] = InputClimate['Rainfall_mm']/1000
    InputData1['PG'] = InputPG[cat]
    InputData1['RET'] = InputRET[cat]
    InputData1.index = pd.to_datetime(InputData1.index) 
    
    ## Run model
    
    HydModInter1 = HydMod_fun(Ks1,SoilThick1,ThetaS1,Thresh1,Lambda1,hbc1,Eta1,InputData1) # Hydraulic model run for each category
    
    ## Store data    
    
    # Fluxes
    OutputHydMod_Cat_Rech1[cat] = HydModInter1['qout']
    OutputHydMod_Cat_Rnff1[cat] = HydModInter1['Runoff']
    OutputHydMod_Cat_Deficit1[cat] = HydModInter1['Deficit']
    OutputHydMod_Cat_AET1[cat] = HydModInter1['AET']
    
    # Variables
    OutputVar_Cat_Kh1[cat] = HydModInter1['Kh']
    OutputVar_Cat_hUnsat1[cat] = HydModInter1['hUnsat']
    OutputVar_Cat_Theta1[cat] = HydModInter1['Theta1']
    
    ############################# LAYER 2 #####################################
    
    InputData2 = pd.DataFrame(columns=['R','RET'])
    InputData2['R'] = HydModInter1['qout']
    InputData2['RET'] = HydModInter1['Deficit']
    InputData2.index = pd.to_datetime(InputData1.index) 
    
    ## Run model
    HydModInter2 = HydMod_fun(Ks2,SoilThick2,ThetaS2,Thresh2,Lambda2,hbc2,Eta2,InputData2) # Hydraulic model run for each category
    
    ##Store data

    # Fluxes
    OutputHydMod_Cat_Rech2[cat] = HydModInter2['qout']
    OutputHydMod_Cat_Deficit2[cat] = HydModInter2['Deficit']
    OutputHydMod_Cat_AET2[cat] = HydModInter2['AET']
    
    # Variables
    OutputVar_Cat_Kh2[cat] = HydModInter2['Kh']
    OutputVar_Cat_hUnsat2[cat] = HydModInter2['hUnsat']
    OutputVar_Cat_Theta2[cat] = HydModInter2['Theta1']
    
#==============================================================================
# CLEAN UP DATA
#==============================================================================

## Recharge 

OutputHydMod_Cat_Rech = OutputHydMod_Cat_Rech2 * 1000.0 # Recharge is converted to mm
OutputHydMod_Cat_Rech[OutputHydMod_Cat_Rech<1] = 0 # Remove negligible values
OutputHydMod_Cat_Rech.index = pd.to_datetime(InputClimate.index) # Indexes are set

## Runoff

OutputHydMod_Cat_Rnff = OutputHydMod_Cat_Rnff1 * 1000.0 # Runoff is converted to mm
OutputHydMod_Cat_Rnff.index = pd.to_datetime(InputClimate.index) # Indexes are set
OutputHydMod_Cat_Rnff_Yrly = OutputHydMod_Cat_Rnff.resample("A").sum() # Resample at yearly interval
OutputHydMod_Cat_Rnff_Yrly.index = OutputHydMod_Cat_Rnff_Yrly.index.year # Set indexes as years

## Deficit

OutputHydMod_Cat_Deficit = OutputHydMod_Cat_Deficit2 * 1000.0 # Deficit is converted to mm
OutputHydMod_Cat_Deficit.index = pd.to_datetime(InputClimate.index) # Indexes are set

## AET

OutputHydMod_Cat_AET = OutputHydMod_Cat_AET1 + OutputHydMod_Cat_AET2 
OutputHydMod_Cat_AET.index = pd.to_datetime(InputClimate.index) # Indexes are set
OutputHydMod_Cat_AET_Yrly = OutputHydMod_Cat_AET.resample("A").sum() # Resample at yearly interval
OutputHydMod_Cat_AET_Yrly.index = OutputHydMod_Cat_AET_Yrly.index.year # Set indexes as years

## Separate natural recharge and return flow

Cf = np.divide(OutputHydMod_Cat_Rech,InputPG[InputPG>0]*1000) # Calculate Return Flow Coefficient Cf
CfCorr = Cf[Cf>InputCf]-InputCf # Obtain coefficients to multiply pumping and obtain natural recharge
OutputHydMod_Cat_RechCorr = CfCorr * InputPG * 1000
OutputHydMod_Cat_RechCorr[OutputHydMod_Cat_RechCorr.isnull()] = 0
a = OutputHydMod_Cat_RechCorr[[i for i, x in enumerate((InputPG != 0).any(axis=0)) if x]] # Intermediate step to obtain corrected values for pumped areas
b = OutputHydMod_Cat_Rech[[i for i, x in enumerate(~(InputPG != 0).any(axis=0)) if x]]# Intermediate step to obtain raw values for unpumped areas
OutputHydMod_Cat_RechCorr = pd.DataFrame(a.join(b)[OutputHydMod_Cat_Rech.columns]).astype(float)

OutputHydMod_Cat_RechCorr_Yrly = OutputHydMod_Cat_RechCorr.resample("A").sum() # Resample at yearly interval
OutputHydMod_Cat_RechCorr_Yrly.index = OutputHydMod_Cat_RechCorr_Yrly.index.year # Set indexes as years

## Variables

# LAYER 1
#Hydraulic conductivity
OutputVar_Cat_Kh_mask1 = ma.masked_array(OutputVar_Cat_Kh1, mask=OutputVar_Cat_Kh1<1e-5) # Remove negligible values
OutputVar_Cat_Kh_mean1 = pd.DataFrame(hmean(OutputVar_Cat_Kh_mask1,axis=0)).transpose() # Harmonic mean because variables exist across several orders of magnitude
OutputVar_Cat_Kh_mean1.columns = InputCategoriesUnique.index
#Hydraulic pressure
OutputVar_Cat_hUnsat_mask1 = ma.masked_array(OutputVar_Cat_hUnsat1, mask=OutputVar_Cat_hUnsat1<-1e+3) # Remove negligible values
OutputVar_Cat_hUnsat_mean1 = pd.DataFrame(OutputVar_Cat_hUnsat_mask1.mean(axis=0)).transpose() # Arithmetic mean 
OutputVar_Cat_hUnsat_mean1.columns = InputCategoriesUnique.index
#Soil moisture
OutputVar_Cat_Theta_mean1 = pd.DataFrame(OutputVar_Cat_Theta1.mean(axis=0)).transpose()

#LAYER 2
#Hydraulic conductivity
OutputVar_Cat_Kh_mask2 = ma.masked_array(OutputVar_Cat_Kh2, mask=OutputVar_Cat_Kh2<1e-5) # Remove negligible values
OutputVar_Cat_Kh_mean2 = pd.DataFrame(hmean(OutputVar_Cat_Kh_mask2,axis=0)).transpose() # Harmonic mean because variables exist across several orders of magnitude
OutputVar_Cat_Kh_mean2.columns = InputCategoriesUnique.index
#Hydraulic pressure
OutputVar_Cat_hUnsat_mask2 = ma.masked_array(OutputVar_Cat_hUnsat2, mask=OutputVar_Cat_hUnsat2<-1e+3) # Remove negligible values
OutputVar_Cat_hUnsat_mean2 = pd.DataFrame(OutputVar_Cat_hUnsat_mask2.mean(axis=0)).transpose() # Arithmetic mean 
OutputVar_Cat_hUnsat_mean2.columns = InputCategoriesUnique.index
#Soil moisture
OutputVar_Cat_Theta_mean2 = pd.DataFrame(OutputVar_Cat_Theta2.mean(axis=0)).transpose()

#==============================================================================
# SPATIALIZE DATA 
#==============================================================================

# Initialize dataframes

OutputHydMod_GridRaw_Rech = pd.DataFrame(index=InputCategories.index, columns=OutputHydMod_Cat_RechCorr_Yrly.index)
OutputHydMod_GridRaw_Rnff = pd.DataFrame(index=InputCategories.index, columns=OutputHydMod_Cat_Rnff_Yrly.index)
OutputHydMod_GridRaw_AET = pd.DataFrame(index=InputCategories.index, columns=OutputHydMod_Cat_AET_Yrly.index)
OutputHydMod_GridRaw_PG = pd.DataFrame(index=InputCategories.index, columns=OutputHydMod_Cat_Rnff_Yrly.index)

OutputVar_GridRaw_Kh1 = pd.DataFrame(index=InputCategories.index, columns=["Kh"]) # Layer 1
OutputVar_GridRaw_hUnsat1 = pd.DataFrame(index=InputCategories.index, columns=["hUnsat"])
OutputVar_GridRaw_Theta1 = pd.DataFrame(index=InputCategories.index, columns=["Theta"])
OutputVar_GridRaw_Kh2 = pd.DataFrame(index=InputCategories.index, columns=["Kh"]) # Layer 2
OutputVar_GridRaw_hUnsat2 = pd.DataFrame(index=InputCategories.index, columns=["hUnsat"])
OutputVar_GridRaw_Theta2 = pd.DataFrame(index=InputCategories.index, columns=["Theta"])

# Obtain appropriate recharge for soil and land use combination at each grid cell

for j in range(len(InputCategories.index)): 

    MatchingIndex = InputCategoriesUnique.index[[i for i, x in enumerate((InputCategoriesUnique['SoilClass']==InputCategories['SoilClass'][j]) & (InputCategoriesUnique['LandUse']==InputCategories['LandUse'][j])) if x]]
    
    #Fluxes
    MatchingData = OutputHydMod_Cat_RechCorr_Yrly[MatchingIndex].transpose()
    MatchingDataRnff = OutputHydMod_Cat_Rnff_Yrly[MatchingIndex].transpose()  
    MatchingDataAET = OutputHydMod_Cat_AET_Yrly[MatchingIndex].transpose() 
    MatchingDataPG = InputPG_Yrly[MatchingIndex].transpose() 

    OutputHydMod_GridRaw_Rech.loc[j] = np.asarray(MatchingData)
    OutputHydMod_GridRaw_Rnff.loc[j] = np.asarray(MatchingDataRnff)
    OutputHydMod_GridRaw_AET.loc[j] = np.asarray(MatchingDataAET)
    OutputHydMod_GridRaw_PG.loc[j] = np.asarray(MatchingDataPG)
    
    #Variables
    MatchingDataVarKh1 = OutputVar_Cat_Kh_mean1[MatchingIndex] # Layer 1
    MatchingDataVarhUnsat1 = OutputVar_Cat_hUnsat_mean1[MatchingIndex]
    MatchingDataVarTheta1 = OutputVar_Cat_Theta_mean1[MatchingIndex]
    MatchingDataVarKh2 = OutputVar_Cat_Kh_mean2[MatchingIndex] # Layer 2
    MatchingDataVarhUnsat2 = OutputVar_Cat_hUnsat_mean2[MatchingIndex]
    MatchingDataVarTheta2 = OutputVar_Cat_Theta_mean2[MatchingIndex]
    
    OutputVar_GridRaw_Kh1.loc[j] = np.asarray(MatchingDataVarKh1) # Layer 1
    OutputVar_GridRaw_hUnsat1.loc[j] = np.asarray(MatchingDataVarhUnsat1)
    OutputVar_GridRaw_Theta1.loc[j] = np.asarray(MatchingDataVarTheta1)
    OutputVar_GridRaw_Kh2.loc[j] = np.asarray(MatchingDataVarKh2) # Layer 2
    OutputVar_GridRaw_hUnsat2.loc[j] = np.asarray(MatchingDataVarhUnsat2)
    OutputVar_GridRaw_Theta2.loc[j] = np.asarray(MatchingDataVarTheta2)    
    
# Fluxes    
OutputHydMod_Grid_Rech = InputCoord.join(OutputHydMod_GridRaw_Rech)
OutputHydMod_Grid_Rnff = InputCoord.join(OutputHydMod_GridRaw_Rnff)
OutputHydMod_Grid_AET = InputCoord.join(OutputHydMod_GridRaw_AET)
OutputHydMod_Grid_PG = InputCoord.join(OutputHydMod_GridRaw_PG)

#Variables
OutputVar_Grid_Kh1 = InputCoord.join(OutputVar_GridRaw_Kh1) # Layer 1
OutputVar_Grid_hUnsat1 = InputCoord.join(OutputVar_GridRaw_hUnsat1)
OutputVar_Grid_Theta1 = InputCoord.join(OutputVar_GridRaw_Theta1)
OutputVar_Grid_Kh2 = InputCoord.join(OutputVar_GridRaw_Kh2) # Layer 2
OutputVar_Grid_hUnsat2 = InputCoord.join(OutputVar_GridRaw_hUnsat2)
OutputVar_Grid_Theta2 = InputCoord.join(OutputVar_GridRaw_Theta2)

##==============================================================================
## INPUT OUTPUT BALANCE
##==============================================================================

# Input
BalancePG_Yrly = OutputHydMod_GridRaw_PG.mean()*1000.0
BalanceR_Yrly = InputData1['R'].resample("A").sum()*1000
BalanceR_Yrly.index = BalancePG_Yrly.index

# Output
BalanceRnff_Yrly = OutputHydMod_GridRaw_Rnff.mean()
BalanceRech_Yrly = OutputHydMod_GridRaw_Rech.mean()
BalanceAET_Yrly = OutputHydMod_GridRaw_AET.mean()*1000.0

#Balance
BalanceTab = pd.DataFrame(columns = ['Runoff','TotRech','AET','Rainfall','PG'])
BalanceTab['Runoff']  = BalanceRnff_Yrly
BalanceTab['TotRech'] = BalanceRech_Yrly
BalanceTab['AET'] = BalanceAET_Yrly
BalanceTab['Rainfall'] = BalanceR_Yrly
BalanceTab['PG'] = BalancePG_Yrly

BalanceTotalDiff_Yrly = BalancePG_Yrly + BalanceR_Yrly - (BalanceAET_Yrly + BalanceRnff_Yrly + BalanceRech_Yrly)
BalanceTotalDiff = BalanceTotalDiff_Yrly[2:].sum()
#plt.plot(BalanceTotalDiff_Yrly)

##==============================================================================
## SAVE OUTPUTS
##==============================================================================

OutputHydMod_Mean_Rnff = pd.DataFrame(OutputHydMod_Cat_Rnff.mean(axis=1))

os.chdir(OutputDir) # Sets working directory

writer = pd.ExcelWriter(ScenName + "_OUTPUT_HydModel.xlsx") # Excel file will all data
OutputHydMod_Grid_Rech.to_excel(writer,'HydModYearly', index=True)
OutputHydMod_Cat_Rech.to_excel(writer,'HydModPerCat', index=True)
OutputHydMod_Mean_Rnff.to_excel(writer,'MeanRnffTimeSeries', index=True)
BalanceTab.to_excel(writer,'BalanceRaw', index=True)
writer.save()

#Variables (yearly)
writer = pd.ExcelWriter(ScenName + "_OUTPUT_VariablesMean.xlsx") # Excel file will all data
OutputVar_Grid_Kh1.to_excel(writer,'Kh_L1', index=True)
OutputVar_Grid_hUnsat1.to_excel(writer,'hUnsat_L1', index=True)
OutputVar_Grid_Theta1.to_excel(writer, 'Theta_L1', index=True,header='False')
OutputVar_Grid_Kh2.to_excel(writer,'Kh_L2', index=True)
OutputVar_Grid_hUnsat2.to_excel(writer,'hUnsat_L2', index=True)
OutputVar_Grid_Theta2.to_excel(writer, 'Theta_L2', index=True,header='False')
writer.save()

#Variables (daily)
writer = pd.ExcelWriter(ScenName + "_OUTPUT_VariablesDaily.xlsx") # Excel file will all data
OutputVar_Cat_Kh1.to_excel(writer,'Kh_L1', index=True)
OutputVar_Cat_hUnsat1.to_excel(writer,'hUnsat_L1', index=True)
OutputVar_Cat_Theta1.to_excel(writer, 'Theta_L1', index=True,header='False')
OutputVar_Cat_Kh2.to_excel(writer,'Kh_L2', index=True)
OutputVar_Cat_hUnsat2.to_excel(writer,'hUnsat_L2', index=True)
OutputVar_Cat_Theta2.to_excel(writer, 'Theta_L2', index=True,header='False')
writer.save()

## Csv file for upscaling
OutputHydMod_Grid_export = OutputHydMod_Grid_Rech
OutputHydMod_Grid_export.columns.values[2::] = ['Y' + str(col) for col in OutputHydMod_Grid_export.columns[2::]] # Adding letters to column names for posterior GIS operation
OutputHydMod_Grid_export.to_csv(ScenName + "_HydModGrid.csv",index=False) # Csv file with gridded values  

#==============================================================================
# Print total elapsed time
#==============================================================================

elapsed = time.time() - t 
print ('Time elapsed: ' + '%.0f' % (elapsed/60) + ' min ' + '%.0f' % ((elapsed/60-int(elapsed/60))*60) + ' sec')

winsound.Beep(1000, 500) # Makes beep

#==============================================================================
# PLOT
#==============================================================================

while(True):
    Cond = raw_input("\nDo you want to plot figures?\n\n(Type y or n)\n\n")
    if (Cond == 'y') | (Cond == 'n'):
        # Break out of loop and continue code 
        break
    else: # Error message for invalid input
        print "\nERROR: Invalid Input"
        continue # Or "continue" to rerun loop, and ask user to re-enter data

if Cond == 'y':
    
    shpFilePath = InputDir + "\GIS\watershed.shp"
    os.chdir(r"C:\Users\Madeleine\Desktop\MAHESHWARAM\01 - SOIL MOISTURE MODEL\01 - MODELS\02 - HydMod\02 - Figures" )
    
    ## RECHARGE AND RUNOFF
    
    for yr in OutputHydMod_GridRaw_Rech.columns[5:]:
            
        ## RECHARGE
        plt.figure('HydMod' + str(yr))
        plt.suptitle(yr, fontsize=18)
        Val = OutputHydMod_Grid_Rech["Y" + str(yr)]
        
        listx=[]
        listy=[]
        test = shapefile.Reader(shpFilePath)
        for sr in test.shapeRecords():
            for xNew,yNew in sr.shape.points:
                listx.append(xNew)
                listy.append(yNew)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.plot(listx,listy,color='k')
        
        plt.title(ScenName + "_Rech")
        plt.xlim((min(listx)-400,max(listx)+400))
        plt.ylim((min(listy)-400,max(listy)+400))
        plt.scatter(OutputHydMod_Grid_Rech['X'],OutputHydMod_Grid_Rech['Y'],c=Val, cmap='Blues',s=100,vmin=-0,vmax=+1500,marker='s',edgecolor='face') # Default is jet
        plt.colorbar()
        
        plt.savefig(str(yr) + '_' + ScenName + '_HydMod.png', dpi = 900)
        plt.show()
        
        ## RUNOFF
        plt.figure('HydModRunoff' + str(yr))
        plt.suptitle(yr, fontsize=18)
        Val = OutputHydMod_Grid_Rnff[yr]
        
        listx=[]
        listy=[]
        test = shapefile.Reader(shpFilePath)
        for sr in test.shapeRecords():
            for xNew,yNew in sr.shape.points:
                listx.append(xNew)
                listy.append(yNew)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.plot(listx,listy,color='k')
        
        plt.title("Runoff")
        plt.xlim((min(listx)-400,max(listx)+400))
        plt.ylim((min(listy)-400,max(listy)+400))
        plt.scatter(OutputHydMod_Grid_Rnff['X'],OutputHydMod_Grid_Rnff['Y'],c=Val, cmap='Blues',s=100,vmin=-0,vmax=+1500,marker='s',edgecolor='face') # Default is jet
        plt.colorbar()
        
        plt.savefig(str(yr) + '_' + ScenName + '_HydMod_Runoff.png', dpi = 900)
        plt.show()

