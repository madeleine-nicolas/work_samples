# -*- coding: utf-8 -*-
"""
Author: Madeleine NICOLAS       
        madeleine.nicolas@univ-rennes1.fr

Purpose: Physically based hydraulic model from Dewandel et al. 2008 applied to a fine grid accounting for Land Use and
         Soil class, parameters are obtaines from De Condappa (2005)

Output: _Yrly recharge at each grid cell

Requirements: HydModel_inputs is run to obtain input data

Remarks: 

Date:   April 2018
"""

import os
import pandas as pd
import numpy as np
import numpy.ma as ma
from scipy.stats import hmean, gmean
import matplotlib.pyplot as plt
import winsound
import shapefile   
import time 

t = time.time()

os.chdir(r"C:\Users\Madeleine\Desktop\MAHESHWARAM\01 - SOIL MOISTURE MODEL\03 - PYTHON CODES") # Sets working directory
InputDir = r"C:\Users\Madeleine\Desktop\MAHESHWARAM\01 - SOIL MOISTURE MODEL\02 - INPUT DATA"
OutputDir = r"C:\Users\Madeleine\Desktop\MAHESHWARAM\01 - SOIL MOISTURE MODEL\01 - MODELS\02 - HydMod\01 - Outputs"

from HydModel__def import HydMod_fun
import HydModel_FineGrid_inputs as HydModel_inputs

#==============================================================================
# INPUT PARAM TABLES 
#==============================================================================

SoilTypesTable = ['Alfisols 1','Alfisols 2','Inceptisols','Entisols','Tank']
#Kp = 0.9 # Regional evaporation coeficient

## GENERAL
ScenName = "1LayMean_2cm_long"
SoilThickTab = [1.4, 1.6, 1.2, 2.4, 2.7] # Soil thickness(m)
ThreshTab = [0.12, 0.02] # Height difference between top of soil and edge bordering field (m)

## DEFAULT
KsTab = np.multiply([5.79E-06, 6.69E-06, 1.09E-05, 3.17E-06, 5.47E-07],86400) # Saturated Hydraulic conductivity (m/s)
ThetaSTab = [0.36, 0.37, 0.35, 0.35, 0.56] # Saturated soil moisture content (m3/m3)
LambdaTab = [0.091, 0.112, 0.099, 0.083, 0.048] # Pore-size-distribution parameter
hbcTab = [-0.048, -0.060, -0.118, -0.036, -0.103] # Air-entry pressure head (m)
EtaTab = [27.4, 21.8, 23.7, 27.6, 45.8] # Texture-dependent conductivity shape parameter

## PADDY
KsPaddy = 2.50E-7 * 86400
ThetaSPaddy = 0.36
LambdaPaddy = 0.099
hbcPaddy = -0.05
EtaPaddy = 23.2

## FOREST
KsTabForest = np.multiply([1.17E-05, 1.24E-05, 1.65E-05, 6.97E-06, 5.47E-07],86400)
ThetaSTabForest = [0.26, 0.36, 0.27, 0.35, 0.56]
LambdaTabForest = [0.090, 0.112, 0.099, 0.083, 0.048]
hbcTabForest = [-0.122, -0.132, -0.193, -0.087, -0.103]
EtaTabForest = [27.5, 21.8, 23.7, 27.6, 45.8]

## SCRUB
KsTabScrub = np.multiply([1.74E-05, 1.80E-05, 2.20E-05, 1.27E-05, 5.47E-07],86400)
ThetaSTabScrub = [0.26, 0.26, 0.27, 0.34, 0.56]
LambdaTabScrub = [0.090, 0.112, 0.099, 0.083, 0.048]
hbcTabScrub = [-0.192, -0.202, -0.267, -0.140, -0.103]
EtaTabScrub = [27.5, 21.8, 23.7, 27.6, 45.8]

#==============================================================================
# INPUT DATA 
#==============================================================================

#Reference recharge from Adil
InputRech_Yrly = pd.read_csv(InputDir + "\INPUT_ref-rech_2002-2015.txt",sep='\t',header=0,index_col=0)

#Coordinates
InputCoord = pd.read_csv(InputDir + "\INPUT_coord_finegrid.txt",sep='\t',header=0,index_col=0)
InputCoordRef = pd.read_csv(InputDir + "\INPUT_coord.txt",sep='\t',header=0,index_col=2)

#Inputs (PG + Rainfall)
InputClimate = HydModel_inputs.InputClimate
InputPG = HydModel_inputs.InputPG
InputPG_Yrly = InputPG.resample("A").sum()
InputPG_Yrly.index = InputPG_Yrly.index.year

#Evapotranspiration
InputRET = HydModel_inputs.InputRET # meters
#InputRET = InputRET * Kp

#Partition coef between IRF and natural Recharge
InputCf = HydModel_inputs.InputCf 

#Input categories for discretizing basin
InputCategories = HydModel_inputs.InputCategories
InputCategoriesUnique = HydModel_inputs.InputCategoriesUnique

#==============================================================================
# FUNCTION RUN 
#==============================================================================

## Initialize dataframes
OutputHydMod_Cat_Rech = pd.DataFrame(columns=InputCategoriesUnique.index,index = range(len(InputClimate)))
OutputHydMod_Cat_Rnff = pd.DataFrame(columns=InputCategoriesUnique.index,index = range(len(InputClimate)))
OutputHydMod_Cat_AET = pd.DataFrame(columns=InputCategoriesUnique.index,index = range(len(InputClimate)))
OutputHydMod_Cat_Deficit = pd.DataFrame(columns=InputCategoriesUnique.index,index = range(len(InputClimate)))

OutputVar_Cat_Kh = pd.DataFrame(columns=InputCategoriesUnique.index,index = range(len(InputClimate)))
OutputVar_Cat_hUnsat = pd.DataFrame(columns=InputCategoriesUnique.index,index = range(len(InputClimate)))
OutputVar_Cat_Theta = pd.DataFrame(columns=InputCategoriesUnique.index,index = range(len(InputClimate)))

for cat in InputCategoriesUnique.index:
    
    ## Select parameters for specific category
    
    SoilIndex = SoilTypesTable.index(InputCategoriesUnique['SoilClass'][cat])

    SoilThick = SoilThickTab[SoilIndex]
    if (InputCategoriesUnique['LandUse'][cat]=='PaddyK') | (InputCategoriesUnique['LandUse'][cat]=='PaddyR'):
        Thresh = ThreshTab[0]
        Ks = KsPaddy
        ThetaS = ThetaSPaddy
        Lambda = LambdaPaddy
        hbc = hbcPaddy
        Eta = EtaPaddy
    else:
        Thresh = ThreshTab[1]
        if  (InputCategoriesUnique['LandUse'][cat] == 'Forest'):   
            Ks = KsTabForest[SoilIndex]
            ThetaS = ThetaSTabForest[SoilIndex]
            Lambda = LambdaTabForest[SoilIndex]
            hbc = hbcTabForest[SoilIndex]
            Eta = EtaTabForest[SoilIndex]
        elif (InputCategoriesUnique['LandUse'][cat] == 'Scrub'):
            Ks = KsTabScrub[SoilIndex]
            ThetaS = ThetaSTabScrub[SoilIndex]
            Lambda = LambdaTabScrub[SoilIndex]
            hbc = hbcTabScrub[SoilIndex]
            Eta = EtaTabScrub[SoilIndex]
        else:
            Ks = KsTab[SoilIndex]
            ThetaS = ThetaSTab[SoilIndex]
            Lambda = LambdaTab[SoilIndex]
            hbc = hbcTab[SoilIndex]
            Eta = EtaTab[SoilIndex]
    
    ## Get input data for specific category
    
    InputData = pd.DataFrame(columns=['R','PG','RET'])
    InputData['R'] = InputClimate['Rainfall_mm']/1000
    InputData['PG'] = InputPG[cat]
    InputData['RET'] = InputRET[cat]
    InputData.index = pd.to_datetime(InputData.index) 
    
    ## Run model
    
    HydModInter = HydMod_fun(Ks,SoilThick,ThetaS,Thresh,Lambda,hbc,Eta,InputData) # Hydraulic model run for each category
    
    ## Store data
    
    #Fluxes
    OutputHydMod_Cat_Rech[cat] = HydModInter['qout']
    OutputHydMod_Cat_Rnff[cat] = HydModInter['Runoff']
    OutputHydMod_Cat_Deficit[cat] = HydModInter['Deficit']
    OutputHydMod_Cat_AET[cat] = HydModInter['AET']
    
    #Variables
    OutputVar_Cat_Kh[cat] = HydModInter['Kh']
    OutputVar_Cat_hUnsat[cat] = HydModInter['hUnsat']
    OutputVar_Cat_Theta[cat] = HydModInter['Theta1']
    
#==============================================================================
# CLEAN UP DATA
#==============================================================================

## Recharge

OutputHydMod_Cat_Rech = OutputHydMod_Cat_Rech * 1000 # Recharge is converted to mm
OutputHydMod_Cat_Rech[OutputHydMod_Cat_Rech<1] = 0 # Remove negligible values
OutputHydMod_Cat_Rech.index = pd.to_datetime(InputClimate.index) # Indexes are set


## Runoff

OutputHydMod_Cat_Rnff = OutputHydMod_Cat_Rnff * 1000 # Runoff is converted to mm
OutputHydMod_Cat_Rnff.index = pd.to_datetime(InputClimate.index) # Indexes are set
OutputHydMod_Cat_Rnff_Yrly = OutputHydMod_Cat_Rnff.resample("A").sum() # Resample at yearly interval
OutputHydMod_Cat_Rnff_Yrly.index = OutputHydMod_Cat_Rnff_Yrly.index.year # Set indexes as years

## Deficit

OutputHydMod_Cat_Deficit = OutputHydMod_Cat_Deficit * 1000 # Deficit is converted to mm
OutputHydMod_Cat_Deficit.index = pd.to_datetime(InputClimate.index) # Indexes are set

## AET

OutputHydMod_Cat_AET.index = pd.to_datetime(InputClimate.index) # Indexes are set
OutputHydMod_Cat_AET_Yrly = OutputHydMod_Cat_AET.resample("A").sum() # Resample at _Yrly interval
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

#Hydraulic conductivity
OutputVar_Cat_Kh_mask = ma.masked_array(OutputVar_Cat_Kh, mask=OutputVar_Cat_Kh<1e-5) # Remove negligible values
OutputVar_Cat_Kh_mean = pd.DataFrame(hmean(OutputVar_Cat_Kh_mask,axis=0)).transpose() # Harmonic mean because variables exist across several orders of magnitude
OutputVar_Cat_Kh_mean.columns = InputCategoriesUnique.index
#Hydraulic pressure
OutputVar_Cat_hUnsat_mask = ma.masked_array(OutputVar_Cat_hUnsat, mask=OutputVar_Cat_hUnsat<-1e+3) # Remove negligible values
OutputVar_Cat_hUnsat_mean = pd.DataFrame(OutputVar_Cat_hUnsat_mask.mean(axis=0)).transpose() # Arithmetic mean 
OutputVar_Cat_hUnsat_mean.columns = InputCategoriesUnique.index
#Soil moisture
OutputVar_Cat_Theta_mean = pd.DataFrame(OutputVar_Cat_Theta.mean(axis=0)).transpose()

#==============================================================================
# SPATIALIZE DATA 
#==============================================================================

## Initialize dataframes

OutputHydMod_GridRaw_Rech = pd.DataFrame(index=InputCategories.index, columns=OutputHydMod_Cat_RechCorr_Yrly.index)
OutputHydMod_GridRaw_Rnff = pd.DataFrame(index=InputCategories.index, columns=OutputHydMod_Cat_Rnff_Yrly.index)
OutputHydMod_GridRaw_AET = pd.DataFrame(index=InputCategories.index, columns=OutputHydMod_Cat_AET_Yrly.index)
OutputHydMod_GridRaw_PG = pd.DataFrame(index=InputCategories.index, columns=OutputHydMod_Cat_Rnff_Yrly.index)

OutputVar_GridRaw_Kh = pd.DataFrame(index=InputCategories.index, columns=["Kh"])
OutputVar_GridRaw_hUnsat = pd.DataFrame(index=InputCategories.index, columns=["hUnsat"])
OutputVar_GridRaw_Theta = pd.DataFrame(index=InputCategories.index, columns=["Theta"])

## Obtain appropriate recharge for soil and land use combination at each grid cell

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
    MatchingDataVarKh = OutputVar_Cat_Kh_mean[MatchingIndex]
    MatchingDataVarhUnsat = OutputVar_Cat_hUnsat_mean[MatchingIndex]
    MatchingDataVarTheta = OutputVar_Cat_Theta_mean[MatchingIndex]
    
    OutputVar_GridRaw_Kh.loc[j] = np.asarray(MatchingDataVarKh)
    OutputVar_GridRaw_hUnsat.loc[j] = np.asarray(MatchingDataVarhUnsat)
    OutputVar_GridRaw_Theta.loc[j] = np.asarray(MatchingDataVarTheta)

## Merge with coordinates

#Fluxes
OutputHydMod_Grid_Rech = InputCoord.join(OutputHydMod_GridRaw_Rech)
OutputHydMod_Grid_Rnff = InputCoord.join(OutputHydMod_GridRaw_Rnff)

#Variables
OutputVar_Grid_Kh = InputCoord.join(OutputVar_GridRaw_Kh)
OutputVar_Grid_hUnsat = InputCoord.join(OutputVar_GridRaw_hUnsat)
OutputVar_Grid_Theta = InputCoord.join(OutputVar_GridRaw_Theta)

##==============================================================================
## INPUT OUTPUT BALANCE
##==============================================================================

# Input
BalancePG_Yrly = OutputHydMod_GridRaw_PG.mean()*1000.0
BalanceR_Yrly = InputData['R'].resample("A").sum()*1000
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

#==============================================================================
# SAVE OUTPUTS
#==============================================================================

OutputHydMod_Mean_Rnff = pd.DataFrame(OutputHydMod_Cat_Rnff.mean(axis=1))

os.chdir(OutputDir) # Sets working directory

## Excel file for data processing
#Fluxes
writer = pd.ExcelWriter(ScenName + "_OUTPUT_HydModel.xlsx") # Excel file will all data
OutputHydMod_Grid_Rech.to_excel(writer,'HydMod_Yrly', index=True)
OutputHydMod_Cat_Rech.to_excel(writer,'HydModPerCat', index=True)
OutputHydMod_Mean_Rnff.to_excel(writer,'MeanRnffTimeSeries', index=True)
InputCategoriesUnique.to_excel(writer, 'Categories', index=True,header='False')
BalanceTab.to_excel(writer,'BalanceRaw', index=True)
writer.save()

#Variables (yearly)
writer = pd.ExcelWriter(ScenName + "_OUTPUT_VariablesMean.xlsx") # Excel file will all data
OutputVar_Grid_Kh.to_excel(writer,'Kh', index=True)
OutputVar_Grid_hUnsat.to_excel(writer,'hUnsat', index=True)
OutputVar_Grid_Theta.to_excel(writer, 'Theta', index=True,header='False')
writer.save()
#Variables (daily)
writer = pd.ExcelWriter(ScenName + "_OUTPUT_VariablesDaily.xlsx") # Excel file will all data
OutputVar_Cat_Kh.to_excel(writer,'Kh', index=True)
OutputVar_Cat_hUnsat.to_excel(writer,'hUnsat', index=True)
OutputVar_Cat_Theta.to_excel(writer, 'Theta', index=True,header='False')
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
    os.chdir(r"C:\Users\Madeleine\Desktop\MAHESHWARAM\01 - SOIL MOISTURE MODEL\01 - MODELS\02 - HydMod\02 - Figures")
    
    ## RECHARGE AND RUNOFF
    
    for yr in OutputHydMod_GridRaw_Rech.columns[2:]:
          
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
        plt.scatter(OutputHydMod_Grid_Rech['X'],OutputHydMod_Grid_Rech['Y'],c=Val, cmap='Blues',s=100,vmin=-0,vmax=+800,marker='s',edgecolor='face') # Default is jet
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
        plt.scatter(OutputHydMod_Grid_Rnff['X'],OutputHydMod_Grid_Rnff['Y'],c=Val, cmap='spring',s=100,vmin=-0,vmax=+1000,marker='s',edgecolor='face') # Default is jet
        plt.colorbar()
        
        plt.savefig(str(yr) + '_' + ScenName + '_HydMod_Runoff.png', dpi = 900)
        plt.show()
      
#        ## HYDRAULIC CONDUCTIVITY (Kh)                  
#        plt.figure('Kh')
#        plt.suptitle("Kh (m/s)" + " - " + ScenName, fontsize=18)
#        Val = OutputVar_Grid_Kh['Kh']/86400
#        listx=[]
#        listy=[]
#        test = shapefile.Reader(shpFilePath)
#        for sr in test.shapeRecords():
#            for xNew,yNew in sr.shape.points:
#                listx.append(xNew)
#                listy.append(yNew)
#        plt.gca().set_aspect('equal', adjustable='box')
#        plt.plot(listx,listy,color='k')
#        plt.xlim((min(listx)-400,max(listx)+400))
#        plt.ylim((min(listy)-400,max(listy)+400))
#        plt.scatter(OutputVar_Grid_Kh['X'],OutputVar_Grid_Kh['Y'],c=Val, cmap='BrBG',vmin=1e-10,vmax=1e-8, s=100,marker='s',edgecolor='face',norm=colors.LogNorm()) # Default is jet
#        plt.colorbar()
#        plt.savefig('Kh' + '_' + ScenName + '_HydMod.png', dpi = 900)
#        plt.show()
#        
#        ## HYDRAULIC PRESSURE (hUnsat)
#        plt.figure('hUnsat')
#        plt.suptitle("hUnsat (cm)" + " - " + ScenName, fontsize=18)
#        Val = - OutputVar_Grid_hUnsat['hUnsat']
#        listx=[]
#        listy=[]
#        test = shapefile.Reader(shpFilePath)
#        for sr in test.shapeRecords():
#            for xNew,yNew in sr.shape.points:
#                listx.append(xNew)
#                listy.append(yNew)
#        plt.gca().set_aspect('equal', adjustable='box')
#        plt.plot(listx,listy,color='k')
#        plt.xlim((min(listx)-400,max(listx)+400))
#        plt.ylim((min(listy)-400,max(listy)+400))
#        plt.scatter(OutputVar_Grid_hUnsat['X'],OutputVar_Grid_hUnsat['Y'],c=Val, cmap='PRGn',vmin=Val.min(),vmax=Val.max(), s=100,marker='s',edgecolor='face') # Default is jet
#        plt.colorbar()
#        plt.savefig('hUnsat' + '_' + ScenName + '_HydMod.png', dpi = 900)
#        plt.show()
#        
#        ## SOIL MOISTURE (Theta)
#        plt.figure('Theta')
#        plt.suptitle("Theta (m3/m3)" + " - " + ScenName, fontsize=18)
#        Val = OutputVar_Grid_Theta['Theta']
#        listx=[]
#        listy=[]
#        test = shapefile.Reader(shpFilePath)
#        for sr in test.shapeRecords():
#            for xNew,yNew in sr.shape.points:
#                listx.append(xNew)
#                listy.append(yNew)
#        plt.gca().set_aspect('equal', adjustable='box')
#        plt.plot(listx,listy,color='k')
#        plt.xlim((min(listx)-400,max(listx)+400))
#        plt.ylim((min(listy)-400,max(listy)+400))
#        plt.scatter(OutputVar_Grid_Theta['X'],OutputVar_Grid_Theta['Y'],c=Val, cmap='winter_r',vmin=0,vmax=0.5, s=100,marker='s',edgecolor='face') # Default is jet
#        plt.colorbar()
#        plt.savefig('Theta' + '_' + ScenName + '_HydMod.png', dpi = 900)
#        plt.show()
        
        ## REFERENCE RECHARGE
        
        #for yr in InputRech_Yrly[InputRech_Yrly.columns[3:]].columns:
        #    
        #    Val = InputRech_Yrly[yr]
        #        
        #    plt.figure('RefRech' + str(yr)) #figsize=(50, 50)
        #    plt.suptitle(yr, fontsize=18)
        #    
        #    listx=[]
        #    listy=[]
        #    test = shapefile.Reader(shpFilePath)
        #    for sr in test.shapeRecords():
        #        for xNew,yNew in sr.shape.points:
        #            listx.append(xNew)
        #            listy.append(yNew)
        #    plt.gca().set_aspect('equal', adjustable='box')
        #    plt.plot(listx,listy,color='k')
        #    
        #    plt.title("Ref recharge")
        #    plt.xlim((min(listx)-400,max(listx)+400))
        #    plt.ylim((min(listy)-400,max(listy)+400))
        #    plt.scatter(InputCoordRef['X'],InputCoordRef['Y'],c=Val, cmap='Blues',s=686.66,vmin=-0,vmax=+500,marker='s',edgecolor='face') # Default is jet
        #    plt.colorbar()
        #    
        #    plt.savefig(str(yr) + '_Ref.png', dpi = 900)
        #    
        #    plt.show()

