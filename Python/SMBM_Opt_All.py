# -*- coding: utf-8 -*-
"""
Author: Madeleine NICOLAS       
        madeleine.nicolas@univ-rennes1.fr

Purpose: Simultaneous parameter optimization for SMBM model

Output: Ideal set of parameters and error 

Requirements: 

Remarks: 

Date:   January 2018
"""

import os
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from SMBM_def import SMBM_fun

os.chdir("C:\Users\Madeleine\Desktop\Soil moisture model\Data\Inputs") # Sets working directory

########## INPUT DATA ##########

# SOIL TYPE PERCENTAGE PER GRID CELL
InputSoilType = pd.read_csv("INPUT_soil-type.txt",sep='\t',header=0,index_col=0)

# SOIL TYPE PROPERTIES
InputSoilAWC = pd.DataFrame(0,index=['Alfi1','Alfi2','Incep','Enti','Tank'],columns=['SoilAWC'])
InputSoilAWC['SoilAWC'] = [80,np.NaN,100,np.NaN,180]

# REFERENCE RECHARGE
InputRechYearly = pd.read_csv("INPUT_ref-rech_2002-2015.txt",sep='\t',header=0,index_col=0) 
InputRechMean = InputRechYearly.mean(axis=0)
InputRechMean.index = pd.to_datetime(InputRechMean.index).year

### INPUT CLIMATE DATA
#InputClimate = pd.read_csv("INPUT_climate_2011-2015.txt",sep='\t',header=0,index_col=0)
#InputClimate.index=pd.to_datetime(InputClimate.index)

IWC = 10 # Initial Water Content (mm)
MWC = 30 # Maximum surface storage (mm)

########## OUT OF LOOP ##########

## Obtain recharge for each soil type THAT DON'T NEED OPTIMIZING
SMBMSoilType = pd.DataFrame(columns=['Alfi1','Alfi2','Incep','Enti','Tank'])
for c in np.array(['Alfi1','Incep','Tank']):
    AWC = InputSoilAWC['SoilAWC'][c]
    SMBMSoilType[c] = SMBM_fun(IWC,MWC,AWC)[0]['Rech']
    print(c)
    
################################
#####    DEFINE FUNCTION   #####
################################

def SMBMSoilTypeOpt_fun(x):
    
    """Parameter optimization for SMBM model
    Inputs: AWC for Alfi2 and Enti soils
    Fixed params: MWC and IWC for all soils, AWC for Alfi1, Incep and Tank
    """
    
    InputSoilAWC['SoilAWC']['Alfi2'] = x[0] # Input parameters, Available Water Contents for Alfi2 and Enti soils
    InputSoilAWC['SoilAWC']['Enti'] = x[1]
    
    ## Obtain recharge for each soil type THAT NEEDS OPTIMIZING
    for c in np.array(['Alfi2','Enti']):
        AWC = InputSoilAWC['SoilAWC'][c]
        SMBMSoilType[c] = SMBM_fun(IWC,MWC,AWC)[0]['Rech']
        print(c)
    
    ## Obtain for each grid cell using soil type percentages    
    SMBMDiscr = pd.DataFrame(0,index=SMBMSoilType.index, columns=InputSoilType.index) # Initialize dataframe
    for gr in InputSoilType.index: # Discretized model
        for st in SMBMSoilType.columns:
           SMBMDiscr[gr] += SMBMSoilType[st]*InputSoilType[st+'Per'][gr]/100

    SMBMDiscrYearly = SMBMDiscr.resample("A").sum().transpose() # Yearly average
    SMBMDiscrYearly.columns = SMBMDiscrYearly.columns.year
    
    ## Error calculation
    ErrTab = pd.DataFrame(0,index=InputRechYearly.index, columns=InputRechYearly.columns)
    
    for y in InputRechYearly.columns:
        #ErrTab[str(y)] = abs((SMBMDiscrYearly[int(y)]-InputRechYearly[str(y)])/InputRechYearly[str(y)])*100
        ErrTab[str(y)] = abs((SMBMDiscrYearly[int(y)]-InputRechYearly[str(y)]))
        
    ErrTabAv = np.mean(ErrTab,axis=0)
    ErrTabTot = np.mean(ErrTabAv,axis=0)    
    
    print('done')
    
    #return (SMBMDiscrYearly,ErrTab,ErrTabTot)
    return (ErrTabTot)

#########################
#####    OPTIMIZE   #####
#########################

x0 = np.array([100,40]) 
bnds = ((90, 140), (20, 60))

res = minimize(SMBMSoilTypeOpt_fun, x0, method='L-BFGS-B', bounds=bnds,
                options={'maxiter': 20, 'disp': True, 'maxls':15})
                
print(res.x)