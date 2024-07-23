# -*- coding: utf-8 -*-
"""
Author: Madeleine NICOLAS
        madeleine.nicolas@univ-rennes1.fr

Purpose: Calibrates SMBM params (AWC and MWC) at each grid cell; pumping, RET and Cf coefficients are estimated by weighting data
         according to land use surface MINIMIZES RMSE

Output: AWC and MWC parameters at each grid cell

Requirements: SMBM_Inputs

Remarks:

Date:   May 2018
"""

import os
import numpy as np
import pandas as pd
from scipy.optimize import minimize
import time 
import sys
sys.path.append(r"C:\Users\Madeleine\Desktop\Soil moisture model\03 - PYTHON CODES")

import SMBM_Inputs
from SMBM_Def import SMBM_Cell_fun

os.chdir(r"C:\Users\Madeleine\Desktop\Soil moisture model\02 - INPUT DATA") # Sets working directory

######################
###   INPUT DATA   ###
######################

# Coordinates
InputCoord = pd.read_csv(r"C:\Users\Madeleine\Desktop\Soil moisture model\02 - INPUT DATA\INPUT_coord.txt",sep='\t',header=0,index_col=2)

# Inputs (PG + rainfall)
InputRainfall = pd.DataFrame(SMBM_Inputs.InputClimate['Rainfall_mm']) # Rainfall in mm
InputPG = SMBM_Inputs.InputPG_Discr # Pumping in mm, obtained at each grid cell by weighting according to land use
InputTotal = pd.DataFrame(index=InputRainfall.index, columns= InputPG.columns) # Sum of rainfall and pumping in mm
for lu in InputPG.columns:
    InputTotal[lu] = InputRainfall['Rainfall_mm'] + InputPG[lu]

# Evapotranspiration  
InputRET = SMBM_Inputs.InputRET_Discr

# Partition coef between IRF and natural Recharge, obtained at each grid cell by weighting according to land use
InputCf = SMBM_Inputs.InputCf_Discr # Obtained at each grid cell by weighting according to land use

# Reference recharge
InputRechYearly = pd.read_csv("INPUT_ref-rech_2002-2015.txt",sep='\t',header=0,index_col=0).transpose() 
InputRechYearly.index = pd.to_datetime(InputRechYearly.index).year # Inputs into datetime format

# Parameters
IWC = 10 # Initial Water Content (mm)

####################
###   FUNCTION   ###
####################


def SMBM_CellByCell_Opt(x):
    
    """Model info...
    """
    
    AWC = x[0]
    MWC = x[1]
    
    Err = SMBM_Cell_fun(gr, AWC, MWC, IWC, InputTotal, InputRET, InputCf, InputRechYearly)[1]
    
    return Err
    
########################
###   FUNCTION RUN   ###
########################

ErrTab = pd.DataFrame(index=InputRechYearly.columns, columns = ['AWC','MWC'])

x0 = np.array([200,20]) 
bnds = ((5, 600), (5, 600))

for gr in InputRechYearly.columns[0:]:
    
    t = time.time()
    
    res = minimize(SMBM_CellByCell_Opt, x0, method='L-BFGS-B', bounds=bnds, options={'maxiter': 15, 'disp': False, 'maxls':15})
    
    elapsed = time.time() - t
    
    print(res.x)
    print ('Time elapsed: ' + '%.0f' % (elapsed/60) + ' min ' + '%.0f' % ((elapsed/60-int(elapsed/60))*60) + ' sec')
    
    ErrTab.at[gr,'AWC'] = res.x[0]
    ErrTab.at[gr,'MWC'] = res.x[1]
    