# -*- coding: utf-8 -*-
"""
Author: Madeleine NICOLAS       
        madeleine.nicolas@univ-rennes1.fr

Purpose: Estimates sensitivity of HydMod model to all input paramters (Ks, SoilThick, ThetaS, Thresh, Lambda, hbc, Eta)

Output: Sensitivity per parameter

Requirements: HydModel_inputs is run to obtain input climate data

Remarks: 

Date:   July 2018
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time 

t = time.time()

os.chdir(r"C:\Users\Madeleine\Desktop\MAHESHWARAM\01 - SOIL MOISTURE MODEL\03 - PYTHON CODES") # Sets working directory

from HydModel__def import HydMod_fun
import HydModel_FineGrid_inputs

NumIter = 50 # Number of values tested for each parameter

#==============================================================================
# INPUT PARAM TABLES 
#==============================================================================

KsTab = np.multiply(np.logspace(-7,-4,NumIter),86400) # Saturated Hydraulic conductivity (m/day)
SoilThickTab = np.linspace(0.04,3,NumIter) # Soil thickness (m)
ThetaSTab = np.linspace(0.01,1,NumIter) # Saturated soil moisture content (m3/m3)
ThreshTab = np.linspace(0.01,0.2,NumIter) # Height difference between top of soil and edge bordering field (m)
LambdaTab = np.linspace(0.01,2,NumIter) # Pore-size-distribution parameter
hbcTab = -np.linspace(0.01,3,NumIter) # Air-entry pressure head (m) 
EtaTab = np.linspace(1,100,NumIter) # Tortuosity parameter

#########################################
## JOIN PARAMETER TABLES IN DICTIONARY ##
#########################################

dct_tables = {'Ks': KsTab, 'SoilThick': SoilThickTab, 'ThetaS': ThetaSTab, 
              'Thresh': ThreshTab, 'Lambda': LambdaTab, 'hbc': hbcTab, 'Eta': EtaTab}

#==============================================================================
# INPUT DATA 
#==============================================================================

# Input data (Rainfall + RET)
InputData = pd.DataFrame(HydModel_FineGrid_inputs.InputClimate['Rainfall_mm']/1000) # in m
InputData['RET'] = HydModel_FineGrid_inputs.InputRET['B'] # Scrub/Alfisols1
InputData.columns = ['R', 'RET']

# Partition coef between IRF and natural Recharge
InputCf = HydModel_FineGrid_inputs.InputCf 

#==============================================================================
# FUNCTION RUN
# Runs function varying each parameter one at the time (the others are default) 
#==============================================================================

# Reference parameters (Alfisols 1)
dct_ref = {'Ks': 1.85E-5 * 86400, 'SoilThick': 0.20, 'ThetaS': 0.32, 'Thresh': 0.12, 'Lambda': 0.175, 'hbc': -0.089, 'Eta': 14.8}
dct_ref_s = sorted(dct_ref.iteritems()) # Sort dictionary
EtaRef, KsRef, LambdaRef, ThetaSRef, SoilThickRef, ThreshRef, hbcRef = [v[1] for v in dct_ref_s] # Separate dictionary into variables

# Reference output
HydModRef = pd.DataFrame(HydMod_fun(KsRef,SoilThickRef,ThetaSRef,ThreshRef,LambdaRef,hbcRef,EtaRef,InputData)['qout'])

# Initialize dataframe
EmptyTab = pd.DataFrame(columns=range(NumIter),index = HydModRef.index)
HydModOutput = {'Ks': EmptyTab.copy(), 'SoilThick': EmptyTab.copy(), 'ThetaS': EmptyTab.copy(), 'Thresh': EmptyTab.copy(), 'Lambda': EmptyTab.copy(), 'hbc': EmptyTab.copy(), 'Eta': EmptyTab.copy()}

for param in dct_ref.keys():
    
    print(param) 
    
    for i in range(NumIter):
        
        print(i)
        
        # Set parameters to ref, except parameter of interest
        dct_inter = dct_ref.copy()
        dct_inter[param] = dct_tables[param][i] # Change parameter of interest
        dct_inter_s = sorted(dct_inter.iteritems()) # Sort dictionary
        Eta, Ks, Lambda, ThetaS, SoilThick, Thresh, hbc = [v[1] for v in dct_inter_s] # Separate dictionary into variables
        
        HydModInter = pd.DataFrame(HydMod_fun(Ks,SoilThick,ThetaS,Thresh,Lambda,hbc,Eta,InputData)['qout']) # Intermediate step
        HydModOutput[param][i] = HydModInter.copy() # Fill total output table

#==============================================================================
# SENSITIVITY TEST
#==============================================================================

HydModSensitivity = {'Ks': EmptyTab.copy(), 'SoilThick': EmptyTab.copy(), 'ThetaS': EmptyTab.copy(), 'Thresh': EmptyTab.copy(), 'Lambda': EmptyTab.copy(), 'hbc': EmptyTab.copy(), 'Eta': EmptyTab.copy()}

for param in dct_ref.keys():    
    for i in range(NumIter-1):
        HydModSensitivity[param][i]=dct_tables[param][i]*(HydModOutput[param][i+1]-HydModOutput[param][i])/(dct_tables[param][i+1]-dct_tables[param][i])

MeanSensitivity = pd.DataFrame(columns=['Sensitivity'],index=dct_ref.keys())
MedianSensitivity = pd.DataFrame(columns=['Sensitivity'],index=dct_ref.keys())
MaxSensitivity = pd.DataFrame(columns=['Sensitivity'],index=dct_ref.keys())

for param in dct_ref.keys():
    MeanSensitivity['Sensitivity'][param]=abs(HydModSensitivity[param]).mean().mean()
    MedianSensitivity['Sensitivity'][param]=abs(HydModSensitivity[param]).mean().median()
    MaxSensitivity['Sensitivity'][param]=abs(HydModSensitivity[param]).mean().max(skipna=True)
    
MeanSensitivity = MeanSensitivity.sort(columns='Sensitivity')
MedianSensitivity = MedianSensitivity.sort(columns='Sensitivity')
MaxSensitivity = MaxSensitivity.sort(columns='Sensitivity')
#==============================================================================
# PLOT
#==============================================================================

os.chdir(r"C:\Users\Madeleine\Desktop\MAHESHWARAM\01 - SOIL MOISTURE MODEL\01 - MODELS\Figures\Raw")

for param in dct_ref.keys():
    
    ## Normal plot
    plt.figure(param + 'Sensitivity')
    plt.ylabel('Sens')
    plt.xlabel('Param Value')
    plt.title(param,fontsize=18)
    plt.plot(dct_tables[param],abs(HydModSensitivity[param]).mean(),'g-o')    
    plt.savefig(param + '_Sensitivity.png', dpi = 900)
    plt.show()

    ## Side by side plot
    fig, ax1 = plt.subplots()
    ax1.plot(dct_tables[param], HydModOutput[param].sum(), 'b-')
    ax1.set_xlabel('Param Value')
    ax1.set_ylabel('SumRech', color='b')
    ax1.tick_params('y', colors='b')
    ax2 = ax1.twinx()
    ax2.plot(dct_tables[param], abs(HydModSensitivity[param]).mean(), 'r-o')
    ax2.set_ylabel('Sens', color='r')
    ax2.tick_params('y', colors='r')
    plt.title(param,fontsize=18)
    fig.tight_layout()
    plt.savefig(param + '_Sensitivity2.png', dpi = 900)
    plt.show()
        
# Print total elapsed time
elapsed = time.time() - t 
print ('Time elapsed: ' + '%.0f' % (elapsed/60) + ' min ' + '%.0f' % ((elapsed/60-int(elapsed/60))*60) + ' sec')