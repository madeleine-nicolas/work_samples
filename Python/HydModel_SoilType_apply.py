# -*- coding: utf-8 -*-
"""
Author: Madeleine NICOLAS       
        madeleine.nicolas@univ-rennes1.fr

Purpose: Physically based soil moisture model based on Dewandel et al., 2008

Output: Recharge rates

Requirements: Input ETR, R and PG data

Remarks: 

Date:   February 2018
"""

import os
import pandas as pd
import numpy as np
import datetime
import matplotlib.pyplot as plt

from HydModel_def import HydMod_fun

os.chdir(r"C:\Users\Madeleine\Desktop\Soil moisture model\02 - INPUT DATA") # Sets working directory

########## INPUT PARAM TABLES ##########

SoilTypesTable = ['Alfi1','Alfi2','Incep','Enti','Tank']
KsTab = [1.87E-5, 2.48E-5, 1.24E-5, 4.23E-6, 2.5E-7] # Saturated Hydraulic conductivity (m/s)
KsTab = np.multiply(KsTab,86400) # (m/day)
SoilThickTab = [0.20,  0.20, 0.10, 0.20, 0.50] # Soil thickness (m)
SigmaSTab = [0.32, 0.36, 0.35, 0.36, 0.60] # Saturated soil moisture content (m3/m3)
Thresh = 100 # Height difference between top of soil and edge bordering field (m)

# Brooks and Corey parameters
LambdaTab = [0.175, 0.193, 0.129, 0.117, 0.033] # Pore-size-distribution parameter
hbcTab = np.divide([0.089, 0.111, 0.060, 0.037, 0.312],100) # Air-entry pressure head (m)
EtaTab = [14.8, 13.3, 19.3, 20.2, 63.7] # Tortuosity parameter 

########## SELECTED PARAMS ##########

SoilNum = input('Choose soil type:\n1:Alfi1\n2:Alfi2\n3:Incep\n4:Enti\n5:Tank\n')

Ks = KsTab[SoilNum-1]
SoilThick = SoilThickTab[SoilNum-1]
SigmaS = SigmaSTab[SoilNum-1]
Lambda = LambdaTab[SoilNum-1]
hbc = hbcTab[SoilNum-1]
Eta = EtaTab[SoilNum-1]
         
########## INPUT DATA ##########

InputData = pd.read_csv("INPUT_climate_2002-2015.txt",sep='\t',header=0,index_col=0) # Potential evaporation & rainfall (in mm)
InputData['PG'] = 0 # No pumping data so far
InputData.columns = [["PET","R","PG"]]
InputData.index=pd.to_datetime(InputData.index)
InputData = InputData/1000 # (in m)

InputData['Kc'] = 0.0
for yr in range(min(InputData.index.year),max(InputData.index.year)):
    InputData.at[(InputData.index.date>=datetime.date(yr, 06, 15))*(InputData.index.date<=datetime.date(yr, 07, 15)),'Kc'] = 0.4
    InputData.at[(InputData.index.date>=datetime.date(yr, 07, 16))*(InputData.index.date<=datetime.date(yr, 11, 25)),'Kc'] = 0.9        
    InputData.at[(InputData.index.date>=datetime.date(yr, 11, 26))*(InputData.index.date<=datetime.date(yr, 12, 15)),'Kc'] = 0.85

InputData['RET']=InputData['PET']
InputData['RET'][InputData['Kc']>0] = InputData['Kc'] * InputData['PET'] * 0.9

########## HYDRAULIC MODEL ##########

OutputData = HydMod_fun(Ks,SoilThick,SigmaS,Thresh,Lambda,hbc,Eta,InputData)
       
########### SAVE DATA ##########

os.chdir(r"C:\Users\Madeleine\Desktop\Soil moisture model\01 - MODELS\02 - HydMod\Outputs") # Sets working directory

writer = pd.ExcelWriter(SoilTypesTable[SoilNum-1] + "_OUTPUT_HydModel.xlsx")
OutputData.to_excel(writer,'Hydraulic Model', index=True)
writer.save()

########### PLOT ##########   

plt.plot(OutputData['Date'], OutputData['qout']*1000,'k', OutputData['Date'], OutputData['qout']*1000,'^',ms=10)
plt.ylim([0,300])
  