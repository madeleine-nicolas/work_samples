# -*- coding: utf-8 -*-
"""
Author: Madeleine NICOLAS       
        madeleine.nicolas@univ-rennes1.fr

Purpose: Physically based soil moisture model based on Dewandel et al., 2008, test run with same data as 'RF-2 layers.xls'

Output: Recharge rates

Requirements: Input ETR, R and PG data

Remarks: 

Date:   2018
"""

import os
import pandas as pd
import datetime

os.chdir(r"C:\Users\Madeleine\Desktop\MAHESHWARAM\01 - SOIL MOISTURE MODEL\03 - PYTHON CODES") # Sets working directory
from HydModel__def import HydMod_fun

os.chdir(r"C:\Users\Madeleine\Desktop\MAHESHWARAM\01 - SOIL MOISTURE MODEL\02 - INPUT DATA") # Sets working directory

########## INPUT PARAM TABLES ##########

Ks = 8.69581991E-06 # Saturated Hydraulic conductivity (m/s)
Ks = Ks * 86400 # (m/day)
SoilThick = 0.17 # Soil thickness (m)
ThetaS = 0.358 # Saturated soil moisture content (m3/m3)
Thresh =  0.06 # Height difference between top of soil and edge bordering field (m)
Lambda = 0.1166137722219530 # Pore-size-distribution parameter
hbc = -0.2000/100 # Air-entry pressure head (m)
Eta = 20.1506329131809 # Tortuosity parameter 

########## INPUT DATA ##########

InputData = pd.read_csv("INPUT_climate_Benoit.txt",sep='\t',header=0,index_col=0) # Potential evaporation & rainfall (in mm)
InputData['PG'] = 0 # No pumping data
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

TestRun = pd.DataFrame(HydMod_fun(Ks,SoilThick,ThetaS,Thresh,Lambda,hbc,Eta,InputData))
       
########### SAVE DATA ##########

#os.chdir(r"C:\Users\Madeleine\Desktop\Soil moisture model\01 - MODELS\02 - HydMod\Outputs") # Sets working directory
#
#writer = pd.ExcelWriter("TestRun" + "_OUTPUT_HydModel.xlsx")
#TestRun.to_excel(writer,'Hydraulic Model', index=True)
#writer.save()
