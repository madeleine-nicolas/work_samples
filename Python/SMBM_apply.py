# -*- coding: utf-8 -*-
"""
Author: Madeleine NICOLAS
        madeleine.nicolas@univ-rennes1.fr

Purpose: Applies SMBM at each grid using actual soil proportions

Output: Modelled yearly recharge at each grid cell and error relative to observed recharge

Requirements: Input climate (P and ET data), Soil type properties and cell percentages

Remarks:

Date:   January 2018
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from SMBM_def import SMBMSoilType_fun

os.chdir(r"C:\Users\Madeleine\Desktop\Soil moisture model\02 - INPUT DATA") # Sets working directory

########## INPUT DATA ##########

# COORDINATES
InputCoord = pd.read_csv("INPUT_coord.txt",sep='\t',header=0,index_col=2) # Coordinates of reference points

# SOIL TYPE PERCENTAGE PER GRID CELL
InputSoilType = pd.read_csv("INPUT_soil-type.txt",sep='\t',header=0,index_col=0)

# SOIL TYPE PROPERTIES
InputSoilAWC = pd.DataFrame(0,index=['Alfi1','Alfi2','Incep','Enti','Tank'],columns=['SoilAWC']) # Initialize Available Water Content (mm)
InputSoilMWC = pd.DataFrame(0,index=['Alfi1','Alfi2','Incep','Enti','Tank'],columns=['SoilMWC']) # Initialize Maximum surface storage (mm)
SMBMSoilType = pd.DataFrame(columns=['Alfi1','Alfi2','Incep','Enti','Tank'])

# REFERENCE RECHARGE
InputRechYearly = pd.read_csv("INPUT_ref-rech_2002-2015.txt",sep='\t',header=0,index_col=0)
InputRechMean = InputRechYearly.mean(axis=0)
InputRechMean.index = pd.to_datetime(InputRechMean.index).year

IWC = 10 # Initial Water Content (mm)

#############################
#####    USE FUNCTION   #####
#############################

AWC = [39,62,45,40,8]
MWC = [18,22,24,49,47]
#AWC = [52,52,52,52,52]
#MWC = [31,31,31,31,31]

InputSoilAWC['SoilAWC'] = AWC
InputSoilMWC['SoilMWC'] = MWC

SMBMDiscrYearly, ErrTab, ErrTabTot = SMBMSoilType_fun(InputSoilAWC,InputSoilMWC,IWC, InputSoilType, SMBMSoilType, InputRechYearly)

print(ErrTabTot)

#####################
#####    PLOT   #####
#####################

## Scatter plot with yearly values averaged for all basin

Obs = pd.DataFrame(InputRechMean)
Obs.columns = ['Obs']
Sim = pd.DataFrame(np.mean(SMBMDiscrYearly,axis=0))
Sim.columns = ['Sim']
TabScatter = pd.merge(Obs, Sim, left_index=True, right_index=True)
#TabScatter = TabScatter[3:]
labels = np.asarray(TabScatter.index)

fig, ax = plt.subplots()
ax.scatter(TabScatter['Obs'],TabScatter['Sim'],s = 200, marker ='o',alpha=0.5,c='c')
for i, txt in enumerate(labels):
    ax.annotate(txt, (np.asarray(TabScatter['Obs'])[i]-5,np.asarray(TabScatter['Sim'])[i]+8))
plt.title('Correlation Ref/SMBM')
plt.axis([0,400,0,400])
plt.gca().set_aspect('equal', adjustable='box')
plt.plot([0, 400], [0, 400], 'k-', lw=2)

############ SAVE DATA ##########
#
#os.chdir("C:\Users\Madeleine\Desktop\Soil moisture model\Data\Outputs") # Sets working directory
#
#InputSoilAWC['Error'] = ErrTabTot
#SMBMDiscrYearly = SMBMDiscrYearly[SMBMDiscrYearly.index.isin(np.asarray(InputCoord.index))]
#
#writer = pd.ExcelWriter("OUTPUT_SMBM.xlsx")
#InputCoord.join(SMBMDiscrYearly).to_excel(writer,'Data', index=True)
#InputCoord.join(ErrTab).to_excel(writer,'Error', index=True)
#InputSoilAWC.to_excel(writer,'Params',index = True)
#writer.save()