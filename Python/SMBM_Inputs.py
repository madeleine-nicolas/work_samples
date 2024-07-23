# -*- coding: utf-8 -*-
"""
Author: Madeleine NICOLAS       
        madeleine.nicolas@univ-rennes1.fr

Purpose: Creates inputs for SMBM for Adil grid

Output:     - Crop coefficient Kc
            - RET (from ET and Kc)
            - Pumping abstraction rates PG
            - Return flow coefficients (Cf)

Requirements:  - Climate data (rainfall and PET)
               - Land Use data (LU percentage per grid cell) for PG, Kc and Cf

Remarks: 

Date:   May 2018
"""

import os
import pandas as pd
import datetime
import matplotlib.pyplot as plt
import shapefile

os.chdir(r"C:\Users\Madeleine\Desktop\Soil moisture model\02 - INPUT DATA") # Sets working directory

########## INPUT DATA ##########

# Land use percentages per grid cell
InputLandUseType = pd.read_csv("INPUT_landuse-type.txt",sep='\t',header=0,index_col=0)  

# Input rainfall and PET
InputClimate = pd.read_csv("INPUT_climate_2002-2015.txt",sep='\t',header=0,index_col=0) 
InputClimate.index=pd.to_datetime(InputClimate.index)

# Input recharge
InputRechYearly = pd.read_csv("INPUT_ref-rech_2002-2015.txt",sep='\t',header=0,index_col=0) 
InputRechMean = InputRechYearly.mean(axis=0)
InputRechMean.index = pd.to_datetime(InputRechMean.index).year

# Initialize tables
InputKc = pd.DataFrame(columns=['Scrub', 'Forest', 'Fruit', 'Tank', 'Urban', 'PaddyR', 'PaddyK'],index = pd.to_datetime(InputClimate.index))
InputPG = pd.DataFrame(columns=['Scrub', 'Forest', 'Fruit', 'Tank', 'Urban', 'PaddyR', 'PaddyK'],index = pd.to_datetime(InputClimate.index))
InputCf = pd.DataFrame(columns=['Scrub', 'Forest', 'Fruit', 'Tank', 'Urban', 'PaddyR', 'PaddyK'],index = pd.to_datetime(InputClimate.index))

# Coordinates
InputCoord = pd.read_csv(r"C:\Users\Madeleine\Desktop\Soil moisture model\02 - INPUT DATA\INPUT_coord.txt",sep='\t',header=0,index_col=2)

# PUMPING (mm)

InputPG['Scrub'] = 0
InputPG['Tank'] = 0
InputPG['Forest'] = 0
InputPG['Fruit'] = 0
InputPG['Urban'] = 0
InputPG['Vegetable'] = 7.7

for yr in range(min(InputPG.index.year),max(InputPG.index.year)+1):
    
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 6, 1))*(InputPG.index.date<=datetime.date(yr, 6, 30)),'PaddyR'] = 10.1*0.1 # Nursery (Kharif)     
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 7, 1))*(InputPG.index.date<=datetime.date(yr, 9, 30)),'PaddyR'] = 10.1 # Irrigation (Kharif) 
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 10, 1))*(InputPG.index.date<=datetime.date(yr, 10, 15)),'PaddyR'] = 0 # Harvesting (Kharif) 
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 10, 16))*(InputPG.index.date<=datetime.date(yr, 12, 15)),'PaddyR'] = 15.2*0.1 # Maintenance + Nursery
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 12, 15))*(InputPG.index.date<=datetime.date(yr, 12, 31)),'PaddyR'] = 15.2 # Irrigation (Rabi)
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 1, 1))*(InputPG.index.date<=datetime.date(yr, 4, 15)),'PaddyR'] = 15.2 # Irrigation (Rabi)
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 4, 16))*(InputPG.index.date<=datetime.date(yr, 4, 30)),'PaddyR'] = 0 # Harvesting (Rabi)
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 5, 1))*(InputPG.index.date<=datetime.date(yr, 5, 31)),'PaddyR'] = 15.2*0.1 # Maintenance + Nursery (Rabi)
                
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 6, 1))*(InputPG.index.date<=datetime.date(yr, 6, 30)),'PaddyK'] = 10.1*0.1 # Nursery (Kharif)     
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 7, 1))*(InputPG.index.date<=datetime.date(yr, 9, 30)),'PaddyK'] = 10.1 # Irrigation (Kharif) 
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 10, 1))*(InputPG.index.date<=datetime.date(yr, 12, 31)),'PaddyK'] = 0 # Harvesting (Kharif) + No pumping
    InputPG.at[(InputPG.index.date>=datetime.date(yr, 1, 1))*(InputPG.index.date<=datetime.date(yr, 5, 31)),'PaddyK'] = 0 # No pumping

# CROP COEFFICIENT

InputKc['Scrub'] = 0.7
InputKc['Forest'] = 1.2      
InputKc['Fruit'] = 0.8  
InputKc['Tank'] = 1  
InputKc['Urban'] = 1  
    
for yr in range(min(InputKc.index.year),max(InputKc.index.year)+1):

    InputKc.at[(InputKc.index.date>=datetime.date(yr, 6, 01))*(InputKc.index.date<=datetime.date(yr, 6, 30)),'Vegetable'] = 0.95
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 7, 01))*(InputKc.index.date<=datetime.date(yr, 7, 31)),'Vegetable'] = 0.7
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 8, 01))*(InputKc.index.date<=datetime.date(yr, 11, 25)),'Vegetable'] = 1.05
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 11, 26))*(InputKc.index.date<=datetime.date(yr, 12, 31)),'Vegetable'] = 0.95
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 01, 01))*(InputKc.index.date<=datetime.date(yr, 1, 10)),'Vegetable'] = 0.95
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 1, 11))*(InputKc.index.date<=datetime.date(yr, 2, 10)),'Vegetable'] = 0.7       
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 2, 11))*(InputKc.index.date<=datetime.date(yr, 5, 31)),'Vegetable'] = 0.85

    InputKc.at[(InputKc.index.date>=datetime.date(yr, 6, 01))*(InputKc.index.date<=datetime.date(yr, 6, 30)), 'PaddyR'] = 1.05
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 7, 01))*(InputKc.index.date<=datetime.date(yr, 10, 31)), 'PaddyR'] = 1.2
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 11, 01))*(InputKc.index.date<=datetime.date(yr, 11, 15)), 'PaddyR'] = 1
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 11, 16))*(InputKc.index.date<=datetime.date(yr, 12, 15)), 'PaddyR'] = 1.05
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 12, 16))*(InputKc.index.date<=datetime.date(yr, 12, 31)), 'PaddyR'] = 1.2
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 1, 1))*(InputKc.index.date<=datetime.date(yr, 3, 25)), 'PaddyR'] = 1.2      
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 3, 26))*(InputKc.index.date<=datetime.date(yr, 5, 31)), 'PaddyR'] = 1
    
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 6, 01))*(InputKc.index.date<=datetime.date(yr, 6, 30)), 'PaddyK'] = 1.05
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 7, 01))*(InputKc.index.date<=datetime.date(yr, 10, 31)), 'PaddyK'] = 1.2
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 11, 01))*(InputKc.index.date<=datetime.date(yr, 12, 31)), 'PaddyK'] = 1
    InputKc.at[(InputKc.index.date>=datetime.date(yr, 1, 1))*(InputKc.index.date<=datetime.date(yr, 5, 31)), 'PaddyK'] = 1     

# REAL EVAPOTRANSPIRATION

InputRET = InputKc.multiply(InputClimate['PET_mm'],axis="index")

# PERCENTAGE OF RECHARGE EQUAL TO RETURN FLOW

InputCf['Scrub'] = 1
InputCf['Forest'] = 1      
InputCf['Fruit'] = 1  
InputCf['Tank'] = 1  
InputCf['Urban'] = 1
InputCf['Vegetable'] = 0.24

for yr in range(min(InputCf.index.year),max(InputCf.index.year)+1):
    
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 6, 1))*(InputCf.index.date<=datetime.date(yr, 6, 30)), 'PaddyR']  = 0.51 # Nursery (Kharif)     
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 7, 1))*(InputCf.index.date<=datetime.date(yr, 9, 30)), 'PaddyR']  = 0.51 # Irrigation (Kharif) 
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 10, 1))*(InputCf.index.date<=datetime.date(yr, 10, 15)), 'PaddyR']  = 1 # Harvesting (Kharif) 
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 10, 16))*(InputCf.index.date<=datetime.date(yr, 12, 15)), 'PaddyR']  = 0.48 # Maintenance + Nursery
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 12, 15))*(InputCf.index.date<=datetime.date(yr, 12, 31)), 'PaddyR']  = 0.48 # Irrigation (Rabi)
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 1, 1))*(InputCf.index.date<=datetime.date(yr, 4, 15)), 'PaddyR']  = 0.48 # Irrigation (Rabi)
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 4, 16))*(InputCf.index.date<=datetime.date(yr, 4, 30)), 'PaddyR']  = 1 # Harvesting (Rabi)
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 5, 1))*(InputCf.index.date<=datetime.date(yr, 5, 31)), 'PaddyR']  = 0.48 # Maintenance + Nursery (Rabi)
     
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 6, 1))*(InputCf.index.date<=datetime.date(yr, 6, 30)), 'PaddyK']  = 0.51 # Nursery (Kharif)     
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 7, 1))*(InputCf.index.date<=datetime.date(yr, 9, 30)), 'PaddyK']  = 0.51 # Irrigation (Kharif) 
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 10, 1))*(InputCf.index.date<=datetime.date(yr, 12, 31)), 'PaddyK']  = 1 # Harvesting (Kharif) + No pumping
    InputCf.at[(InputCf.index.date>=datetime.date(yr, 1, 1))*(InputCf.index.date<=datetime.date(yr, 5, 31)), 'PaddyK']  = 1 # No pumping

    InputCf.at[(InputCf.index.date>=datetime.date(yr, 6, 1))*(InputCf.index.date<=datetime.date(yr, 10, 31)),'Vegetable'] = 0.26 # Kharif
     
########## WEIGHTED DATA ##########

## Obtain data for each grid cell using land use percentages
    
InputPG_Discr = pd.DataFrame(0,index=InputClimate.index, columns=InputLandUseType.index) # Initialize dataframes
InputRET_Discr = pd.DataFrame(0,index=InputClimate.index, columns=InputLandUseType.index)
InputCf_Discr = pd.DataFrame(0,index=InputClimate.index, columns=InputLandUseType.index)

for gr in InputLandUseType.index:
    for lu in InputLandUseType.columns:
       InputRET_Discr[gr] += InputRET[lu[:-3]]*InputLandUseType[lu][gr]/100
       InputPG_Discr[gr] += InputPG[lu[:-3]]*InputLandUseType[lu][gr]/100
       InputCf_Discr[gr] += InputCf[lu[:-3]]*InputLandUseType[lu][gr]/100





#################
####   PLOT   ###
#################
#
#shpFilePath = r"C:\Users\Madeleine\Desktop\Soil moisture model\02 - INPUT DATA\GIS\watershed.shp"
#
#os.chdir(r"C:\Users\Madeleine\Desktop\Soil moisture model\01 - MODELS\Figures\SMBM\Inputs")
#
### PUMPING 
#       
#ValPG = pd.DataFrame(InputPG_Discr.resample("A").sum())
#ValPG.index = ValPG.index.year
#ValPG = ValPG.transpose()
#ValPG = pd.DataFrame(ValPG[2002][ValPG[2002].index[ValPG[2002].index.isin(InputRechYearly.index)]]) # Get matching indexes, All years are the same
#ValPG.columns = ['PG']
#
#plt.figure('Pumping') #figsize=(50, 50)
#
#listx=[]
#listy=[]
#test = shapefile.Reader(shpFilePath)
#for sr in test.shapeRecords():
#    for xNew,yNew in sr.shape.points:
#        listx.append(xNew)
#        listy.append(yNew)
#plt.gca().set_aspect('equal', adjustable='box')
#plt.plot(listx,listy,color='k')
#
#plt.title("Pumping")
#plt.xlim((min(listx)-400,max(listx)+400))
#plt.ylim((min(listy)-400,max(listy)+400))
#
#plt.scatter(InputCoord['X'],InputCoord['Y'],c=ValPG, cmap='Blues',s=686.66,vmin=-0,vmax=+2500,marker='s',edgecolor='face') # Default is jet
#plt.colorbar()
#
#plt.savefig('Pumping.png',dpi = 900)
#plt.show()
#
### REAL EVAPOTRANSPIRATION
#       
#ValRET = pd.DataFrame(InputRET_Discr.resample("A").sum())
#ValRET.index = ValRET.index.year
#ValRET = ValRET.transpose()
#ValRET = pd.DataFrame(ValRET[2002][ValRET[2002].index[ValRET[2002].index.isin(InputRechYearly.index)]]) #  Get matching indexes, All years are the same
#ValRET.columns = ['RET']
#
#plt.figure('RET') #figsize=(50, 50)
#
#listx=[]
#listy=[]
#test = shapefile.Reader(shpFilePath)
#for sr in test.shapeRecords():
#    for xNew,yNew in sr.shape.points:
#        listx.append(xNew)
#        listy.append(yNew)
#plt.gca().set_aspect('equal', adjustable='box')
#plt.plot(listx,listy,color='k')
#
#plt.title("RET")
#plt.xlim((min(listx)-400,max(listx)+400))
#plt.ylim((min(listy)-400,max(listy)+400))
#
#plt.scatter(InputCoord['X'],InputCoord['Y'],c=ValRET, cmap='autumn_r',s=686.66,vmin=-0,vmax=+3000,marker='s',edgecolor='face') # Default is jet
#plt.colorbar()
#
#plt.savefig('RET.png',dpi = 900)
#plt.show()
#
### REAL EVAPOTRANSPIRATION
#       
#ValCf = pd.DataFrame(InputCf_Discr.resample("A").sum())/pd.DataFrame(InputCf_Discr.resample("A").count()) # mean doesn't work for some reason
#ValCf.index = ValCf.index.year
#ValCf = ValCf.transpose()
#ValCf = pd.DataFrame(ValCf[2002][ValCf[2002].index[ValCf[2002].index.isin(InputRechYearly.index)]]) #  Get matching indexes, All years are the same
#ValCf.columns = ['Cf']
#
#plt.figure('Cf') #figsize=(50, 50)
#
#listx=[]
#listy=[]
#test = shapefile.Reader(shpFilePath)
#for sr in test.shapeRecords():
#    for xNew,yNew in sr.shape.points:
#        listx.append(xNew)
#        listy.append(yNew)
#plt.gca().set_aspect('equal', adjustable='box')
#plt.plot(listx,listy,color='k')
#
#plt.title("Cf")
#plt.xlim((min(listx)-400,max(listx)+400))
#plt.ylim((min(listy)-400,max(listy)+400))
#
#plt.scatter(InputCoord['X'],InputCoord['Y'],c=ValCf, cmap='YlGnBu',s=686.66,vmin=-0,vmax=+1,marker='s',edgecolor='face') # Default is jet
#plt.colorbar()
#
#plt.savefig('Cf.png',dpi = 900)
#plt.show()
#
#
#######################
####   SAVE INPUTS  ###
#######################
#
#ValClimate = pd.DataFrame(InputClimate.resample("A").sum())
#ValClimate.index = ValClimate.index.year
#ValClimate.columns = ['PET','R']
#
#ValTotal= ValCf.join(ValRET).join(ValPG)
#
#os.chdir(r"C:\Users\Madeleine\Desktop\Soil moisture model\01 - MODELS\01 - SMBM\Outputs") # Sets working directory
#
#writer = pd.ExcelWriter("TOTAL" + "_INPUT_SMBM.xlsx")
#ValClimate.to_excel(writer,'General inputs', index=True)
#ValTotal.to_excel(writer,'Discretized inputs', index=True)
#writer.save()