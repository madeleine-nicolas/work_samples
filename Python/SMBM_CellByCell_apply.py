# -*- coding: utf-8 -*-
"""
Author: Madeleine NICOLAS
        madeleine.nicolas@univ-rennes1.fr

Purpose: Applies SMBM at each grid using AWC and MWC obtained from model optimisation SMBM_Opt_CellByCell

Output: Modelled daily and yearly recharge and runoff at each grid cell and error relative to observed recharge (yearly only)

Requirements: Inputs from SMBM_Inputs

Remarks:

Date:   June 2018
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import shapefile
import numpy as np
import sys
sys.path.append(r"C:\Users\Madeleine\Desktop\Soil moisture model\03 - PYTHON CODES")

import SMBM_Inputs
from SMBM_Def import SMBM_Cell_fun

os.chdir(r"C:\Users\Madeleine\Desktop\Soil moisture model\02 - INPUT DATA") # Sets working directory

#==============================================================================
# INPUT DATA 
#==============================================================================

# Coordinates
InputCoord = pd.read_csv(r"C:\Users\Madeleine\Desktop\Soil moisture model\02 - INPUT DATA\INPUT_coord.txt",sep='\t',header=0,index_col=2)

# Inputs (PG + rainfall)
InputRainfall = pd.DataFrame(SMBM_Inputs.InputClimate['Rainfall_mm']) # Rainfall in mm
InputRainfall.index = pd.to_datetime(InputRainfall.index.strftime('%Y-%m-%d')) # Format indexes to remove hh:mm:ss but keep datetime format
InputPG = SMBM_Inputs.InputPG_Discr # Pumping in mm, obtained at each grid cell by weighting according to land use
InputTotal = pd.DataFrame(index=InputRainfall.index, columns= InputPG.columns) # Sum of rainfall and pumping in mm
for lu in InputPG.columns:
    InputTotal[lu] = InputRainfall['Rainfall_mm'] + InputPG[lu]

# Evapotranspiration
InputRET = SMBM_Inputs.InputRET_Discr # meters

# Partition coef between IRF and natural Recharge, obtained at each grid cell by weighting according to land use
# Obtained at each grid cell by weighting according to land use
InputCf = SMBM_Inputs.InputCf_Discr # meters

# Reference recharge
InputRechYearly = pd.read_csv("INPUT_ref-rech_2002-2015.txt",sep='\t',header=0,index_col=0).transpose() # meters
InputRechYearly.index = pd.to_datetime(InputRechYearly.index).year # Inputs into datetime format

# Parameters
InputParams = pd.read_csv("INPUT_params_cellbycell.txt",sep='\t',header=0,index_col=0).transpose()  # Available Water Content and Maximum Water Content (mm)
IWC = 10 # Initial Water Content (mm)

#==============================================================================
# FUNCTION RUN 
#==============================================================================

while(True):
    Cond = input("\nDo you want to: \n\n\t1. Run the model from scratch?\n\t2. Upload the last run?\n\n(Type 1 or 2)\n\n")
    if (Cond== 1) | (Cond==2):
        # Break out of loop and continue code 
        break
    else: # Error message for invalid input
        print "\nERROR: Invalid Input"
        exit(0) # Or "continue" to rerun loop, and ask user to re-enter data
    
if Cond == 1: #Run the model from scratch
    
    print("\n\nRunning model...\n\n")
    # Initialize tables
    SMBM_Rech_Nat = pd.DataFrame(index=InputTotal.index, columns= InputTotal.columns)
    SMBM_Rech_Tot = pd.DataFrame(index=InputTotal.index, columns= InputTotal.columns)
    SMBM_Runoff = pd.DataFrame(index=InputTotal.index, columns= InputTotal.columns)
    SMBM_ErrTab = pd.DataFrame(index=InputTotal.columns,columns = ['Err'])
    
    for gr in SMBM_Rech_Nat.columns[SMBM_Rech_Nat.columns.isin(InputParams.columns)]:
        
        AWC = InputParams[gr]['AWC']
        MWC = InputParams[gr]['MWC']
        
        SMBM_Cell_Inter = SMBM_Cell_fun(gr, AWC, MWC, IWC, InputTotal, InputRET, InputCf, InputRechYearly)
        SMBM_Rech_Nat[gr] = SMBM_Cell_Inter[0]['NatRech']
        SMBM_Rech_Tot[gr] = SMBM_Cell_Inter[0]['TotRech']
        SMBM_Runoff[gr] = SMBM_Cell_Inter[0]['Runoff']
        
        SMBM_ErrTab['Err'][gr] = SMBM_Cell_Inter[1]
        
        print(gr)
    
    # Remove cells for which no params are available        
    SMBM_Rech_Nat = SMBM_Rech_Nat[SMBM_Rech_Nat.columns[SMBM_Rech_Nat.columns.isin(InputParams.columns)]]
    SMBM_Rech_Tot = SMBM_Rech_Tot[SMBM_Rech_Tot.columns[SMBM_Rech_Tot.columns.isin(InputParams.columns)]]
    SMBM_Runoff = SMBM_Runoff[SMBM_Runoff.columns[SMBM_Runoff.columns.isin(InputParams.columns)]]
    
    # Yearly averages
    SMBM_Rech_Nat_Yearly = pd.DataFrame(SMBM_Rech_Nat.resample("A").sum())
    SMBM_Rech_Nat_Yearly.index = SMBM_Rech_Nat_Yearly.index.year
    SMBM_Rech_Nat_Yearly = SMBM_Rech_Nat_Yearly.transpose()
    
    SMBM_Rech_Tot_Yearly = pd.DataFrame(SMBM_Rech_Tot.resample("A").sum())
    SMBM_Rech_Tot_Yearly.index = SMBM_Rech_Tot_Yearly.index.year
    SMBM_Rech_Tot_Yearly = SMBM_Rech_Tot_Yearly.transpose()
    
    SMBM_Runoff_Yearly = pd.DataFrame(SMBM_Runoff.resample("A").sum())
    SMBM_Runoff_Yearly.index = SMBM_Runoff_Yearly.index.year
    SMBM_Runoff_Yearly = SMBM_Runoff_Yearly.transpose()
    
    #######################
    ###   SAVE OUTPUTS  ###
    #######################
    
    os.chdir(r"C:\Users\Madeleine\Desktop\Soil moisture model\01 - MODELS\01 - SMBM\Outputs") # Sets working directory
    
    writer = pd.ExcelWriter("TOTAL" + "_OUTPUT_SMBM.xlsx")
    SMBM_Rech_Nat.to_excel(writer,'SMBMDailyRechNat', index=True)
    SMBM_Rech_Nat_Yearly.to_excel(writer,'SMBMYearlyRechNat', index=True)
    SMBM_Rech_Tot.to_excel(writer,'SMBMDailyRechTot', index=True)
    SMBM_Rech_Tot_Yearly.to_excel(writer,'SMBMYearlyRechTot', index=True)
    SMBM_Runoff.to_excel(writer,'SMBMDailyRunoff', index=True)
    SMBM_Runoff_Yearly.to_excel(writer,'SMBMYearlyRunoff', index=True)
    SMBM_ErrTab.to_excel(writer,'SMBMError', index=True)
    InputRechYearly.transpose().to_excel(writer, 'RefRech', index=True,header='False')
    writer.save()

else: # Loading past model run

    print("\n\nLoading data...\n\n")
    
    os.chdir(r"C:\Users\Madeleine\Desktop\Soil moisture model\01 - MODELS\01 - SMBM\Outputs") # Sets working directory
    
    SMBM_Rech_Nat = pd.read_excel("TOTAL_OUTPUT_SMBM.xlsx",sheetname='SMBMDailyRechNat',sep='\t',header=0,index_col=0)
    SMBM_Rech_Nat_Yearly = pd.read_excel("TOTAL_OUTPUT_SMBM.xlsx",sheetname='SMBMYearlyRechNat',sep='\t',header=0,index_col=0)
    SMBM_Runoff = pd.read_excel("TOTAL_OUTPUT_SMBM.xlsx",sheetname='SMBMDailyRunoff',sep='\t',header=0,index_col=0)
    SMBM_Runoff_Yearly = pd.read_excel("TOTAL_OUTPUT_SMBM.xlsx",sheetname='SMBMYearlyRunoff',sep='\t',header=0,index_col=0)
    SMBM_ErrTab = pd.read_excel("TOTAL_OUTPUT_SMBM.xlsx",sheetname='SMBMError',sep='\t',header=0,index_col=0)
    InputRechYearly = pd.read_excel("TOTAL_OUTPUT_SMBM.xlsx",sheetname='RefRech',sep='\t',header=0,index_col=0)
    
#==============================================================================
# PLOT
#==============================================================================
    
while(True):
    Cond2 = raw_input("\nDo you want to plot figures?\n\n(Type y or n)\n\n")
    if (Cond2 == 'y') | (Cond2 == 'n'):
        # Break out of loop and continue code 
        break
    else: # Error message for invalid input
        print "\nERROR: Invalid Input"
        continue # Or "continue" to rerun loop, and ask user to re-enter data

if Cond2 == 'y':
    
    print("\n\nPlotting figures...\n\n")
    
    shpFilePath = r"C:\Users\Madeleine\Desktop\Soil moisture model\02 - INPUT DATA\GIS\watershed.shp"
    
    for yr in SMBM_Rech_Nat_Yearly.columns:
        
        ## RECHARGE
        
        os.chdir(r"C:\Users\Madeleine\Desktop\Soil moisture model\01 - MODELS\Figures\SMBM\Recharge")
        Val = SMBM_Rech_Nat_Yearly[yr]
        
        plt.figure(str(yr) + '_' +'Recharge') #figsize=(50, 50)
        plt.suptitle(yr, fontsize=18)
        
        listx=[]
        listy=[]
        test = shapefile.Reader(shpFilePath)
        for sr in test.shapeRecords():
            for xNew,yNew in sr.shape.points:
                listx.append(xNew)
                listy.append(yNew)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.plot(listx,listy,color='k')
        
        plt.title("SMBM Recharge")
        plt.xlim((min(listx)-400,max(listx)+400))
        plt.ylim((min(listy)-400,max(listy)+400))
        
        plt.scatter(InputCoord['X'],InputCoord['Y'],c=Val, cmap='Blues',s=686.66,vmin=-0,vmax=+500,marker='s',edgecolor='face') # Default is jet
        plt.colorbar()
        
        plt.savefig(str(yr) + '_Recharge-SMBM.png',dpi = 900)
        plt.show()
        
    #    ## RUNOFF
    #    
    #    os.chdir(r"C:\Users\Madeleine\Desktop\Soil moisture model\01 - MODELS\Figures\SMBM\Runoff")
    #    Val2 = SMBM_Runoff_Yearly[yr]
    #    
    #    plt.figure(str(yr) + '_' +'Runoff') #figsize=(50, 50)
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
    #    plt.title("SMBM Runoff")
    #    plt.xlim((min(listx)-400,max(listx)+400))
    #    plt.ylim((min(listy)-400,max(listy)+400))
    #    
    #    plt.scatter(InputCoord['X'],InputCoord['Y'],c=Val2, cmap='Greens',s=686.66,vmin=-0,vmax=+500,marker='s',edgecolor='face') # Default is jet
    #    
    #    plt.colorbar()
    #    
    #    plt.savefig(str(yr) + '_Runoff-SMBM.png',dpi = 900)
    #    
    #    plt.show()
    
    os.chdir(r"C:\Users\Madeleine\Desktop\Soil moisture model\01 - MODELS\Figures\SMBM")
    
    ## AVAILABLE WATER CONTENT (AWC)
    
    InputParamsPlot = InputParams.transpose()
    
    Val3 = InputParamsPlot['AWC']
    Val3[np.isnan(Val3)] = -999
    
    plt.figure('AWC') #figsize=(50, 50)
    plt.suptitle("AWC", fontsize=18)
    
    listx=[]
    listy=[]
    test = shapefile.Reader(shpFilePath)
    for sr in test.shapeRecords():
        for xNew,yNew in sr.shape.points:
            listx.append(xNew)
            listy.append(yNew)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.plot(listx,listy,color='k')
    
    plt.xlim((min(listx)-400,max(listx)+400))
    plt.ylim((min(listy)-400,max(listy)+400))
    
    palette2 = plt.cm.cool
    palette2.set_under('k',1)
    
    plt.scatter(InputCoord['X'],InputCoord['Y'],c=Val3, cmap=palette2,s=686.66,vmin=0,vmax=+600,marker='s',edgecolor='face') # Default is jet
    plt.colorbar()
    
    plt.savefig('All' + '_AWC.png',dpi = 900)
    plt.show()
    
    ## MAX WATER CONTENT (MWC)
    
    Val4 = InputParamsPlot["MWC"]
    Val4[np.isnan(Val4)] = -999
    
    plt.figure('MWC') #figsize=(50, 50)
    plt.suptitle("MWC", fontsize=18)
    
    listx=[]
    listy=[]
    test = shapefile.Reader(shpFilePath)
    for sr in test.shapeRecords():
        for xNew,yNew in sr.shape.points:
            listx.append(xNew)
            listy.append(yNew)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.plot(listx,listy,color='k')
    
    plt.xlim((min(listx)-400,max(listx)+400))
    plt.ylim((min(listy)-400,max(listy)+400))
    
    palette3 = plt.cm.spring_r
    palette3.set_under('k',1)
    
    plt.scatter(InputCoord['X'],InputCoord['Y'],c=Val4, cmap=palette3,s=686.66,vmin=-0,vmax=+50,marker='s',edgecolor='face') # Default is jet
    plt.colorbar()
    
    plt.savefig('All' + '_MWC.png',dpi = 900)
    plt.show()
    
    ## ERROR
    
    Val5 = SMBM_ErrTab['Err']
    Val5 = Val5[Val5.index[Val5.index.isin(InputParams.columns)]] # Keep only matching indexes
    
    plt.figure('Error') #figsize=(50, 50)
    plt.suptitle("RMSE", fontsize=18)
    
    listx=[]
    listy=[]
    test = shapefile.Reader(shpFilePath)
    for sr in test.shapeRecords():
        for xNew,yNew in sr.shape.points:
            listx.append(xNew)
            listy.append(yNew)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.plot(listx,listy,color='k')
    
    plt.xlim((min(listx)-400,max(listx)+400))
    plt.ylim((min(listy)-400,max(listy)+400))
    
    palette4 = plt.cm.RdYlBu_r
    palette4.set_over('k',1)
    
    plt.scatter(InputCoord['X'],InputCoord['Y'],c=Val5, cmap=palette4,s=686.66,vmin=-0,vmax=+100,marker='s',edgecolor='face') # Default is jet
    plt.colorbar()
    
    plt.savefig('All' + '_Error.png',dpi = 900)
    plt.show()