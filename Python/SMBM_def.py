# -*- coding: utf-8 -*-
"""
Author: Madeleine NICOLAS       
        madeleine.nicolas@univ-rennes1.fr

Purpose: Soil moisture model used in SMBM_Mahesh

Output: Recharge for one specific set of soil properties

Requirements: 

Remarks: 

Date:   January 2018
"""
import pandas as pd
import os 
import numpy as np

os.chdir(r"C:\Users\Madeleine\Desktop\Soil moisture model\03 - PYTHON CODES")

################################
#####    SINGLE FUNCTION   #####
################################

def SMBM_fun(IWC,MWC,AWC):
    
    """Model which calculates Soil Moisture Balance for a single soil type
    Inputs: Initial Water Content (IWC), Maximum surface storage (MWC) and Available Water Content (AWC)
    """
     
    # REFERENCE RECHARGE
    InputRechYearly = pd.read_csv("INPUT_ref-rech_2002-2015.txt",sep='\t',header=0,index_col=0) 
    InputRechMean = InputRechYearly.mean(axis=0)
    InputRechMean.index = pd.to_datetime(InputRechMean.index).year
    
    ## INPUT CLIMATE DATA
    InputClimate = pd.read_csv("INPUT_climate_2002-2015.txt",sep='\t',header=0,index_col=0)
    InputClimate.index=pd.to_datetime(InputClimate.index)

    # Initialize dataframe
    SMBMTab = pd.DataFrame(index=InputClimate.index, columns=[["DeltaR","AW","D/E","RET","Rech"]]) 
   
    ## Stock Variations
    SMBMTab["DeltaR"][InputClimate["Rainfall_mm"]>=MWC] = MWC - InputClimate["PET_mm"]
    SMBMTab["DeltaR"][InputClimate["Rainfall_mm"]<MWC] = InputClimate["Rainfall_mm"] - InputClimate["PET_mm"]
        
    for d in range(0,len(InputClimate)):
        
        ## Initialize data frame
        if d==0:  
            prev = IWC
        else:
            prev = SMBMTab["AW"][d-1]
            
        ## Available Water
        if prev + SMBMTab["DeltaR"][d] > AWC:
            SMBMTab["AW"][d] = AWC
        elif prev + SMBMTab["DeltaR"][d] < 0:
            SMBMTab["AW"][d] = 0
        else:
            SMBMTab["AW"][d] = prev + SMBMTab["DeltaR"][d]
        
        ## Excess or deficit
        SMBMTab["D/E"][d] = SMBMTab["DeltaR"][d] + prev - SMBMTab["AW"][d]
        
        # Real evaporation and Recharge
        if SMBMTab["D/E"][d] > 0:
            SMBMTab["RET"][d] = InputClimate["PET_mm"][d]
            SMBMTab["Rech"][d] = SMBMTab["D/E"][d]
        else:
            SMBMTab["RET"][d] = InputClimate["PET_mm"][d] + SMBMTab["D/E"][d]
            SMBMTab["Rech"][d] = 0
    
    SMBMTabYearly = SMBMTab["Rech"].resample("A").sum()
        
    SMBMTabYearly.index = SMBMTabYearly.index.year
    Err=np.sqrt(((InputRechYearly - SMBMTabYearly) ** 2).mean())
    
    return SMBMTab, Err

####################################
#####    AGGREGATED FUNCTION   #####
####################################

def SMBMSoilType_fun(InputSoilAWC,InputSoilMWC,IWC, InputSoilType, SMBMSoilType, InputRechYearly):

    """Model run which aggregates Soil Moisture Balance for different soil types
    Inputs: array of AWCs and MWCs (must be same length), Initial Water Content (IWC), 
    Soil Type percentages per grid cell (InputSoilType), initialized table with correct colnames (SMBMSoilType),
    Reference recharge at each grid cell per year (InputRechYearly)
    """

    ## Obtain recharge for each soil type
    for c in SMBMSoilType.columns:
        AWC = InputSoilAWC['SoilAWC'][c]
        MWC = InputSoilMWC['SoilMWC'][c]
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
        ErrTab[str(y)] = (SMBMDiscrYearly[int(y)]-InputRechYearly[str(y)])

    ErrTabAv = np.mean(ErrTab,axis=0)
    ErrTabTot = np.sum(ErrTabAv,axis=0)

    print('done')

    return (SMBMDiscrYearly,ErrTab,ErrTabTot)

#################################
#####    CORRECT FUNCTION   #####
#################################


def SMBM_Cell_fun(gr, AWC, MWC, IWC, InputTotal, InputRET, InputCf, InputRechYearly):
    
    """Model which calculates Soil Moisture Balance for a single grid cell considering pumping and return flow
    Inputs: Obtained from SMBM_Inputs
    Outputs: SMBM for all period for grid cell, yearly values and error compared to obs
    """

    InputTotalCell = pd.DataFrame(InputTotal[gr]) # Rainfall and pumping
    InputTotalCell.columns = ['Inputs_mm']
    
    InputRETCell = pd.DataFrame(InputRET[gr]) # Real evapotranspiration
    InputRETCell.columns = ['RET_mm']
    
    # Initialize dataframe
    SMBMTabCell = pd.DataFrame(index=InputTotalCell.index, columns=[["DeltaR","AW","D/E","RET","TotRech","NatRech","Runoff"]]) 
   
    ## Stock Variations
    SMBMTabCell["DeltaR"][InputTotalCell["Inputs_mm"]>=MWC] = MWC - InputRETCell["RET_mm"]
    SMBMTabCell["DeltaR"][InputTotalCell["Inputs_mm"]<MWC] = InputTotalCell["Inputs_mm"] - InputRETCell["RET_mm"]
    
    SMBMTabCell["Runoff"][InputTotalCell["Inputs_mm"]<=MWC] = 0
    SMBMTabCell["Runoff"][InputTotalCell["Inputs_mm"]>MWC] = InputTotalCell["Inputs_mm"] - MWC
        
    for d in range(0,len(InputTotalCell)):
        
        ## Initialize data frame
        if d==0:  
            prev = IWC
        else:
            prev = SMBMTabCell["AW"][d-1]
            
        ## Available Water
        if prev + SMBMTabCell["DeltaR"][d] > AWC:
            SMBMTabCell["AW"][d] = AWC
        elif prev + SMBMTabCell["DeltaR"][d] < 0:
            SMBMTabCell["AW"][d] = 0
        else:
            SMBMTabCell["AW"][d] = prev + SMBMTabCell["DeltaR"][d]
        
        ## Excess or deficit
        SMBMTabCell["D/E"][d] = SMBMTabCell["DeltaR"][d] + prev - SMBMTabCell["AW"][d]
        
        # Real evaporation and Recharge
        if SMBMTabCell["D/E"][d] > 0:
            SMBMTabCell["RET"][d] = InputRETCell["RET_mm"][d]
            SMBMTabCell["TotRech"][d] = SMBMTabCell["D/E"][d]
        else:
            SMBMTabCell["RET"][d] = InputRETCell["RET_mm"][d] + SMBMTabCell["D/E"][d]
            SMBMTabCell["TotRech"][d] = 0
    
    ## Partition natural recharge from recharge flow with RF Coef (Cf)
    SMBMTabCell["NatRech"] = pd.DataFrame(InputCf[gr]).multiply(SMBMTabCell["TotRech"],axis="index")
    
    ## Yearly values
    #Sim
    SMBMRechCellYearly = pd.DataFrame(SMBMTabCell["NatRech"].resample("A").sum())
    SMBMRechCellYearly.index = SMBMRechCellYearly.index.year
    SMBMRechCellYearly.columns = ['RechCalc']
    #Obs
    RefRechCellYearly = pd.DataFrame(InputRechYearly[gr])
    RefRechCellYearly.columns = ['RechObs']
    
    RechCellYearly = RefRechCellYearly.join(SMBMRechCellYearly)
    
    # Error calculation   
    Err = np.sqrt(((RechCellYearly['RechCalc'][0:] - RechCellYearly['RechObs'][0:]) ** 2).mean())

    return SMBMTabCell, Err, RechCellYearly
    