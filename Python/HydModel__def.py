# -*- coding: utf-8 -*-
"""
Author: Madeleine NICOLAS       
        madeleine.nicolas@univ-rennes1.fr

Purpose: Physically based soil moisture model based on Dewandel et al., 2008 used in HydModel_apply

Requirements: Input ETR, R and PG data; hydraulic parameters 

Remarks: 

Date:   April 2018
"""

import numpy as np
import pandas as pd

########## HYDRAULIC MODEL ##########

def HydMod_fun(Ks,SoilThick,ThetaS,Thresh,Lambda,hbc,Eta,InputData):
    
    """Physically based hydraulic model to calculate recharge
    
    Inputs: Dataframe containing Potential Evapotranspiration (PET), Rainfall (R), 
    Pumping (PG), and Real Evapotranspiration (RET) at a daily time-step
    
    Outputs: Dataframe containing different components of hydraulic model (namely qout) at a daily time-step
    """
  
    #OutputData = pd.DataFrame()
    OutputData = pd.DataFrame(0.0,index=InputData.index, columns=['Rd','Runoff','Theta1','hUnsat','Kh','qh','Theta2','hs','qs','qout','DeltaS','Deficit','AET']).reset_index()
    
    for i in range(1,len(OutputData)):
        
        # Flux entering field at beggining of timestep (Rd)
        if 'PG' in InputData.columns:
            Rdbis = InputData['PG'][i] + InputData['R'][i] + OutputData['DeltaS'][i-1] * (OutputData['DeltaS'][i-1]>=0)
        else: 
            Rdbis = InputData['R'][i] + OutputData['DeltaS'][i-1] * (OutputData['DeltaS'][i-1]>=0)
            
        if Rdbis <= Thresh:
            OutputData.at[i,'Rd'] = Rdbis
            OutputData.at[i,'Runoff'] = 0
        else:
            OutputData.at[i,'Rd'] = Thresh
            OutputData.at[i,'Runoff'] = Rdbis-Thresh
            
        # Soil moisture at beggining of timestep (Theta1)
        Theta1bis = OutputData['Theta2'][i-1]+(OutputData['Rd'][i]-InputData['RET'][i])/SoilThick
        if  Theta1bis >= ThetaS:
            OutputData.at[i,'Theta1'] = ThetaS
        elif Theta1bis > 0:        
            OutputData.at[i,'Theta1'] = Theta1bis
        else:
            OutputData.at[i,'Theta1'] = 1E-10 # To allow calculation of hUnsat..?
        
        # Hydraulic pressure in unsaturated zone (hUnsat) (Brooks and Corey)
        if OutputData['Theta1'][i] > 0:
            OutputData.at[i,'hUnsat'] = hbc*np.power((ThetaS/OutputData['Theta1'][i]),(1/Lambda)) #in m
            
        # Unsaturated hydraulic conductivity (Kh) (Brooks and Corey)
        if OutputData['Theta1'][i] < ThetaS:
           OutputData.at[i,'Kh'] = Ks*np.power((OutputData['Theta1'][i]/ThetaS),Eta) #in m/day
        
        # Flux from unsaturated zone (qh)
        qmax = OutputData['Rd'][i] - InputData['RET'][i] + OutputData['Theta2'][i-1]*SoilThick
        qhbis = -OutputData['Kh'][i]*(OutputData['hUnsat'][i]/SoilThick-1)
        if (OutputData['Kh'][i] > 0):
            if (qhbis < qmax):
                OutputData.at[i,'qh'] = qhbis
            elif qmax > 0:
                OutputData.at[i,'qh'] = qmax
        
        # Saturated thickness (hs)
        if OutputData['Theta1'][i] == ThetaS:
           OutputData.at[i,'hs'] = OutputData['Rd'][i] - InputData['RET'][i] - (ThetaS - OutputData['Theta2'][i-1]) * SoilThick + SoilThick
        
        # Flux from saturated zone (qs)
        qsbis = (OutputData['hs'][i])/SoilThick*Ks
        if (OutputData['qh'][i] == 0) & (OutputData['Theta1'][i] == ThetaS) & (OutputData['hs'][i] >= 0):
            if qsbis <= qmax:
                OutputData.at[i,'qs'] = qsbis
            elif qsbis > qmax:
                OutputData.at[i,'qs'] = qmax
                
        # Total outflow  (qout)
        OutputData.at[i,'qout'] = OutputData['qs'][i] + OutputData['qh'][i]
        
        # Saturated media budget (DeltaS)
        if (OutputData['qh'][i] == 0) & (OutputData['qs'][i] >= 0):
            OutputData.at[i,'DeltaS'] = OutputData['Rd'][i] - InputData['RET'][i] - (ThetaS - OutputData['Theta2'][i-1])*SoilThick - OutputData['qout'][i]
            
        # Soil moisture at end of timestep (Theta2)
        Theta2_1 = 0
        if (OutputData['qh'][i] >= 0) & (OutputData['Theta1'][i] < ThetaS): # If media is unsaturated
            Theta2_1 = OutputData['Theta1'][i] - OutputData['qh'][i]/SoilThick 
        
        Theta2_2bis = OutputData['Theta1'][i] + (OutputData['Rd'][i] - InputData['RET'][i] - (ThetaS - OutputData['Theta2'][i-1])*SoilThick - OutputData['qout'][i])/SoilThick
        
        if (Theta2_1 == 0)  & (OutputData['qs'][i] >= 0) & (OutputData['Theta1'][i] == ThetaS) & (Theta2_2bis > 0) & (OutputData['DeltaS'][i] < 0):  # If media is saturated
            Theta2_2 = Theta2_2bis
        else:
            Theta2_2 = 1E-10
            
        if (Theta2_1 == 0)  & (OutputData['qs'][i] >= 0) & (OutputData['Theta1'][i] == ThetaS) & (OutputData['DeltaS'][i] >= 0):
            Theta2_2 += ThetaS - 1E-10            
        
        if (Theta2_1 >= 0) & (Theta2_2 >= 0):
            OutputData.at[i,'Theta2'] = Theta2_1 + Theta2_2
        else:
            OutputData.at[i,'Theta2'] = 0
            
        # Soil moisture deficit
        if (OutputData['Theta1'][i] == 1E-10) & (Theta1bis < 0):
            OutputData.at[i,'Deficit'] = -(OutputData['Theta2'][i-1]*SoilThick + OutputData['Rd'][i]- InputData['RET'][i])
        
        # Actual evapotranspiration
        OutputData.at[i,'AET'] = InputData['RET'][i] - OutputData['Deficit'][i]
         
    return(OutputData)   
    
### APPLICATION EXAMPLE
#
### Inputs
#import HydModel_FineGrid_inputs_THESE as HydModel_inputs
#
#InputClimate = HydModel_inputs.InputClimate
#InputPG = HydModel_inputs.InputPG 
#InputPG_Yrly = InputPG.resample("A").sum()
#InputPG_Yrly.index = InputPG_Yrly.index.year
#InputRET = HydModel_inputs.InputRET
#InputCf = HydModel_inputs.InputCf 
#    
### hE
#Thresh = 0.03
#Ks = np.multiply(1.21E-05,86400) # Saturated Hydraulic conductivity (m/s)
#ThetaS = 0.35 # Saturated soil moisture content (m3/m3)
#Lambda = 0.129 # Pore-size-distribution parameter
#hbc = -0.06 # Air-entry pressure head (m)
#Eta = 19.3 # Texture-dependent conductivity shape parameter
#SoilThick = 0.4 # Soil thickness(m)
#
#HydMod_fun(Ks,SoilThick,ThetaS,Thresh,Lambda,hbc,Eta,InputData)
#
### hB
#Ks2 = np.multiply(1.08E-5,86400)
#ThetaS2 =  0.35
#Lambda2 = 0.096
#hbc2 = -0.123
#Eta2 = 24.1
#SoilThick2 = 1.13 # Soil thickness (hB) (m)

