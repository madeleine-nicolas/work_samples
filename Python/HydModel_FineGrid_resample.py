# -*- coding: utf-8 -*-
"""
Author: Madeleine NICOLAS       
        madeleine.nicolas@univ-rennes1.fr

Purpose: Resample HydMod outputs at larger scale using ArcGis, to be compared with Mizan recharge estimations

Output: Resampled gridded values, plots of said values

Requirements: Must run HydModel_FineGrid_apply or HydModel_FineGrid_apply-2lay before to generate inputs HydModGrid.csv 

Remarks: 

Date:   August 2018
"""

ScenName = '2Lay_2cm'

import os
import arcpy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shapefile

arcpy.CheckOutExtension('Spatial')

InputDir = r"C:\Users\Madeleine\Desktop\MAHESHWARAM\01 - SOIL MOISTURE MODEL\01 - MODELS\02 - HydMod\01 - Outputs\05 - MRO DATA"
OutputDir = r"C:\Users\Madeleine\Desktop\MAHESHWARAM\01 - SOIL MOISTURE MODEL\02 - INPUT DATA"

HydModGrid = pd.read_csv(os.path.join(InputDir, ScenName + "_HydModGrid.csv"))

InputRechYearly = pd.read_csv(os.path.join(OutputDir,"INPUT_ref-rech_2002-2015.txt"),sep='\t',header=0,index_col=0) # Reference recharge
InputCoordRef = pd.read_csv(OutputDir + "\INPUT_coord.txt",sep='\t',header=0,index_col=2)

#==============================================================================
#  SAVE CSV VALUES AS SHAPEFILE
#==============================================================================

try:       
    in_Table = os.path.join(InputDir, ScenName + '_HydModGrid.csv')
    x_coords = "X"
    y_coords = "Y"
    out_Layer = "HydModGrid_points"
            
    spRef = r"Coordinate Systems\Projected Coordinate Systems\Utm\WGS 1984\WGS 1984 UTM Zone 44N.prj"
    arcpy.env.overwriteOutput = True
    arcpy.MakeXYEventLayer_management(in_Table, x_coords, y_coords, out_Layer, spRef)
    arcpy.FeatureClassToFeatureClass_conversion(out_Layer, os.path.join(OutputDir,'GIS'),"HydModGrid_points")
    
except Exception as err:
    
    print(err.args[0])  

#==============================================================================
# RESAMPLE TO LARGE GRID
#==============================================================================

##  SPATIAL JOIN
target_features = os.path.join(OutputDir,"GIS", "DST_grid.shp")
join_features = os.path.join(OutputDir,"GIS", "HydModGrid_points.shp")
out_feature_class = os.path.join(OutputDir,"GIS", "HydModGrid_resampled.shp")

fieldmappings = arcpy.FieldMappings()
fieldmappings.addTable(target_features)
fieldmappings.addTable(join_features)
for findex in HydModGrid.columns[2:]:
    FieldIndex = fieldmappings.findFieldMapIndex(findex)
    fieldmap = fieldmappings.getFieldMap(FieldIndex)
    fieldmap.mergeRule = "mean"
    fieldmappings.replaceFieldMap(FieldIndex, fieldmap)
arcpy.SpatialJoin_analysis(target_features, join_features, out_feature_class,"#", "#",fieldmappings)

#==============================================================================
# SAVE DATA
#==============================================================================

lyr = os.path.join(OutputDir,"GIS", "HydModGrid_resampled.shp") # Get data from attribute table
HydModGrid_resampled = pd.DataFrame(arcpy.da.FeatureClassToNumPyArray(lyr, np.concatenate((['TARGET_FID'],HydModGrid.columns[2:]),axis=0) )) 
HydModGrid_resampled.columns = np.concatenate((['CellNumber'],[str(col)[1:] for col in HydModGrid_resampled.columns][1:]),axis=0) # Rename columns
HydModGrid_resampled = HydModGrid_resampled.loc[HydModGrid_resampled.index.isin(InputRechYearly.index)] # Get only rows with matching indexes to reference recharge
HydModGrid_resampled.to_csv(os.path.join(OutputDir,"INPUT_HydModResamp_1995-2015_" + ScenName + ".txt"),index=False,sep="\t") # Save csv file with gridded values       

#==============================================================================
# PLOT AND SAVE FIGURES
#==============================================================================
    
shpFilePath = OutputDir + "\GIS\watershed.shp"

os.chdir(r"C:\Users\Madeleine\Desktop\MAHESHWARAM\01 - SOIL MOISTURE MODEL\01 - MODELS\02 - HydMod\02 - Figures\05 - MRO DATA")

for yr in HydModGrid_resampled.columns[6:]:
    
    Val = HydModGrid_resampled[yr]
        
    plt.figure('HydModResamp' + str(yr)) #figsize=(50, 50)
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
    
    plt.title("HydMod Resamp recharge")
    plt.xlim((min(listx)-400,max(listx)+400))
    plt.ylim((min(listy)-400,max(listy)+400))
    plt.scatter(InputCoordRef['X'],InputCoordRef['Y'],c=Val, cmap='Blues',s=686.66,vmin=-0,vmax=+1500,marker='s',edgecolor='face') # Default is jet
    plt.colorbar()
    
    plt.savefig(str(yr) + '_' + ScenName + '_HydModResamp.png')
    plt.show()
    
#==============================================================================
# CALCULATE AND PLOT DIFFERENCES TO REFERENCE
#==============================================================================

Diff = InputRechYearly - HydModGrid_resampled

## MEAN DIFFERENCE
Val = Diff.mean(axis=1)
    
plt.figure('DiffRech') #figsize=(50, 50)
plt.suptitle("Recharge diff (mm)", fontsize=18)

listx=[]
listy=[]
test = shapefile.Reader(shpFilePath)
for sr in test.shapeRecords():
    for xNew,yNew in sr.shape.points:
        listx.append(xNew)
        listy.append(yNew)
plt.gca().set_aspect('equal', adjustable='box')
plt.plot(listx,listy,color='k')

plt.title("MeanDiff = " +'%.1f'%(Val.mean()))
plt.xlim((min(listx)-400,max(listx)+400))
plt.ylim((min(listy)-400,max(listy)+400))
plt.scatter(InputCoordRef['X'],InputCoordRef['Y'],c=Val, cmap='RdYlBu',s=686.66,vmin=-150,vmax=+150,marker='s',edgecolor='face') # Default is jet
plt.colorbar()

plt.savefig('Diff' + '_' + ScenName + '_HydModResamp.png', dpi = 900)
plt.show()

## DIFFERENCE PER YEAR
for yr in InputRechYearly.columns:
    
    plt.figure('DiffRech' + str(yr)) #figsize=(50, 50)
    plt.suptitle('Recharge diff ' + str(yr) +  ' (mm)', fontsize=18)
    
    Val = Diff[yr]
    
    listx=[]
    listy=[]
    test = shapefile.Reader(shpFilePath)
    for sr in test.shapeRecords():
        for xNew,yNew in sr.shape.points:
            listx.append(xNew)
            listy.append(yNew)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.plot(listx,listy,color='k')
    
    plt.title("MeanDiff = " +'%.1f'%(Val.mean()))
    plt.xlim((min(listx)-400,max(listx)+400))
    plt.ylim((min(listy)-400,max(listy)+400))
    plt.scatter(InputCoordRef['X'],InputCoordRef['Y'],c=Val, cmap='RdYlBu',s=686.66,vmin=-300,vmax=+300,marker='s',edgecolor='face') # Default is jet
    plt.colorbar()
    
    plt.savefig('Diff' + str(yr) + '_' + ScenName + '_HydModResamp.png', dpi = 900)
    plt.show()
