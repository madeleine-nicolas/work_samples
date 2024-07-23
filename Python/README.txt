PYTHON SCRIPTS

SMBM_Def.py
	
	Defines the Soil Moisture Balance function for a single soil type (SMBM_fun) and for a combination of soil types (SMBMSoilType_fun)

SMBM_Inputs.py
	
	Creates inputs for SMBM for Adil grid
	
SMBM_Opt_All.py
	
	Simulataneously optimizes soil properties for each soil type(AWC and MWC) using InputSoilType percentages and reference recharge data provided by Adil
	(works well for value averaged over all catchment but not for distribution)

SMBM_Opt_CellByCell.py
	
	Calibrates SMBM params (AWC and MWC) at each grid cell; pumping, RET and Cf coefficients are estimated by weighting data
    according to land use surface
	
SMBM_Apply.py
	
	Uses the above defined equations to perform SMBM using InputSoilType percentages, and soil properties (AWC and MWC) for soil types
	Also provides Error table per grid per year using reference recharge data provided by Adil

SMBM_Apply_CellByCell.py

	Uses the set of parameters obtained in SMBM_Opt_CellByCell and inputs from SMBM_inputs
	
SMBM_PlotError.py
	
	To be run after SMBM_apply.py. Plots spatialized error from SMBM model considering reference recharge


#==============================================================================
# HYD MODEL	
#==============================================================================

###### GENERAL ######

HydModel__def.py
	
	Physically based soil moisture model based on Dewandel et al., 2008 used in HydModel_apply

HydModel__sensitivity-test.py

	Estimates sensitivity of HydMod model to all input paramters (Ks, SoilThick, ThetaS, Thresh, Lambda, hbc, Eta)
	
###### TEST RUN ######

HydModel__TestRun.py

	Physically based soil moisture model based on Dewandel et al., 2008, test run with same data as 'RF-2 layers.xls'
	
###### FINE GRID ######

HydModel_FineGrid_inputs.py
	
	Calculates PG, Kc and RET for each fine grid cell for HydMod

HydModel_FineGrid_apply.py
	
	Uses the above defined function to use HydMod for all soil and land use combinations, then distributed across the catchment FOR SINGLE SOIL HORIZON (hE)

HydModel_FineGrid_apply-2lay.py
	
	Uses the above defined function to use HydMod for all soil and land use combinations, then distributed across the catchment FOR  BOTH SOIL HORIZONS (hE and hB)
	
HydModel_FineGrid_resample.py

	Resample HydMod outputs at larger scale using ArcGis, to be compared with Mizan recharge estimations TO BE RUN AFTER HydModel_FineGrid_apply.py or HydModel_FineGrid_apply-2layer.py

HydModel_FineGrid_plot-inputs.py
	
	Plot inputs at a fine grid scale used for HydMod (PG (single map), RET (one per year), Cf (single map)
	
##### SOIL TYPE ######

HydModel_SoilType_apply.py
	
	Uses the above defined function to use HydMod on one soil type (chosen) for a given input climate data.
	
HydModel_SoilTypePer_apply.py
	
	Uses the above defined function to use HydMod for all soil types which is then is weighted according to soil type percentage for each soil grid; SINGLE LAND USE (no pumping, variable Kc)

	
