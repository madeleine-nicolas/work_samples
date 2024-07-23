"""
Author: Madeleine NICOLAS       
        madeleine.nicolas@univ-rennes1.fr

Purpose: MODFLOW 

Output: 

Requirements: MODFLOW-NWT_64.exe, pip in cmd window, Label_lines.py

Remarks: 

Date:   October 2017
"""
##EXECUTE Label_lines.py FIRST FOR GRAPHIC OUTPUT

########## LOAD PACKAGES ##########
import os
import numpy as np
import pandas as pd
import flopy
import flopy.utils.binaryfile as bf  
from math import sqrt
import matplotlib.pyplot as plt

os.chdir("C:\Users\Madeleine\Documents\Doctorat\ARTICLE\Data\Groundwater mound - Numerical") # Sets working directory

########## MODEL GEOMETRY ##########
namemodel = "Basic"
modelnum=1

ztop = 30               # Top of aquifer
zbot = 0                # Bottom of aquifer
nlay = 6                # Number of layers
Ly = 2000               # Length of the model
Lx = Ly                 # Width of the model
delc = 20               # Length of one cell ~function of DEM
delr = 20               # Width of one cell  ~function of DEM
nrow = Ly/delc          # Number of rows
ncol = Lx/delr          # Number of columns
botm = np.linspace(ztop, zbot, nlay + 1) # Layer positions
delv = (ztop - zbot) / nlay              # Thickness of the model layers

initHH = 5              # Initial hydraulic head
rech_val = 1E-6         # Recharge (m/s)

xcoord = np.linspace(delc,Lx,ncol,dtype=np.int32)
ycoord = np.linspace(delr,Ly,nrow,dtype=np.int32)

########## TIME PARAMETERS ##########
timestep=1
nper=1                  # Number of stress periods
nstp = 300              # Number of timesteps/days
perlen = nstp*24*3600   # Period length(s)
steady = False            # Is the simulation steady-state or not
times_num = np.linspace(0,nstp-1,nstp,dtype=np.int32)

########## HYDRAULIC PARAMETERS ##########
K = 1E-4                 # Horizontal hydraulic conductivity
S = 1E-2                  # Storativity 
hk=K # in [m/t]
ss = 1e-4                 # Specific storage
sy = S
laytyp = np.ones((nlay), dtype=np.int32)  # If =0 confined, =1 unconfined, =2 confined/unconfined (T=const), =3 confined/unconfined (T varies)
laytyp[1:] = 3
layvka = 1    # If layvka=0 vka=vertical hydraulic conductivity, if not vka=ratio horizontal to vertical conductivity
vka = 10

# Leakage
vcont=2/((delv/K)+(delv/K))

########## WETTING PARAMETERS ##########
laywet = 1      # Indicates if wetting is active =0 inactive =non-zero active
wetfct = 0.5  # factor used when cell is converted from dry to wet
iwetit = 3   # iteration interval in wetting/drying algorithm
wetdry = -0.1 # combination of the wetting threshold and a flag to indicate which neighboring cells can cause a cell to become wet

## Variables for the BAS package
ibound = np.ones((nlay, nrow, ncol), dtype=np.int32) # Initialize array
ibound[len(ibound)-1][:,len(ibound[1][1])-1]=[-1]*nrow # Boundary conditions
ibound[len(ibound)-1][len(ibound[1])-1,:]=[-1]*ncol # Boundary conditions

strt=np.ones((nlay, nrow, ncol), dtype=np.float32)*initHH

## Variables for the recharge package
rech = np.zeros((nrow, ncol), dtype=np.float32) # Initialize array
rech[0:2,0:2] = rech_val # Recharge

########## FLOPY OBJECTS ##########
modelname = '%s_%d' % (namemodel,modelnum)

mf = flopy.modflow.Modflow(modelname, exe_name='MODFLOW-NWT_64', version='mfnwt')
dis = flopy.modflow.ModflowDis(mf, nlay, nrow, ncol, delr=delr, delc=delc,
                           top=ztop, botm=botm[1:],
                           nper=nper, perlen=perlen, nstp=nstp, steady=steady) # Discretization File (required in all models)
bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt) # Basic package (required in all models)

########## WETTING CAPABILITY PACKAGE ##########
wet = flopy.modflow.ModflowBcf(mf, laycon=laytyp, iwdflg=1, wetfct=wetfct, iwetit=iwetit, hy=hk, vcont=vcont, sf1=ss, sf2=sy, wetdry=wetdry) # Block-Centered Flow Package

########## SOLVER ##########
pcg=flopy.modflow.ModflowPcg(mf,mxiter=100,iter1=100,hclose=0.0001,rclose=0.0001) # Preconditioned Conjugate-Gradient Package

#upw = flopy.modflow.ModflowUpw(mf, hk=hk, vka=vka, sy=sy, ss=ss,laytyp=laytyp,layvka=layvka) # Upstream Weighting Package
#nwt= flopy.modflow.ModflowNwt(mf) # Newton Solver (can only be used with the UPW package)

########## RECHARGE PACKAGE ##########
nrchop = 3 # 1: Recharge to top grid layer only, 2: recharge to layer defined in irch, 3: recharge to highest active cell
rch = flopy.modflow.ModflowRch(mf, nrchop=nrchop, rech=rech)


########## OUTPUT CONTROL ##########
oc = flopy.modflow.ModflowOc(mf) # Output Control Option
mf.write_input() # Write the model input files

########## RUN THE MODEL ##########
success, mfoutput = mf.run_model(silent=False, pause=False) 
if not success:
    raise Exception('MODFLOW did not terminate normally.')
    
########## DATA PROCESSING ##########
#obs_dist = np.array([0,200,400,600])
obs_dist = np.linspace(0,Lx,Lx/100+1) # Observation points every 300 meters from center of basin
obs_x = obs_dist/sqrt(2)
obs_y = obs_dist/sqrt(2)
nobs=len(obs_x)

## Find nearest coordinates in grid
def find_nearest(array,value): # Define a function
    idx = (np.abs(array-value)).argmin()
    return array[idx]

obs_x_near = []
obs_y_near = []
obs_num = np.linspace(0,len(obs_x)-1,len(obs_x),dtype=np.int32)

for i in obs_num:
    obs_x_near.append(find_nearest(xcoord, obs_x[i]))
    obs_y_near.append(find_nearest(ycoord, obs_y[i]))

##################################################################
#plt.figure(1)
#plt.scatter(obs_x_near,obs_y_near)
#for label, x, y in zip(obs_num+1, obs_x_near, obs_y_near):
#    plt.annotate(
#        label,
#        xy=(x, y), xytext=(-20, 20),
#        textcoords='offset points', ha='right', va='bottom',
#        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
#
#plt.show()
##################################################################

## Create the headfile object
headobj = bf.HeadFile(modelname+'.hds')
head_total = headobj.get_alldata(mflay=nlay-1)
#times_tab = np.linspace(0,nstp,nstp/50+1,dtype=np.int32)

for n in obs_num:       
    head_inter=pd.DataFrame(head_total[:,np.ndarray.tolist(ycoord).index(obs_y_near[n]), 
                                       np.ndarray.tolist(xcoord).index(obs_x_near[n])])
    head_inter.columns=[str(int(obs_dist[n]))]
    
    if n == 0:
        head=head_inter
    else:
        head=head.join(head_inter)
        
########## PLOT FIGURES ##########
plt.figure(2)
for i in obs_num:
    plt.plot(head[str(int(obs_dist[i]))])
labelLines(plt.gca().get_lines(),zorder=2.5)

plt.figure(3)
plt.plot(obs_dist,np.transpose(np.asarray(head.iloc[[nstp-1]])))


########### SAVE DATA ##########
#writer = pd.ExcelWriter(modelname +'.xlsx')
#head.to_excel(writer,'Data', index=True)
#
#Params=[[0 for x in range(3)] for y in range(2) ]
#Params[0][:]=[K, S, Lx,ss]
#Params[1][:]=[wetfct,iwetit,wetdry]
#Params=pd.DataFrame(Params)
#
#Params.to_excel(writer,'Params', index=False,header=False)
#writer.save()
