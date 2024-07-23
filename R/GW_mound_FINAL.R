#########################################################################################################################################
#
# Author: Madeleine NICOLAS       
#         madeleine.nicolas@univ-rennes1.fr
# 
# Purpose: Implements functions from GW_Hantush_functions to perform optimisation of parameters for all boreholes
#          
# Date:   August 2017
#
#########################################################################################################################################

#########################################################################################################################################
#PRELUDE
#########################################################################################################################################

tic()

################
### PACKAGES ###      
################

#install.packages("beepr")
library(beepr)

#install.packages("pracma")
require(pracma)


#install.packages("emdbook")
require(emdbook)

#install.packages("plotly")
require(plotly)

#install.packages("xlsx")
require(xlsx) #Make sure Java version is compatible, if necessary install 64 bit

###################
### USER INPUTS ###      
###################

## Choose borewell from list

b_num=3 
borewell_list=c('CH01','CH02','CH03','CH08','CH11','CH15','CH16','CH18')
borewell=borewell_list[b_num]

## Choose period to fit over 

nfit=1
fit=c("1","2","all")

## Set parameters

L=115 # Length of basin
W=40 # Width of basin

Impulse_duration=180 # Duration of one impulse in days

####################### Lists of data ########################
X_list=c(50,3,215,247,219,177,167,211)
Y_list=c(87,66,305,300,334,293,328,346)
b_mean_list=c(8.4,8.8,1.0,1.6,1.1,3.6,10.0,4.1)
#b_total_list=c(23.1,24.1,17.6,18.5,17.9,19.0,27.2,13.5)
b_list=c(7.4,7.9,0.0,0.6,0.1,2.4,7.5,1)
fisslay_list=c(341.7,341.3,347.2,346.6,347.1,345.0,337.0,343.6)

X=X_list[b_num] # X coordinate from center of rectangle
Y=Y_list[b_num] # Y coordinate from center of rectangle
b_mean=b_mean_list[b_num] # Average aquifer thickness
b=b_list[b_num] # Original aquifer thickness
Fiss_lay_depth=fisslay_list[b_num] # Fissured layer depth
###############################################################


## Input data names and paths

dir="C:/Users/Madeleine/Documents/Doctorat/ARTICLE/Data_Groundwater mound/Data"
dir_out="C:/Users/Madeleine/Documents/Doctorat/ARTICLE/Data_Groundwater mound/Output"

# Infiltration data
inf_file_name="Data_inf.txt"
inf_file_dir=paste(dir,"/",inf_file_name,sep="")
Tab_inf <- read.csv(inf_file_dir,sep='\t', header=TRUE)
Tab_inf$Date=as.Date(Tab_inf$Date) # Time data
Tab_inf$Time_days=as.integer(Tab_inf$Date-Tab_inf$Date[1])
Tab_inf$Time_sec=Tab_inf$Time_days*24*60*60
Tab_inf$Inf_m_s=Tab_inf$Inf_mm_day/(1000*60*60*24)

#########################################################################################################################################
#LOOP
#########################################################################################################################################

## Hydrodynamic properties range of investigation
Srange=c(seq(1,9) %o% 10^(-4:-2))[5:27] # Investigated range of parameters
Krange=c(seq(1,9) %o% 10^(-8:-3))[5:46]

for (b_num in 1:length(borewell_list)){
  
  ## PARAMETERS ##
  
  ## Borewell
  borewell=borewell_list[b_num]
  X=X_list[b_num] # Coordinates from center of rectangle
  Y=Y_list[b_num]
  b_mean=b_mean_list[b_num] # Average aquifer thickness
  b=b_list[b_num] # Original aquifer thickness
  Fiss_lay_depth=fisslay_list[b_num]
  
  ## Load observed data
  obs_file_name=paste("Data_",borewell,"_obs.txt",sep="") # Observed hydraulic head
  obs_file_dir=paste(dir,"/",obs_file_name,sep="")
  h_obs<- read.csv(obs_file_dir,sep='\t', header=TRUE)
  h_obs[paste(borewell,"_WL",sep="")]=h_obs[paste(borewell,"_HH",sep="")]-Fiss_lay_depth

  
  ## OPTIMISATION ##
  
  ## Initialize tables
  Tab_RMSE=matrix(0, length(Srange),length(Krange)) 
  Tab_RMSE_1=matrix(0, length(Srange),length(Krange))
  Tab_RMSE_2=matrix(0, length(Srange),length(Krange))
  
  for (n in 1:length(Srange)){
    
    for (m in 1:length(Krange)){
      
      res=Hantush_optimisation(c(Srange[n],Krange[m]))
      
      Tab_RMSE[n,m]=res[1]
      Tab_RMSE_1[n,m]=res[2]
      Tab_RMSE_2[n,m]=res[3]
      
    }
    
  }
  
  ## Row and column names
  Tab_rownames=sprintf("%11.1e",Srange) # Set row names as S values (truncated)
  Tab_colnames=sprintf("%11.1e",Krange) # Set col names as K values (truncated)
  rownames(Tab_RMSE)=Tab_rownames
  colnames(Tab_RMSE)=Tab_colnames
  rownames(Tab_RMSE_1)=Tab_rownames
  colnames(Tab_RMSE_1)=Tab_colnames
  rownames(Tab_RMSE_2)=Tab_rownames
  colnames(Tab_RMSE_2)=Tab_colnames
  
  ## Optimum parameters
  Par_opt_all=which(Tab_RMSE == min(Tab_RMSE), arr.ind = TRUE)
  Par_opt_1=which(Tab_RMSE_1 == min(Tab_RMSE_1), arr.ind = TRUE)
  Par_opt_2=which(Tab_RMSE_2 == min(Tab_RMSE_2), arr.ind = TRUE)
  S_opt_all=Srange[Par_opt_all[1]]
  K_opt_all=Krange[Par_opt_all[2]]
  S_opt_1=Srange[Par_opt_1[1]]
  K_opt_1=Krange[Par_opt_1[2]]
  S_opt_2=Srange[Par_opt_2[1]]
  K_opt_2=Krange[Par_opt_2[2]]
  
  RMSE_all=Tab_RMSE[Par_opt_all]
  RMSE_1=Tab_RMSE_1[Par_opt_1]
  RMSE_2=Tab_RMSE_2[Par_opt_2]
  
  #########################################################################################################################################
  #SAVE
  #########################################################################################################################################
  
  setwd(dir=dir_out)
  write.xlsx(Tab_RMSE,file=paste(dir_out,"/",borewell_list[b_num],"_Tab_RMSE.xlsx",sep=""),sheetName="RMSE", row.names=TRUE)
  write.xlsx(Tab_RMSE_1,file=paste(dir_out,"/",borewell_list[b_num],"_Tab_RMSE.xlsx",sep=""),sheetName="RMSE1", row.names=TRUE,append=TRUE)
  write.xlsx(Tab_RMSE_2,file=paste(dir_out,"/",borewell_list[b_num],"_Tab_RMSE.xlsx",sep=""),sheetName="RMSE2", row.names=TRUE,append=TRUE)
  
  Par_tab=rbind(c(K_opt_1,K_opt_2,K_opt_all),c(S_opt_1,S_opt_2,S_opt_all))
  write.xlsx(Par_tab,file=paste(dir_out,"/",borewell_list[b_num],"_Par_tab.xlsx",sep=""),sheetName="Optimal params", row.names=FALSE)
  
  for (nfit in 1:length(fit)){
    
    Kfit=eval(as.name(paste("K_opt_",fit[nfit],sep="")))
    Sfit=eval(as.name(paste("S_opt_",fit[nfit],sep="")))
    
    Hantush_solution(c(Sfit,Kfit,nfit))
  }
  
  beepr::beep(2)
  
}

beepr::beep(8)
toc()
  