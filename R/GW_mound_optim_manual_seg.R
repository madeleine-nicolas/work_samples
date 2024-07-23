#########################################################################################################################################
#
# Author: Madeleine NICOLAS       
#         madeleine.nicolas@univ-rennes1.fr
# 
# Purpose: Analytical solution rectangular infiltration basin (Hantush 1967), optimisation and segmentation
#          
# Date:   June 2017
#
#########################################################################################################################################

#########################################################################################################################################
#PRELUDE
#########################################################################################################################################

tic()

################
### PACKAGES ###      
################

#install.packages("pracma")
require(pracma)

#install.packages("emdbook")
require(emdbook)

#install.packages("plotly")
require(plotly)

###################
### USER INPUTS ###      
###################

## FILE LOCATION ##

dir="C:/Users/Madeleine/Desktop/Groundwater mound/Data"

inf_file_name="Data_inf.txt" # Infiltration 
obs_file_name="Data_CH01_obs.txt" # Observed hydraulic head

##  PARAMETERS ##

L=115 # Length of basin
W=40 # Width of basin

b_mean=24 # Average aquifer thickness
b=8 # Original aquifer thickness

X=50 # Coordinates from center of rectangle
Y=87

Impulse_duration=180 # Duration of one impulse in days

Fiss_lay_depth= 340.9 

## Optimisation ranges

Srange=c(seq(1,9) %o% 10^(-4:-1))[1:32] # Investigated range of parameters
Krange=c(seq(1,9) %o% 10^(-8:-3))[1:50]

#################
### LOAD DATA ###      
#################

## Observed data

obs_file_dir=paste(dir,"/",obs_file_name,sep="")
h_obs<- read.csv(obs_file_dir,sep='\t', header=TRUE)
h_obs$CH01_WL=h_obs$CH01_HH-Fiss_lay_depth

## Infiltration data

inf_file_dir=paste(dir,"/",inf_file_name,sep="")
Tab_inf <- read.csv(inf_file_dir,sep='\t', header=TRUE)
Tab_inf$Date=as.Date(Tab_inf$Date) # Time data

Tab_inf$Time_days=as.integer(Tab_inf$Date-Tab_inf$Date[1])+1
Tab_inf$Time_sec=Tab_inf$Time_days*24*60*60

Tab_inf$Inf_m_s=Tab_inf$Inf_mm_day/(1000*60*60*24)

len_Tab=dim(Tab_inf)[1]

#########################################################################################################################################
#FUNCTION
#########################################################################################################################################

Hantush_solution <- function(x,y){
  
   # x are input flow parameters
   # y are input time parameters
  
  S=x[1]
  K=x[2]
  
  day=y
  
  alpha= K*b/S
  
  ########################
  ### IMPULSE RESPONSE ###      
  ########################
  
  #Groundwater level response to R=1 (impulse response to be multiplied by actual recharge rate)
  #Impulse duration set as same as input length (we do not study further period)
  
  Impulse=as.data.frame(Tab_inf$Time_sec) 
  colnames(Impulse)=c('Time_sec')
  Impulse$d=sqrt(4*alpha*Impulse$Time_sec)
  
  Impulse$dZdt=b_mean/(2*S)*(erf((L/2+X)/Impulse$d)+erf((L/2-X)/Impulse$d))*
    (erf((W/2+Y)/Impulse$d)+erf((W/2-Y)/Impulse$d))
  
  
  ###################
  ### CONVOLUTION ###      
  ###################
  
  len=length(Impulse$Time_sec)
  time_step=Tab_inf$Time_sec[1]
  
  Output=rep(0, len)
  
  for (i in 2:len){
    
    c=0
    
    i2=min(Impulse_duration,i-1)
    
    for (i_offset in 1:i2){ # Offset timestep
      
      c=c+Tab_inf$Inf_m_s[i-i_offset]*Impulse$dZdt[i_offset]
      
      Output[i]=time_step*c
      
    }
  }
  
  Output=as.data.frame(Output)
  h=sqrt(Output+b^2)
  
  FINAL=as.data.frame(c(h_obs[,c(2,4)],h))
  
  #############
  ### ERROR ###      
  #############
  
  RMSE <- sqrt(mean((FINAL$CH01_WL[1:day]-FINAL$Output[1:day])^2))
  
  RMSE
  
}

#########################################################################################################################################
#########################################################################################################################################



#########################################################################################################################################
#MANUAL OPTIMISATION
#########################################################################################################################################

## DOUBLE OPTIMISATION

seg=40 # Number of days in segmentation

Tab_RMSE=array(0,dim=c(length(Srange),length(Krange),round(len_Tab/seg))) # Initialize tables
Tab_S=array(0,round(len_Tab/seg))
Tab_K=array(0,round(len_Tab/seg))
RMSE_opt=array(0,round(len_Tab/seg))
  
Tab_rownames=sprintf("%11.1e",Srange) # Set row names as S values (truncated)
Tab_colnames=sprintf("%11.1e",Krange) # Set col names as K values (truncated)

rownames(Tab_RMSE)=Tab_rownames
colnames(Tab_RMSE)=Tab_colnames

k=1

for (day in seq(seg,len_Tab,seg)){
  
  cat("k: ", k, "\n")
  
  i=1
  j=1
  
  for (S in Srange){
    
    for (K in Krange){
      
      Tab_RMSE[i,j,k]=Hantush_solution(c(S,K),day)
      j=j+1
      
    }
    
    i=i+1
    j=1
    
  }
  
  Par_opt=which(Tab_RMSE[,,k] == min(Tab_RMSE[,,k]), arr.ind = TRUE)
  
  Tab_S[k]=Srange[Par_opt[1]]
  Tab_K[k]=Krange[Par_opt[2]]
  
  RMSE_opt[k]=Tab_RMSE[Par_opt[1],Par_opt[2],k]
  
  k=k+1
}

RMSE_opt=as.data.frame(RMSE_opt)
Tab_S=as.data.frame(Tab_S)
Tab_K=as.data.frame(Tab_K)

toc()

#########################################################################################################################################
#OUTPUT
#########################################################################################################################################

output_name="Segmented_optimization.csv"

Tab_output=as.data.frame(c(RMSE_opt,Tab_K,Tab_S))

file_out=paste(dir,"/",output_name,sep="") # Chemin en sortie

write.table(Tab_output,file=file_out,quote=FALSE,row.names=FALSE,col.names=FALSE,sep=",")