#########################################################################################################################################
#
# Author: Madeleine NICOLAS       
#         madeleine.nicolas@univ-rennes1.fr
# 
# Purpose: Analytical solution rectangular infiltration basin (Hantush 1967) and optimisation
#          
# Date:   June 2017
#
#########################################################################################################################################

#########################################################################################################################################
#PRELUDE
#########################################################################################################################################

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

tic()

###################
### USER INPUTS ###      
###################

## Choose borewell from list

b_num=1 
borewell_list=c('CH01','CH02','CH03','CH08','CH11','CH15','CH16','CH18')
borewell=borewell_list[b_num]

## Set parameters

L=115 # Length of basin
W=40 # Width of basin

Impulse_duration=180 # Duration of one impulse in days

####################### Lists of data ###########################

X_list=c(50,3,215,247,219,177,167,211)
Y_list=c(87,66,305,300,334,293,328,346)
b_mean_list=c(8.4,8.8,1.0,1.6,1.1,3.6,10.0,4.1)
b_list=c(7.4,7.9,0.0,0.6,0.1,2.4,7.5,1)
fisslay_list=c(341.7,341.3,347.2,346.6,347.1,345.0,337.0,343.6)

X=X_list[b_num] # X coordinate from center of rectangle
Y=Y_list[b_num] # Y coordinate from center of rectangle
b_mean=b_mean_list[b_num] # Average aquifer thickness
b=b_list[b_num] # Original aquifer thickness
Fiss_lay_depth=fisslay_list[b_num] # Fissured layer depth

##################################################################

## Input data names and paths

dir="C:/Users/Madeleine/Documents/Doctorat/ARTICLE/Data_Groundwater mound/Data"
dir_out="C:/Users/Madeleine/Documents/Doctorat/ARTICLE/Data_Groundwater mound/Output"

# Hydraulic head
obs_file_name=paste("Data_",borewell,"_obs.txt",sep="")
obs_file_dir=paste(dir,"/",obs_file_name,sep="") 
HH_obs<- read.csv(obs_file_dir,sep='\t', header=TRUE)
HH_obs[paste(borewell,"_WL",sep="")]=HH_obs[paste(borewell,"_HH",sep="")]-Fiss_lay_depth

# Infiltration data
inf_file_name="Data_inf.txt"
inf_file_dir=paste(dir,"/",inf_file_name,sep="")
Tab_inf <- read.csv(inf_file_dir,sep='\t', header=TRUE)
Tab_inf$Date=as.Date(Tab_inf$Date) # Time data
Tab_inf$Time_days=as.integer(Tab_inf$Date-Tab_inf$Date[1])
Tab_inf$Time_sec=Tab_inf$Time_days*24*60*60
Tab_inf$Inf_m_s=Tab_inf$Inf_mm_day/(1000*60*60*24)

#########################################################################################################################################
#FUNCTION
#########################################################################################################################################

Hantush_solution <- function(x){
  
  S=x[1]
  K=x[2]
  
  alpha= K*b_mean/S
  
  ########################
  ### IMPULSE RESPONSE ###      
  ########################
  
  #Groundwater level response to R=1 (impulse response to be multiplied by actual recharge rate)
  
  Impulse=as.data.frame(Tab_inf$Time_sec)
  colnames(Impulse)=c('Time_sec')
  Impulse$d=sqrt(4*alpha*Impulse$Time_sec)
  
  Impulse$dZdt=b_mean/(2*S)*(erf((L/2+X)/Impulse$d)+erf((L/2-X)/Impulse$d))*
    (erf((W/2+Y)/Impulse$d)+erf((W/2-Y)/Impulse$d))
  
  
  ###################
  ### CONVOLUTION ###      
  ###################
  
  len=length(Impulse$Time_sec)
  time_step=Tab_inf$Time_sec[2]
  
  Output=rep(0, len)
  
  for (i in 2:len){
    
    c=0
    
    i2=min(Impulse_duration,i-1)
    
    for (i_offset in 1:i2){ # Offset timestep
      
      c=c+Tab_inf$Inf_m_s[i-i_offset]*Impulse$dZdt[i_offset]
      
    }
    
    Output[i]=time_step*c
    
  }
  
  Output=as.data.frame(Output)
  h=sqrt(Output+b^2)
  h[paste(borewell,"_Days",sep="")]=Tab_inf$Time_days
  
  Results=merge(h,HH_obs[c(paste(borewell,"_Days",sep=""),paste(borewell,"_WL",sep=""))], by=paste(borewell,"_Days",sep=""),all.x = TRUE)
  
  Results=Results[,c(1,3,2)]

  #############
  ### ERROR ###      
  #############
  
  RMSE <- sqrt(mean((Results[paste(borewell,"_WL",sep="")]-Results$Output)^2,na.rm=TRUE))
  
  RMSE_1 <- sqrt(mean((Results[1:80,paste(borewell,"_WL",sep="")]-Results$Output[1:80])^2,na.rm=TRUE))
  
  c(RMSE,RMSE_1)
  
  }

#########################################################################################################################################
#########################################################################################################################################



#########################################################################################################################################
#OPTIMISATION
#########################################################################################################################################

Srange=c(seq(1,9) %o% 10^(-4:-2))[1:20]
Krange=c(seq(1,9) %o% 10^(-6:-4))

# Srange=c(seq(1,9) %o% 10^(-4:-2))[5:27] # Investigated range of parameters
# Krange=c(seq(1,9) %o% 10^(-8:-3))[5:46]

# Srange=c(seq(1,9) %o% 10^(-3:-2))[5:14] # Investigated range of parameters
# Krange=c(seq(1,9) %o% 10^(-6:-5))[5:14]

Tab_RMSE=matrix(0, length(Srange),length(Krange)) # Initialize tables
Tab_RMSE_1=matrix(0, length(Srange),length(Krange))

for (n in 1:length(Srange)){
  
  for (m in 1:length(Krange)){
    
    res=Hantush_solution(c(Srange[n],Krange[m]))
    
    Tab_RMSE[n,m]=res[1]
    Tab_RMSE_1[n,m]=res[2]

  }
  
}

Tab_rownames=sprintf("%11.1e",Srange) # Set row names as S values (truncated)
Tab_colnames=sprintf("%11.1e",Krange) # Set col names as K values (truncated)

rownames(Tab_RMSE)=Tab_rownames
colnames(Tab_RMSE)=Tab_colnames

rownames(Tab_RMSE_1)=Tab_rownames
colnames(Tab_RMSE_1)=Tab_colnames

Par_opt_all=which(Tab_RMSE == min(Tab_RMSE), arr.ind = TRUE)
Par_opt_1=which(Tab_RMSE_1 == min(Tab_RMSE_1), arr.ind = TRUE)

S_opt_all=Srange[Par_opt_all[1]]
K_opt_all=Krange[Par_opt_all[2]]

S_opt_1=Srange[Par_opt_1[1]]
K_opt_1=Krange[Par_opt_1[2]]

RMSE_all=Tab_RMSE[Par_opt_all]
RMSE_1=Tab_RMSE_1[Par_opt_1]

#########################################################################################################################################
#SAVE
#########################################################################################################################################

# setwd(dir=dir_out)
# write.xlsx(Tab_RMSE,file=paste(dir_out,"/",borewell_list[b_num],"_Tab_RMSE.xlsx",sep=""),sheetName="RMSE", row.names=TRUE)
# write.xlsx(Tab_RMSE_1,file=paste(dir_out,"/",borewell_list[b_num],"_Tab_RMSE.xlsx",sep=""),sheetName="RMSE1", row.names=TRUE,append=TRUE)
# 
# Par_tab=rbind(c(K_opt_1,K_opt_2,K_opt_all),c(S_opt_1,S_opt_2,S_opt_all))
# write.xlsx(Par_tab,file=paste(dir_out,"/",borewell_list[b_num],"_Par_tab.xlsx",sep=""),sheetName="Optimal params", row.names=FALSE)

toc()
beep(2)

#########################################################################################################################################
#3D PLOT
#########################################################################################################################################

# Tab_RMSE_plot=Tab_RMSE
# Tab_RMSE_plot[Tab_RMSE_plot>12]=NA
# 
# Tab_RMSE_plot_1=Tab_RMSE_1
# Tab_RMSE_plot_1[Tab_RMSE_plot_1>12]=NA

#plot_ly(x = Krange, y = Srange, z = Tab_RMSE_plot) %>% add_surface()
#plot_ly(x = Krange, y = Srange, z = Tab_RMSE_plot_1) %>% add_surface()
#plot_ly(x = Krange, y = Srange, z = Tab_RMSE_plot_2) %>% add_surface()