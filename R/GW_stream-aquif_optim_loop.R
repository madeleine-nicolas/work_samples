#########################################################################################################################################
#
# Author: Madeleine NICOLAS       
#         madeleine.nicolas@univ-rennes1.fr
# 
# Purpose: Implementation of analytical solution to model stream/aquifer interaction in fully saturated conditions
#          
# Date: August 2017
#
#########################################################################################################################################

#########################################################################################################################################
#PRELUDE
#########################################################################################################################################

################
### PACKAGES ###      
################

#install.packages("pracma")
require(pracma)

#install.packages("emdbook")
require(emdbook)

#install.packages("xlsx")
require(xlsx) #Make sure Java version is compatible, if necessary install 64 bit

###################
### USER INPUTS ###      
###################

borewell_list=c('CH01','CH02','CH03','CH08','CH11','CH15','CH16','CH18')


####################### Lists of data ########################
distance_list=c(100,67,375,390,400,405,370,340)
niv_const_list=c(346.5,346.5,346.4,346.9,346.5,346.4,344,345)
###############################################################

## Input data names and paths

dir="C:/Users/Madeleine/Documents/Doctorat/ARTICLE/Data_Catherine"
dir_out="C:/Users/Madeleine/Documents/Doctorat/ARTICLE/Data_Catherine/Output"

obs_file_name=paste("INPUTS_",borewell,"_daily.txt",sep="") # Infiltration 
obs_file_dir=paste(dir,"/",obs_file_name,sep="")
HH_obs<- read.csv(obs_file_dir,sep='\t', header=TRUE)
HH_obs$Date=as.Date(HH_obs$Date,format="%d/%m/%Y")

#########################################################################################################################################
#FUNCTION
#########################################################################################################################################

Stream_aquif_optim <- function(x){
  
  D=x

  len_tab=dim(HH_obs)[1]
  Output=matrix(0, len_tab) # Initialize table
  
  for (n in (1:len_tab)-1){ #Loop over stream level variations
    for (k in len_tab:(n+1)){ #Loop over all days of observed period
      if (n==0){
        h1=HH_obs$HH_basin[1]
        h1=356
      }else{
        h1=HH_obs$HH_basin[n]
      }
      h2=HH_obs$HH_basin[n+1]
      Output[k]=Output[k] + (h2-h1)*erfc(distance_list[b_num]/sqrt(4*D*(k-n)*3600*24))
    }
  }
  
  Output=as.data.frame(Output)
  
  Results=as.data.frame(HH_obs$Date)
  colnames(Results)="Date"
  Results["Output"]=Output
  Results["Obs"]=HH_obs[paste("HH_",borewell_list[b_num],sep="")]
  Results["Obs"][Results["Obs"]==9999]=NA
  Results["Sim"]=Output$V1+HH_obs[paste("HH_",borewell_list[b_num],sep="")][1,1]
  
  R_coef=cor(Results$Obs,Results$Sim,use="pairwise.complete.obs")
  RMSE_coef=sqrt(mean((Results$Obs-Results$Sim)^2,na.rm=TRUE))
  
  c(Results,R_coef,RMSE_coef)
}


Drange=c(seq(1,9) %o% 10^(-4:0))[1:38]
Tab_R=matrix(0, length(Drange),length(borewell_list)) 
Tab_RMSE=matrix(0, length(Drange),length(borewell_list)) 

D_opt_R=matrix(0, length(borewell_list)) 
D_opt_RMSE=matrix(0,length(borewell_list)) 
  
for (b_num in 1:length(borewell_list)){
  
  borewell=borewell_list[b_num]
  
  for (i in 1:length(Drange)){
    res=Stream_aquif_optim(Drange[i])
    Tab_R[i,b_num]=res[[5]]
    Tab_RMSE[i,b_num]=res[[6]]
    
  }
  
  D_opt_R[b_num]=Drange[which(Tab_R^2 == max(Tab_R^2), arr.ind = TRUE)][1]
  D_opt_RMSE[b_num]=Drange[which(Tab_RMSE == min(Tab_RMSE), arr.ind = TRUE)][1]
  
}


rownames(Tab_R)=Drange
colnames(Tab_R)=borewell_list
rownames(Tab_RMSE)=Drange
colnames(Tab_RMSE)=borewell_list

