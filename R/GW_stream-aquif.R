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

## Choose borewell from list

b_num=1
borewell_list=c('CH01','CH02','CH03','CH08','CH11','CH15','CH16','CH18')
borewell=borewell_list[b_num]

## Set parameters

D=2

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
  
len_tab=dim(HH_obs)[1]
Output=matrix(0, len_tab) # Initialize table

for (n in (1:len_tab)-1){ #Loop over stream level variations
  for (k in len_tab:(n+1)){ #Loop over all days of observed period
    if (n==0){
      h1=HH_obs$HH_basin[1]
      #h1=356
    }else{
      h1=HH_obs$HH_basin[n]
    }
    h2=HH_obs$HH_basin[n+1]
    Output[k]=Output[k] + (h2-h1)*erfc(distance_list[b_num]/sqrt(4*D*(k-n)*3600*24))
    #Output[k]=Output[k] + (h2-h1)
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

############
### PLOT ###      
############

# Get the range for the x and y axis

obs_range=range(Results["Obs"],na.rm = TRUE)[2]-range(Results["Obs"],na.rm = TRUE)[1]
sim_range=range(Results["Sim"])[2]-range(Results["Sim"])[1]


yrange <- c(min(range(Results["Obs"],na.rm = TRUE)[1],range(Results["Sim"])[1]),
              max(range(Results["Obs"],na.rm = TRUE)[2],range(Results["Sim"])[2]))

xrange <- range(Results$Date)

nlines=2

# Set up the plot

plot(xrange, yrange, type="n", xlab="Date",
     ylab="Hydraulic head (m)" )
colors <- rainbow(nlines)
linetype <- c(1:nlines)
plotchar <- seq(18,18+nlines,1) # Plotting symbol

# Add lines

for (li in 1:nlines) {
  lines(Results[,1], Results[,li+2], type="b", lwd=1.5,
        lty=linetype[li], col=colors[li], pch=plotchar[li])
}