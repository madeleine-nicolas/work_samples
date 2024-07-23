#########################################################################################################################################
#
# Author: Madeleine NICOLAS       
#         madeleine.nicolas@univ-rennes1.fr
# 
# Purpose: Calculates profile of groundwater from infiltration using Hantush analytical solution 
#          
# Date:   June 2017
#
#########################################################################################################################################

#install.packages("pracma")
require(pracma)

#install.packages("xlsx")
library(xlsx)

#########################################################################################################################################
#USER INPUTS
#########################################################################################################################################

## FILE LOCATION ##

dir_in="C:/Users/Madeleine/Documents/Doctorat/ARTICLE/Data_R/Inputs"
dir_out="C:/Users/Madeleine/Documents/Doctorat/ARTICLE/Data_R/Outputs"

inf_file_name="INPUTS_Infiltration.txt" # Infiltration 
profile_file_name="INPUTS_Coordinates.txt" # Profile coordinates

##  PARAMETERS ##

L=115 # Length of basin
W=40 # Width of basin

b_mean=4.8 # Average aquifer thickness
b=3.3 # Original aquifer thickness

Impulse_duration=180 # Duration of one impulse in days

S = 1.9E-2
K = 3.6E-4

#########################################################################################################################################
#LOAD DATA
#########################################################################################################################################

## Profile coordinate data

profile_file_dir=paste(dir_in,"/",profile_file_name,sep="")
profile_coord<- read.csv(profile_file_dir,sep='\t', header=TRUE)

## Infiltration data

inf_file_name="INPUTS_Infiltration.txt"
inf_file_dir=paste(dir_in,"/",inf_file_name,sep="")
Tab_inf <- read.csv(inf_file_dir,sep='\t', header=TRUE)
Tab_inf$Date=as.Date(Tab_inf$Date,format="%d/%m/%Y") # Time data
Tab_inf$N_days=as.integer(Tab_inf$Date-Tab_inf$Date[1])
Tab_inf$Time_sec=Tab_inf$N_days*24*60*60
Tab_inf$Inf_m_s=Tab_inf$Inf_mm_day/(1000*60*60*24)

len_Tab=dim(Tab_inf)[1]

#########################################################################################################################################
#FUNCTION
#########################################################################################################################################

Hantush_solution_profile <- function(x){
  
  X=x[1]
  Y=x[2]
  
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
  
  h

}

#########################################################################################################################################
#APPLICATION
#########################################################################################################################################

Mound_profile=as.data.frame(matrix(0, len_Tab,dim(profile_coord)[1]))

for (point in 1:dim(profile_coord)[1]){
  h <-Hantush_solution_profile(c(profile_coord[point,'X'],profile_coord[point,'Y']))
  Mound_profile[,point]=h
}

Mound_profile=as.data.frame(Mound_profile)
colnames(Mound_profile)=profile_coord$Name

output_name <- readline(prompt="Enter output file name: ")

write.xlsx(Mound_profile, paste(dir_out,"/",output_name,".xlsx",sep=""))

#########################################################################################################################################
#PLOT
#########################################################################################################################################

setwd(dir=dir_out)

pdf(paste(output_name,".pdf",sep=""))

for (i in 1:244){
  
  yrange <- range(Mound_profile)

  xrange <- range(profile_coord$Distance)

  # set up the plot

  plot(xrange, yrange, type="n", xlab="Distance (m)",
       ylab="Groundwater mound (m)" )
  colors <- "#000000"
  linetype <- 1
  plotchar <- 3 # Plotting symbol
  
  lines(profile_coord$Distance, Mound_profile[i,], type="b", lwd=1.5,
        lty=linetype, col=colors, pch=plotchar)
}

dev.off()
