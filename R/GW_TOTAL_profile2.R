#########################################################################################################################################
#
# Author: Madeleine NICOLAS       
#         madeleine.nicolas@univ-rennes1.fr
# 
# Purpose: 
#          
# Date:   
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

lim_P1=60 # Limit phase 1
lim_P2=110 # Limit phase 2

b_mean=4.8 # Average aquifer thickness
b=3.3 # Original aquifer thickness

Impulse_duration=180 # Duration of one impulse in days

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

#########################################################################################################################################
########################################################   P1   #########################################################################
#########################################################################################################################################

Hantush_solution_TOTAL <- function(x){
  
  X=x[1]
  Y=x[2]
  
  len=length(Tab_inf$N_days)
  Output_P1=matrix(0, len) # Initialize table
  
  K_P1=3.6e-4
  S_P1=1.9e-2
  
  ### IMPULSE RESPONSE ###      
  
  alpha= K_P1*b_mean/S_P1
  Impulse_P1=as.data.frame(Tab_inf$Time_sec)
  colnames(Impulse_P1)=c('Time_sec')
  Impulse_P1$d=sqrt(4*alpha*Impulse_P1$Time_sec)
  
  Impulse_P1$dZdt=b_mean/(2*S_P1)*(erf((L/2+X)/Impulse_P1$d)+erf((L/2-X)/Impulse_P1$d))*
    (erf((W/2+Y)/Impulse_P1$d)+erf((W/2-Y)/Impulse_P1$d))
  
  ### CONVOLUTION ###      
  
  time_step=Tab_inf$Time_sec[2]
  
  for (day in 2:len){
    
    c=0
    
    i2=min(Impulse_duration,day-1)
    
    for (i_offset in 1:i2){ # Offset timestep
      
      if (i_offset<lim_P1){
        
        Inf_offset=Tab_inf$Inf_m_s[day-i_offset]
        
      }else{
        
        Inf_offset=0
        
      }
      
      c=c+Inf_offset*Impulse_P1$dZdt[i_offset]
      
    }
    
    Output_P1[day]=time_step*c
    
  }
  
  Output_P1=as.data.frame(Output_P1)
  h_P1=sqrt(Output_P1+b^2)
  h_P1["N_days"]=Tab_inf$N_days
  colnames(h_P1)=c("h","N_days")
  
  #########################################################################################################################################
  ########################################################   P2   #########################################################################
  #########################################################################################################################################
  
  ############################## Function ###############################
  
  Output_P2=matrix(0, len) # Initialize table
  
  ### IMPULSE RESPONSE ###   
  

  K_P2 = 7.3E-7
  S_P2 = 1.4E-3
    
  alpha= K_P2*b_mean/S_P2
  Impulse_P2=as.data.frame(Tab_inf$Time_sec)
  colnames(Impulse_P2)=c('Time_sec')
  Impulse_P2$d=sqrt(4*alpha*Impulse_P2$Time_sec)
  
  Impulse_P2$dZdt=b_mean/(2*S_P2)*(erf((L/2+X)/Impulse_P2$d)+erf((L/2-X)/Impulse_P2$d))*
    (erf((W/2+Y)/Impulse_P2$d)+erf((W/2-Y)/Impulse_P2$d))
  
  ### CONVOLUTION ### 
  
  for (day in (lim_P1+1):len){
    
    c=0
    i2=min(Impulse_duration,day-1)
    
    for (i_offset in (lim_P1+1):i2){ # Offset timestep
      
      if (i_offset<lim_P2){
        
        Inf_offset=Tab_inf$Inf_m_s[day-i_offset+lim_P1]
        
      }else{
        
        Inf_offset=0
        
      }
      
      c=c+Inf_offset*Impulse_P2$dZdt[i_offset]
    }
    
    Output_P2[day]= time_step*c
  }
  
  Output_P2=as.data.frame(Output_P2)
  h_P2=sqrt(as.data.frame(Output_P2)+b^2)
  h_P2["N_days"]=Tab_inf$N_days 
  
  Output_t=as.data.frame(Output_P2+Output_P1)
  h_t=sqrt(Output_t+b^2)
  h_t["N_days"]=Tab_inf$N_days
  colnames(h_t)=c("h","N_days")
  
  h_t
}

#########################################################################################################################################
#APPLICATION
#########################################################################################################################################

Mound_profile=as.data.frame(matrix(0, len_Tab,dim(profile_coord)[1]))

for (point in 1:dim(profile_coord)[1]){
  h <-Hantush_solution_TOTAL(c(profile_coord[point,'X'],profile_coord[point,'Y']))
  Mound_profile[,point]=h$h
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

