#########################################################################################################################################
#
# Author: Madeleine NICOLAS       
#         madeleine.nicolas@univ-rennes1.fr
# 
# Purpose: Implementation of Hantush analytical solution (1967) to model recharge from rectangular infiltration basin 
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

#install.packages("pracma")
require(pracma)

#install.packages("xlsx")
require(xlsx) #Make sure Java version is compatible, if necessary install 64 bit

###################
### USER INPUTS ###      
###################

## Choose borewell from list

b_num=1
borewell_list=c('CH01','CH02','CH03','CH08','CH11','CH15','CH16','CH18')
borewell=borewell_list[b_num]

## Choose period to fit over 

nfit=1
fit=c("1","2","all")

## Set parameters

# K=eval(as.name(paste("K_opt_",fit[nfit],sep="")))
# S=eval(as.name(paste("S_opt_",fit[nfit],sep="")))

K=1.5e-4 # Hydraulic conductivity (m/s)
S=1e-2 # Storativity (-)

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

# X=X_list[b_num] # X coordinate from center of rectangle
# Y=Y_list[b_num] # Y coordinate from center of rectangle
b_mean=b_mean_list[b_num] # Average aquifer thickness
b=b_list[b_num] # Original aquifer thickness
Fiss_lay_depth=fisslay_list[b_num] # Fissured layer depth

###############################################################

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

########################
### IMPULSE RESPONSE ###      
########################

## Groundwater level response to R=1 (impulse response to be multiplied by actual recharge rate)

alpha= K*b_mean/S
Impulse=as.data.frame(Tab_inf$Time_sec)
colnames(Impulse)=c('Time_sec')
Impulse$d=sqrt(4*alpha*Impulse$Time_sec)

Impulse$dZdt=b_mean/(2*S)*(erf((L/2+X)/Impulse$d)+erf((L/2-X)/Impulse$d))*
  (erf((W/2+Y)/Impulse$d)+erf((W/2-Y)/Impulse$d))
  
#plot(Impulse$Time_sec,Impulse$dZdt,xlab='Time (sec)', ylab='Impulse response (m)')

###################
### CONVOLUTION ###      
###################

len_tab=length(Impulse$Time_sec)
time_step=Tab_inf$Time_sec[2]
       
Output=rep(0, len_tab)

for (i1 in 2:len_tab){
  
  c=0
  
  i2=min(Impulse_duration,i1-1)
  
  for (i_offset in 1:i2){ # Offset timestep
    
    c=c+Tab_inf$Inf_m_s[i1-i_offset]*Impulse$dZdt[i_offset]
    
  }
  Output[i1]=time_step*c
}

Output=as.data.frame(Output)
h=sqrt(Output+b^2)
h[paste(borewell,"_Days",sep="")]=Tab_inf$Time_days

Results=merge(h,HH_obs[c(paste(borewell,"_Days",sep=""),paste(borewell,"_WL",sep=""))], by=paste(borewell,"_Days",sep=""),all.x = TRUE)

Results=as.data.frame(c(Tab_inf["Date"],Results[,c(1,3,2)]))

#############
### ERROR ###      
#############

RMSE <- sqrt(mean((Results[paste(borewell,"_WL",sep="")]-Results$Output)^2,na.rm=TRUE))
RMSE_2 <- sqrt(mean((Results[-(1:80),paste(borewell,"_WL",sep="")]-Results$Output[-(1:80)])^2,na.rm=TRUE))

############
### SAVE ###      
############

# write.xlsx(Results,file=paste(dir_out,"/",borewell_list[b_num],"_FINAL_",fit[nfit],".xlsx",sep=""),sheetName="Output",row.names=FALSE)
# write.xlsx(RMSE,file=paste(dir_out,"/",borewell_list[b_num],"_FINAL_",fit[nfit],".xlsx",sep=""),sheetName="RMSE",row.names=FALSE,append=TRUE)

Results=Results[,c(2,3,4)]

############
### PLOT ###      
############

# Create Line Chart

yrange <- c(min(range(Results[paste(borewell,"_WL",sep="")],na.rm = TRUE)[1],range(Results["Output"])[1]),
            max(range(Results[paste(borewell,"_WL",sep="")],na.rm = TRUE)[2],range(Results["Output"])[2]))

xrange <- range(Results[paste(borewell,"_Days",sep="")])

nlines=2

# Set up the plot

plot(xrange, yrange, type="n", xlab="Time (days)",
     ylab="Water level (m)" )
colors <- rainbow(nlines)
linetype <- c(1:nlines)
plotchar <- seq(18,18+nlines,1) # Plotting symbol

# Add lines

for (li in 1:nlines) {
  HH <- Results[c(1,li)]
  lines(Results[,1], Results[,li+1], type="b", lwd=1.5,
        lty=linetype[li], col=colors[li], pch=plotchar[li])
}

# Add a title and subtitle

#title("Analytical groundwater mound solution","Hantush 1967")

# Add a legend

# legend(180, 6, c('Observed','Analytical solution'), cex=0.6, col=colors,
#        pch=plotchar, lty=linetype, title=NULL)
