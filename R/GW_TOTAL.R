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

#########################################################################################################################################
#PRELUDE
#########################################################################################################################################



################
### PACKAGES ###      
################

#install.packages("pracma")
require(pracma)
tic()

#install.packages("emdbook")
require(emdbook)

#install.packages("xlsx")
require(xlsx) #Make sure Java version is compatible, if necessary install 64 bit

###################
### USER INPUTS ###      
###################

## Choose borewell from list

b_num=8
borewell_list=c('CH01','CH02','CH03','CH08','CH11','CH15','CH16','CH18')
borewell=borewell_list[b_num]

## Set parameters

L=115 # Length of basin
W=40 # Width of basin

lim_P1=60 # Limit phase 1
lim_P2=110 # Limit phase 2

Impulse_duration=180 # Duration of one impulse in days

########################### Lists of data ############################

# Hydraulic disconnection

X_list=c(50,3,215,247,219,177,167,211)
Y_list=c(87,66,305,300,334,293,328,346)
b_mean_list=c(8.4,8.8,1.0,1.6,1.1,3.6,10.0,4.1)
b_list=c(7.4,7.9,0.0,0.6,0.1,2.4,7.5,1)
fisslay_list=c(341.7,341.3,347.2,346.6,347.1,345.0,337.0,343.6)
P1_K_list=c(1.5E-4,1.5E-4,1.0E-3,3.0E-4,1.0E-3,2.0E-4,1.0E-4,9.0E-6)
P1_S_list=c(1.0E-2,1.0E-2,4.0E-2,1.0E-2,3.0E-2,2.0E-2,2.0E-2,1.0E-2)
  
X=X_list[b_num] # X coordinate from center of rectangle
Y=Y_list[b_num] # Y coordinate from center of rectangle
b_mean=b_mean_list[b_num] # Average aquifer thickness
b=b_list[b_num] # Original aquifer thickness
fiss_lay_depth=fisslay_list[b_num] # Fissured layer depth
P1_K=P1_K_list[b_num]
P1_S=P1_S_list[b_num]
  
# Hydraulic connection

distance_list=c(100,67,375,390,400,405,370,340)
niv_const_list=c(346.5,346.5,346.4,346.9,346.5,346.4,344,345)

#######################################################################

## Input data names and paths

dir_in="C:/Users/Madeleine/Documents/Doctorat/ARTICLE/Data_R/Inputs"
dir_out="C:/Users/Madeleine/Documents/Doctorat/ARTICLE/Data_R/Outputs"

# Hydraulic head
obs_file_name=paste("INPUTS_",borewell,"_daily.txt",sep="") # Infiltration 
obs_file_dir=paste(dir_in,"/",obs_file_name,sep="")
HH_obs<- read.csv(obs_file_dir,sep='\t', header=TRUE)
HH_obs$Date=as.Date(HH_obs$Date,format="%d/%m/%Y")
HH_obs$N_days=as.integer(HH_obs$Date-HH_obs$Date[1]+9)
HH_obs[paste("WL_",borewell,sep="")]=HH_obs[paste("HH_",borewell,sep="")]-fiss_lay_depth

# Infiltration data
inf_file_name="INPUTS_Infiltration.txt"
inf_file_dir=paste(dir_in,"/",inf_file_name,sep="")
Tab_inf <- read.csv(inf_file_dir,sep='\t', header=TRUE)
Tab_inf$Date=as.Date(Tab_inf$Date,format="%d/%m/%Y") # Time data
Tab_inf$N_days=as.integer(Tab_inf$Date-Tab_inf$Date[1])
Tab_inf$Time_sec=Tab_inf$N_days*24*60*60
Tab_inf$Inf_m_s=Tab_inf$Inf_mm_day/(1000*60*60*24)

# Basin water levels
bas_file_name="INPUTS_Basin_daily.txt"
bas_file_dir=paste(dir_in,"/",bas_file_name,sep="")
Tab_bas <- read.csv(bas_file_dir,sep='\t', header=TRUE)
Tab_bas$Date=as.Date(Tab_bas$Date,format="%d/%m/%Y") # Time data
Tab_bas$N_days=as.integer(Tab_bas$Date-Tab_bas$Date[1])

#########################################################################################################################################
#FUNCTION
#########################################################################################################################################

#########################################################################################################################################
########################################################   P1   #########################################################################
#########################################################################################################################################

len=length(Tab_inf$N_days)
Output_P1=matrix(0, len) # Initialize table

K=P1_K
S=P1_S

### IMPULSE RESPONSE ###      
  
alpha= K*b_mean/S
Impulse_P1=as.data.frame(Tab_inf$Time_sec)
colnames(Impulse_P1)=c('Time_sec')
Impulse_P1$d=sqrt(4*alpha*Impulse_P1$Time_sec)

Impulse_P1$dZdt=b_mean/(2*S)*(erf((L/2+X)/Impulse_P1$d)+erf((L/2-X)/Impulse_P1$d))*
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

P2_opt <- function(x){
  
  Output_P2=matrix(0, len) # Initialize table
  
  ### IMPULSE RESPONSE ###   
  
  S=x[1]
  K=x[2]
  
  alpha= K*b_mean/S
  Impulse_P2=as.data.frame(Tab_inf$Time_sec)
  colnames(Impulse_P2)=c('Time_sec')
  Impulse_P2$d=sqrt(4*alpha*Impulse_P2$Time_sec)
  
  Impulse_P2$dZdt=b_mean/(2*S)*(erf((L/2+X)/Impulse_P2$d)+erf((L/2-X)/Impulse_P2$d))*
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
  
  Output_t=as.data.frame(Output_P2+Output_P1)
  h_t=sqrt(Output_t+b^2)
  h_t["N_days"]=Tab_inf$N_days
  colnames(h_t)=c("h","N_days")

  Results_t=merge(h_t,HH_obs[c("N_days",paste("WL_",borewell,sep=""))], by=paste("N_days"),all.x = TRUE)
  Results_t=as.data.frame(c(Tab_inf["Date"],Results_t[,c(1,3,2)]))
  
  ### ERROR ###      
  
  RMSE_P2 <- sqrt(mean((Results_t[1:lim_P2,paste("WL_",borewell,sep="")]-Results_t$h[1:lim_P2])^2,na.rm=TRUE))

  c(RMSE_P2,Output_P2,as.data.frame(h_t$h),as.data.frame(Impulse_P2$dZdt))
  
}

#######################################################################

############################ Optimisation #############################

# Srange=c(seq(1,9) %o% 10^(-4:-2))[5:27] # Investigated range of parameters
# Krange=c(seq(1,9) %o% 10^(-8:-3))[5:46]
Srange=c(seq(1,9) %o% 10^(-4:-2)) # Investigated range of parameters
Krange=c(seq(1,9) %o% 10^(-7:-5))

Tab_RMSE_P2=matrix(0, length(Srange),length(Krange)) # Initialize tables

for (n in 1:length(Srange)){
  
  for (m in 1:length(Krange)){
    
    res=P2_opt(c(Srange[n],Krange[m]))

    Tab_RMSE_P2[n,m]=res[[1]]
    
  }
}

P2_K=Krange[which(Tab_RMSE_P2 == min(Tab_RMSE_P2), arr.ind = TRUE)[2]] # Optimal parameters
P2_S=Srange[which(Tab_RMSE_P2 == min(Tab_RMSE_P2), arr.ind = TRUE)[1]]

#######################################################################

### TEMPORARY OUTPUT ###  

Res_t=P2_opt(c(P2_S,P2_K))
h_t=as.data.frame(Res_t[[3]])
h_t["N_days"]=Tab_inf$N_days  

Output_P2=as.data.frame(Res_t[[2]])
h_P2=sqrt(as.data.frame(Output_P2)+b^2)
h_P2["N_days"]=Tab_inf$N_days  

Tab_res_t=merge(merge(h_P1,merge(h_t,h_P2,by="N_days"),by="N_days"),HH_obs[c("N_days",paste("WL_",borewell,sep=""))], by="N_days",all.x = TRUE)
Tab_res_t=as.data.frame(c(Tab_inf["Date"],Tab_res_t[,c(1,5,2,4,3)]))
colnames(Tab_res_t)=c("Date","N_days",paste("WL_",borewell,sep=""),"h_P1","h_P2","h_t")

Impulse_t=as.data.frame(Res_t[[4]])
colnames(Impulse_t)="Impulse_t"


#########################################################################################################################################
########################################################   P3   #########################################################################
#########################################################################################################################################

P3_opt <- function(x){

  D=x

  Output_P3=matrix(0, len) # Initialize table

  for (n in (lim_P2:len)-1){ #Loop over stream level variations
    for (k in len:(n+1)){ #Loop over all days of observed period

      if (n==lim_P2){
        h1=Tab_bas$HH_basin[lim_P2]
      }else{
        h1=Tab_bas$HH_basin[n]
      }

      h2=Tab_bas$HH_basin[n+1]
      Output_P3[k]=Output_P3[k] + (h2-h1)*erfc(distance_list[b_num]/sqrt(4*D*(k-n)*3600*24))

    }
  }

  Output_P3=as.data.frame(Output_P3)

  Output_t2=Output_P3
  Output_t2[lim_P2:len,]=Output_t2[lim_P2:len,]+Tab_res_t$h_t[lim_P2]
  # Output_t2[lim_P2:len,]=Output_t2[lim_P2:len,]+HH_obs[paste("WL_",borewell_list[b_num],sep="")][lim_P2,1]

  Results_t2=as.data.frame(Tab_bas$N_days)
  colnames(Results_t2)="N_days"
  Results_t2["Output_P3"]=Output_P3
  Results_t2["h_P3"]=Output_t2
  Results_t2=merge(Results_t2,HH_obs[c(4,5)],by="N_days",all.x = TRUE)

  RMSE_P3=sqrt(mean((Results_t2[lim_P2:len,paste("WL_",borewell,sep="")]-Results_t2$h_P3[lim_P2:len])^2,na.rm=TRUE))

  c(Results_t2,RMSE_P3)
}


# Drange=c(seq(1,9) %o% 10^(-4:0))[1:38]
# Tab_RMSE_P3=matrix(0, length(Drange))
# 
# for (i in 1:length(Drange)){
#   res=P3_opt(Drange[i])
#   Tab_RMSE_P3[i]=res[[5]]
# }
# 
# rownames(Tab_RMSE_P3)=Drange
# P3_D=Drange[which(Tab_RMSE_P3 == min(Tab_RMSE_P3), arr.ind = TRUE)][1]

P3_D=P1_K*b_mean/P1_S

Res_t2=P3_opt(P3_D)

Tab_res_t2=as.data.frame(Res_t2[c(1,3)])
Tab_res_t2=merge(Tab_res_t,Tab_res_t2, by="N_days",all.x = TRUE)

Tab_res_t2$h_t2=Tab_res_t2$h_t
Tab_res_t2$h_t2[lim_P2:len]=Tab_res_t2$h_P3[lim_P2:len]

Tab_res_t2=Tab_res_t2[c(2,1,3,4,5,7,6,8)]

#################
### SAVE DATA ###      
#################

write.xlsx(Tab_res_t2,paste(dir_out,"/",borewell,".xlsx",sep=""),sheetName="Output",row.names=FALSE)
params=t(as.data.frame(c(P2_K,b,P2_S,P3_D)))
colnames(params)=c("K","b","S","D")
write.xlsx(params,paste(dir_out,"/",borewell,".xlsx",sep=""),sheetName="Params",row.names=FALSE,append=TRUE)

############
### PLOT ###      
############

### ALL DATA ###

## Set up the plot

yrange <- c(min(range(Tab_res_t2[paste("WL_",borewell,sep="")],na.rm = TRUE)[1],range(Tab_res_t2["h_P2"])[1]),
            max(range(Tab_res_t2[paste("WL_",borewell,sep="")],na.rm = TRUE)[2],range(Tab_res_t2["h_P2"])[2]))
xrange <- range(Tab_res_t2$N_days)

nlines=ncol(Tab_res_t2)-2

plot(xrange, yrange, type="n", xlab="Time (days)",
     ylab="Water level (m)" )
colors <- rainbow(nlines)
linetype <- c(1:nlines)
plotchar <- seq(18,18+nlines,1) # Plotting symbol

## Add lines

for (li in 1:nlines) {
  lines(Tab_res_t2[,2], Tab_res_t2[,li+2], type="b", lwd=1.5,
        lty=linetype[li], col=colors[li], pch=plotchar[li])
}

### RELEVANT DATA ###

Tab_res_final=Tab_res_t2[c(1,2,3,8)]

#x11()
yrange <- c(min(range(Tab_res_final[paste("WL_",borewell,sep="")],na.rm = TRUE)[1],range(Tab_res_final["h_t2"],na.rm=TRUE)[1]),
            max(range(Tab_res_final[paste("WL_",borewell,sep="")],na.rm = TRUE)[2],range(Tab_res_final["h_t2"],na.rm=TRUE)[2]))
xrange <- range(Tab_res_final$N_days)

nlines=ncol(Tab_res_final)-2

plot(xrange, yrange, type="n", xlab="Time (days)",
     ylab="Water level (m)" )
colors <- rainbow(nlines)
linetype <- c(1:nlines)
plotchar <- seq(18,18+nlines,1)

for (li in 1:nlines) {
  lines(Tab_res_final[,2], Tab_res_final[,li+2], type="b", lwd=1.5,
        lty=linetype[li], col=colors[li], pch=plotchar[li])
}

# beepr::beep(2)
toc()