#########################################################################################################################################
#
# Author: Madeleine NICOLAS       
#         madeleine.nicolas@univ-rennes1.fr
# 
# Purpose: Functions for solution and optimization of basin infiltration analytical solutions
#          
# Date:   August 2017
#
#########################################################################################################################################

Hantush_optimisation <- function(x){
  
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
  
  FINAL=merge(h,h_obs[c(paste(borewell,"_Days",sep=""),paste(borewell,"_WL",sep=""))], by=paste(borewell,"_Days",sep=""),all.x = TRUE)
  
  FINAL=FINAL[,c(1,3,2)]
  
  #############
  ### ERROR ###      
  #############
  
  RMSE <- sqrt(mean((FINAL[paste(borewell,"_WL",sep="")]-FINAL$Output)^2,na.rm=TRUE))
  
  RMSE_1 <- sqrt(mean((FINAL[1:50,paste(borewell,"_WL",sep="")]-FINAL$Output[1:50])^2,na.rm=TRUE))
  
  RMSE_2 <- sqrt(mean((FINAL[-(1:80),paste(borewell,"_WL",sep="")]-FINAL$Output[-(1:80)])^2,na.rm=TRUE))
  
  c(RMSE,RMSE_1,RMSE_2)
  
}


Hantush_solution <- function(x){
  
  S=x[1]
  K=x[2]
  nfit=x[3]
  
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
      
      Output[i]=time_step*c
      
    }
  }
  
  Output=as.data.frame(Output)
  h=sqrt(Output+b^2)
  h[paste(borewell,"_Days",sep="")]=Tab_inf$Time_days
  
  FINAL=merge(h,h_obs[c(paste(borewell,"_Days",sep=""),paste(borewell,"_WL",sep=""))], by=paste(borewell,"_Days",sep=""),all.x = TRUE)
  FINAL=as.data.frame(c(Tab_inf["Date"],FINAL[,c(1,3,2)]))
  
  #############
  ### ERROR ###      
  #############
  
  RMSE <- sqrt(mean((FINAL[paste(borewell,"_WL",sep="")]-FINAL$Output)^2,na.rm=TRUE))
  
  ############
  ### SAVE ###      
  ############
  
  write.xlsx(FINAL,file=paste(dir_out,"/",borewell_list[b_num],"_FINAL_",fit[nfit],".xlsx",sep=""),sheetName="Output",row.names=FALSE)
  write.xlsx(RMSE,file=paste(dir_out,"/",borewell_list[b_num],"_FINAL_",fit[nfit],".xlsx",sep=""),sheetName="RMSE",row.names=FALSE,append=TRUE)
  
}