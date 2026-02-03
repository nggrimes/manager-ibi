### Run a forward simulation using the harvest control rules

library(tidyverse)

load(here::here('data','ins_pol_unknown.Rdata'))

load(here::here('data','no_pol_unknown.Rdata'))

load(here::here('data','rn_pol_2_22.Rdata'))

# make 30 year random weather draws 500 times
set.seed(42)


rnd_weather <- replicate(500, rnorm(30, 0, 0.5))

#replace any values of rnd_Weather less than -1 with -0.99

rnd_weather[rnd_weather < -1] <- -0.99


### Parameters
cshape=1.3 #cost shaping parameter
delta=.96

r=0.8
K=1000

small=K/100
p=5
c=250


ut_mod='cara'
a=0.008 

insurance=1 #use insurance 1, 0 no insurance
gamma=250
trigger=0

mean_theta=0

sigma_theta=0.5
premium=gamma*pnorm(trigger,mean_theta,sigma_theta)


# define growth, profit, and utility functions

growth = function(b,f,theta)
{
  
  escapement=b*(1+theta)-harvest(f=f,biomass=b,theta = theta)
  
  if(escapement>K){
    bnext=K
  }else if(escapement<0){
    bnext=small
  }  else{
    
    bnext = escapement + r*escapement*(1-escapement/K)
  }
  
  return(bnext)
}

harvest <- function(f, biomass, theta){
  harvest = f * biomass * (1+theta)
  
  
  return(harvest)
}

pi = function(f,b,theta)
{
  
  profit = (p*harvest(f=f,biomass=b,theta = theta))-(c*(f^cshape))
  
  return(profit)
}

ut<-function(f,b,theta,trigger,pi,payout,premium,method='cara'){
  
  
  if(theta<trigger){
    pi_out<-pi+payout-premium
  } else{
    pi_out<-pi-premium
  }
  
  
  if(method=='cara'){
    out=1-exp(-a*(pi_out))
  }else if(method=='log'){
    out=log(pi_out+1000)  #Plus 1000 is to make sure we dont have negative values
  }else if (method=='rn'){
    out=pi_out
  } else{
    print('Not a valid utility function. Please use cara: for constant absolute risk aversion, log: for log utility, or rn: for risk neutral')
  }
  
  return(out)
}

b0=500

final_ut<-rep(NA,500)
final_b<-rep(NA,500)

final_ut_no<-rep(NA,500)
final_b_no<-rep(NA,500)

final_ut_rn<-rep(NA,500)
final_b_rn<-rep(NA,500)

for(j in 1:ncol(rnd_weather)){
  
  # make a data frame to store the results
  sim <- data.frame(year = 1:30, 
                    ins_ut = rep(NA, 30), 
                    no_ut = rep(NA, 30),
                    ins_b= rep(NA, 30),
                    no_b = rep(NA, 30),
                    rn_ut = rep(NA, 30),
                    rn_b = rep(NA, 30))
  
  # set the initial population size
  sim$ins_b[1] <- b0
  sim$no_b[1] <- b0
  sim$rn_b[1] <- b0
  
  f_ins<-spline(x=ins_pol_conv$b,y=ins_pol_conv$pol_opt,xout=sim$ins_b[1])$y
  f_no<-spline(x=no_pol_conv$b,y=no_pol_conv$pol_opt,xout=sim$no_b[1])$y
  f_rn<-spline(x=no_pol_conv$b,y=no_pol_conv$pol_opt,xout=sim$rn_b[1])$y
  
  # get insurance utility
  temp_pi<-pi(f_ins,sim$ins_b[1],rnd_weather[1,j])
  
  sim$ins_ut[1]<-ut(f_ins,sim$ins_b[1],rnd_weather[1,j],trigger,temp_pi,gamma,premium,ut_mod)
  
  # get no insurance utility
  temp_pi<-pi(f_no,sim$no_b[1],rnd_weather[1,j])
  
  sim$no_ut[1]<-ut(f_no,sim$no_b[1],rnd_weather[1,j],trigger,temp_pi,0,0,ut_mod)
  
  # get risk neutral utility
  temp_pi<-pi(f_rn,sim$rn_b[1],rnd_weather[1,j])
  sim$rn_ut[1]<-ut(f_rn,sim$rn_b[1],rnd_weather[1,j],trigger,temp_pi,gamma,premium,ut_mod)
  
  # get the next year's population size
  sim$ins_b[1] <- growth(sim$ins_b[1],f_ins,rnd_weather[1,j])
  sim$no_b[1] <- growth(sim$no_b[1],f_no,rnd_weather[1,j])
  sim$rn_b[1] <- growth(sim$rn_b[1],f_rn,rnd_weather[1,j])
  
  # loop over the years
  for(i in 2:30){
    sim$ins_b[i] <- growth(sim$ins_b[i-1],f_ins,rnd_weather[i,j])
    sim$no_b[i] <- growth(sim$no_b[i-1],f_no,rnd_weather[i,j])
    sim$rn_b[i] <- growth(sim$rn_b[i-1],f_rn,rnd_weather[i,j])
    
    
    f_ins<-spline(x=ins_pol_conv$b,y=ins_pol_conv$pol_opt,xout=sim$ins_b[i])$y
    f_no<-spline(x=no_pol_conv$b,y=no_pol_conv$pol_opt,xout=sim$no_b[i])$y
    f_rn<-spline(x=no_pol_conv$b,y=no_pol_conv$pol_opt,xout=sim$rn_b[i])$y
    
    # get insurance utility
    temp_pi<-pi(f_ins,sim$ins_b[i],rnd_weather[i,j])
    sim$ins_ut[i]<-ut(f_ins,sim$ins_b[i],rnd_weather[i,j],trigger,temp_pi,gamma,premium,ut_mod)
    
    # get no insurance utility
    temp_pi<-pi(f_no,sim$no_b[i],rnd_weather[i,j])
    
    sim$no_ut[i]<-ut(f_no,sim$no_b[i],rnd_weather[i,j],trigger,temp_pi,0,0,ut_mod)
    
    # get risk neutral utility
    temp_pi<-pi(f_rn,sim$rn_b[i],rnd_weather[i,j])
    sim$rn_ut[i]<-ut(f_rn,sim$rn_b[i],rnd_weather[i,j],trigger,temp_pi,gamma,premium,ut_mod)
    
    # get the next year's population size
   
    
  }
  
  # discount the utilities
  sim$ins_ut <- sim$ins_ut * delta^(0:29)
  sim$no_ut <- sim$no_ut * delta^(0:29)
  sim$rn_ut <- sim$rn_ut * delta^(0:29)
  
  # store the final biomass and summed discounted utilities in storage vec
  
  final_ut[j] <- sum(sim$ins_ut)
  final_b[j] <- sim$ins_b[30]
  
  final_ut_no[j] <- sum(sim$no_ut)
  
  final_b_no[j] <- sim$no_b[30]
  
  final_ut_rn[j] <- sum(sim$rn_ut)
  final_b_rn[j] <- sim$rn_b[30]
  
  
  
}

df<-data.frame(final_ut, final_b, final_ut_no, final_b_no,final_ut_rn,final_b_rn,sim=1:500,b0=b0) |> 
  drop_na() |> 
  pivot_longer(cols=c(final_ut, final_b, final_ut_no, final_b_no,final_b_rn,final_ut_rn), names_to='variable', values_to='value')

df |> filter(variable %in% c('final_ut', 'final_ut_no','final_ut_rn')) |> 
  ggplot(aes(x=value, fill=variable,color=variable)) +
  geom_density(alpha=0.33,linewidth=2.5) +
  theme_minimal() +
  scale_color_manual(name="",values=c("#003660","#79A540","red"),guide='none')+
  scale_fill_manual(name="",values=c("#003660","#79A540","red"),labels=c('final_ut'="Insurance",'final_ut_no'="No Insurance",'final_ut_rn'="No Ins Policy Fcn\n w/ ins"))+
  theme_classic()+
  theme(legend.text=element_text(size=24))+
  theme(legend.title =element_text(size=28))+
  theme(axis.text =element_text(size=22))+
  theme(axis.title = element_text(size=26))+
  labs(title='', x='Sum of discounted utility', y='Density')+
  scale_y_continuous(expand = c(0,0))+
  guides(fill = guide_legend(byrow=TRUE, override.aes = list(shape = NA, color = NA, fill = c("#003660", "#79A540","red"))))+
  theme(legend.spacing.y = unit(10, "cm"))




df |> filter(variable %in% c('final_b_rn', 'final_b_no')) |> 
  ggplot(aes(x=value, fill=variable, color=variable)) +
  geom_density(aes(y = after_stat(density * length(df$value))),alpha=0.33,linewidth=2.5) +
  theme_minimal() +
  scale_color_manual(name="",values=c("#003660","#79A540"),guide='none')+
  scale_fill_manual(name="",values=c("#003660","#79A540"),labels=c("Insurance","No Insurance"))+
  theme_classic()+
  theme(legend.text=element_text(size=24))+
  theme(legend.title =element_text(size=28))+
  theme(axis.text =element_text(size=22))+
  theme(axis.title = element_text(size=26))+
  labs(title='', x='Final Biomass', y='Density')+
  scale_y_continuous(expand = c(0,0))+
  guides(fill = guide_legend(byrow=TRUE, override.aes = list(shape = NA, color = NA, fill = c("#003660", "#79A540"))))+
  theme(legend.spacing.y = unit(10, "cm"))

# Run a t.test to compare if the distributions are different
t.test(final_ut ~ final_ut_no, data=df |> pivot_wider(names_from='variable', values_from='value'))




plot_df<-rbind(no_pol_conv %>% 
                 mutate(ins='No Insurance'),
               ins_pol_conv %>% 
                 mutate(ins='Insurance'),
               rn_pol_conv %>% 
                 mutate(ins='Risk Neutral'))




ggplot(plot_df |> filter(ins!='Risk Neutral'),aes(x=b,y=pol_opt,color=ins))+
  geom_line(linewidth=2.5)+
  scale_color_manual(name="",values=c("#003660","#79A540"))+
  theme_classic()+
  labs(y='Optimal Harvest Rate',x='Biomass')+
  theme(legend.text=element_text(size=24))+
  theme(legend.title =element_text(size=28))+
  theme(axis.text =element_text(size=22))+
  theme(axis.title = element_text(size=26))


ggplot(plot_df |> filter(ins!='Insurance'),aes(x=b,y=pol_opt,color=ins))+
  geom_line(linewidth=2.5)+
  scale_color_manual(name="",values=c("#79A540","#900C3F"))+
  theme_classic()+
  labs(y='Optimal Harvest Rate',x='Biomass')+
  theme(legend.text=element_text(size=24))+
  theme(legend.title =element_text(size=28))+
  theme(axis.text =element_text(size=22))+
  theme(axis.title = element_text(size=26))

ggplot(plot_df |> filter(ins=='Risk Neutral'),aes(x=b,y=pol_opt,color=ins))+
  geom_line(linewidth=2.5)+
  scale_color_manual(name="",values=c("#900C3F"))+
  theme_classic()+
  labs(y='Optimal Harvest Rate',x='Biomass')+
  theme(legend.text=element_text(size=24))+
  theme(legend.title =element_text(size=28))+
  theme(axis.text =element_text(size=22))+
  theme(axis.title = element_text(size=26))

