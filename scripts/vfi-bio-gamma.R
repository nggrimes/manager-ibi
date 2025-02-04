library(tidyverse)

sizex = 50 #size of the state grid
sizez = 4  #size of the shock grid
T=50 #time horizon for backward induction

cshape=1.3 #cost shaping parameter
delta=.96

r=0.8
K=1000

small=K/100
p=5
c=250


ut_mod='cara'
a=0.008 #risk aversion

#Environmental parameters
mean_theta=0

sigma_theta=0.5

set.seed(42)
sim<-rnorm(1000,mean=mean_theta,sd=sigma_theta)

sim<-sim[which(sim>=-1)]


#### HCR control ####
hcr="two"

#### Payouts and Premiums parameters ####
insurance=1 #use insurance 1, 0 no insurance
gamma=0.5
trigger=-.33

premium=gamma*pnorm(trigger,mean_theta,sigma_theta)





#### Set up state grid ####
xgrid = seq(small,K,length.out=sizex)


#### Optimization functions ####
harvest <- function(f, biomass, theta){
    harvest = f * biomass * (1+theta)

  
  return(harvest)
}

#Growth function
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
    out=log(pi_out+1000)  #Plus 2000 is to make sure we dont have negative values
  }else if (method=='rn'){
    out=pi_out
  } else{
    print('Not a valid utility function. Please use cara: for constant absolute risk aversion, log: for log utility, or rn: for risk neutral')
  }
  
  return(out)
}


Payoff = function(f,b,sim,V)
{
  
  payout=gamma*pi(f,b,theta=0)
  
  premium=payout*pnorm(trigger,mean_theta,sigma_theta)
 
  
  xnext=map_dbl(sim,~growth(b=b,f=f,theta=.x))
  
  temp<-map_dbl(xnext,~spline(x=xgrid,y=V,xout=.x,method='natural')$y)
  
  pi_df=pi(f=f,b=b,theta=sim)
  
  df<-data.frame(f=f,b=b,theta=sim,trigger=trigger,pi=pi_df,payout=payout,premium=premium,method=ut_mod)
  utility<-pmap_dbl(df,ut)
  
  eu=mean(utility+temp*delta,na.rm=TRUE)
  #xnext = max(small,f(h,x))

  
  #browser()
  return(-eu)
}

DFall = data.frame()
Vnext = vector(mode='numeric',length=sizex)
V = rep(0,length.out=sizex)
pol_converged=vector(mode='numeric',length=sizex) #Harvest Policy function
polnow=vector(mode='numeric',length=sizex)

step=1
error="False"
tol=0.001

while(error=="False"){
  t=step
  print(t)
  #Loop over shocks space
    
    for(i in length(xgrid):1){
      b=xgrid[i]
      guess=0.1
      low=.0001
      high= 1
      output = optim(par=guess,fn=Payoff,lower=low,upper=high,b=b,V=V,sim=sim,method='L-BFGS-B')
      fstar = output$par
      Vstar = -output$value
      Vnext[i]=Vstar
      DFnow = data.frame(time=t,b=b,fstar=fstar,Vstar=Vstar)
      DFall = bind_rows(DFall,DFnow)
      polnow[i]=fstar
    }
  
  
  #Check if errors reach tolerance for policy function convergence
  comp<-abs(polnow-pol_converged)
  error_vec<-mean(comp)
  
  
  if(error_vec<tol & step>9){
    error="True"
  }else{
    error="False"
  }
  
  step=step+1
  #browser()
  V=Vnext
  pol_converged=polnow
  
  
  #Check to make sure we stop eventually
  if(step==400){
    print('Maximum number of iterations (400) reached.')
    break
  }
}



ins_pol_conv<-DFall %>% 
  filter(time==max(DFall$time)) %>% 
  group_by(b) %>% 
  summarize(pol_opt=mean(fstar))

ggplot(ins_pol_conv,aes(x=b,y=fstar))+
  geom_line()+
  geom_line(data=no_pol_conv,aes(x=b,y=pol_opt),color='blue')
  

ggplot(rn_pol_conv,aes(x=b,y=pol_opt))+
  geom_line(aes(color='Risk Neutral'),linewidth=3)+
  geom_line(data=ins_pol_conv,aes(x=b,y=pol_opt,color='Insurance'),linewidth=3)+
  geom_line(data=noi_pol_conv,aes(x=b,y=pol_opt,color='Risk Averse'),linewidth=3)+
  scale_color_manual(name="Risk Preferences",values=c('Risk Averse'="blue",'Insurance'="red",'Risk Neutral'="forestgreen"))+
  theme_minimal()+
  labs(y='Optimal Harvest',x='biomass')


combo_df<-cbind(rn_pol_conv,ins_pol_conv$pol_opt,noi_pol_conv$pol_opt)

colnames(combo_df)<-c("b","neutral","insurance","averse")

combo_df_long<-pivot_longer(combo_df,cols=!b,names_to="model",values_to = "f_opt")

combo_df_long[1:3,3]<-0

ggplot(combo_df_long,aes(x=b,y=f_opt,color=model))+
  geom_line(linewidth=1.5)+
  scale_color_manual(name="Risk Preferences",labels=c("Risk Neutral","Insurance","Risk Averse"),values=c("blue","red","forestgreen"))+
  theme_classic()+
  labs(y='Optimal Harvest',x='Biomass')


