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

a=0.0008

#Environmental parameters
mean_theta=0

sigma_theta=0.2


#### HCR control ####
hcr="two"

#### Payouts and Premiums parameters ####
insurance=1 #use insurance 1, 0 no insurance
gamma=250
coverage=0
findNearest<-function(x,y){
  index<-which.min(abs(x-y))
  return(index)
}




#### Set up state grid ####
xgrid = seq(small,K,length.out=sizex)

#### Set up shock grid ####
#eps_list<-c(-.1,0)
#theta_list<-c(.5,-.1)
#d=c(theta_list,eps_list)
#z_list=split(d,ceiling(seq_along(d)/2))

#z_grid=c(0.4,0.6)

n_z=13
n_sims=100000

set.seed(42)

#Get distribution and bins of
theta_raw<-rnorm(n_sims,mean=mean_theta,sd=sigma_theta)
theta_dist<-theta_raw[which(theta_raw>-1)]

thetahist<-hist(theta_dist,breaks=seq(min(theta_dist),max(theta_dist),l=n_z+1))
## Make bins classification
bins=theta_dist
for(k in 2:length(thetahist$breaks)){
  index=which(theta_dist>thetahist$breaks[k-1]&theta_dist<=thetahist$breaks[k])
  bins[index]=k-1
}

theta_df<-data.frame(theta_dist=theta_dist,bins=bins)

means<- theta_df %>% 
  group_by(bins) %>% 
  summarize(mean=mean(theta_dist))

#elminate first observation
means=means[-1,]
theta_prop=thetahist$counts/n_sims
theta_val=means$mean

# Use distribution to get insurance parameters
trigger_index=findNearest(coverage,theta_val)
trigger=theta_val[trigger_index]
trigger_prop=sum(theta_prop[1:trigger_index])


#### Optimization functions ####
harvest <- function(f, biomass, theta){
  if(hcr == "one"){
    #1)Harvest set to according to pre-shock biomass
    harvest = f * biomass
  }else if(hcr == "two"){ 
    #2)estimate of post-shock biomass (no epsilon or beta)
    harvest = f * biomass * (1+theta)
  }else if(hcr == "three"){ 
    #3)perfect knowledge of post-shcok biomass (w/epsilon and beta)
    #harvest = f * biomass * (beta * exp((theta - sigma_theta^2/2))+(1-beta)*exp(epsilon - sigma_eps^2/2))
  }
  
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
  
  #Make Calculations with insurance, insurance defined globally in param space
  if(insurance==1){
    if(theta<=trigger){ #check if we payout, note: constant payment for now so premium and payout easy to calculate
      payout=gamma
    }else{
      payout=0
    }
    
    premium=gamma*trigger_prop
    profit=(p*harvest(f=f,biomass=b,theta = theta))-(c*(f^cshape))+payout-premium
  }else{
    
    profit = (p*harvest(f=f,biomass=b,theta = theta))-(c*(f^cshape))
  }
  
  return(profit)
}

Payoff = function(f,b,theta,V){

  xnext = growth(b,f,theta)
  
  sp_store<-vector(mode="numeric",length=ncol(V))
  

  
  for(n in 1:ncol(V)){
    
    temp<-spline(x=xgrid,y=V[,n],xout=xnext,method="natural")
    sp_store[n]<-temp$y
    
    
  }
  
  Vfuture=delta*sp_store*theta_prop
  #Discounts and add shock probability weights for all splined future Vs
  temp_prof=pi(f,b,theta)
  

  utility=-(1-exp(-a*(temp_prof))+sum(Vfuture))
  
  
  #browser()
  return(utility)
}

options=list("algorithm"="NLOPT_LN_COBYLA",xtol_rel=1e-06)
DFall = data.frame()
Vnext = matrix(ncol = length(theta_prop),nrow =sizex)
V = matrix(0,ncol=length(theta_prop),nrow=sizex)
pol_converged=matrix(0,nrow=sizex,ncol = length(theta_prop)) #Harvest Policy function
polnow=matrix(0,nrow=sizex,ncol = length(theta_prop))

step=1
error="False"
tol=0.001

while(error=="False"){
  t=step
  print(t)
  #Loop over shocks space
  for(j in 1:length(theta_prop)){
    
    temptheta=theta_val[j]
    
    
    for(i in 1:length(xgrid)){
      b=xgrid[i]
      guess=0.1
      low=.0001
      high= 1
      output = nloptr(x0=guess,eval_f=Payoff,lb=0,ub=1,opts=options,b=b,V=V,theta=temptheta)
      fstar = output$solution
      Vstar = -output$objective
      Vnext[i,j]=Vstar
      DFnow = data.frame(time=t,b=b,fstar=fstar,Vstar=Vstar,theta=temptheta)
      DFall = bind_rows(DFall,DFnow)
      polnow[i,j]=fstar
    }
  }
  
  #Check if errors reach tolerance for policy function convergence
  
 
  comp<-abs(polnow-pol_converged)
  error_vec<-colMeans(comp)
  
  
  if(sum(error_vec<tol)==length(error_vec) & step>10){
    error="True"
  }else{
    error="False"
  }
  
  step=step+1
 
  V=Vnext
  pol_converged=polnow
  
  
  #Check to make sure we stop eventually
  if(step==400){
    print('Maximum number of iterations (400) reached.')
    break
  }
}



temp_ins_t<-DFall %>% 
  filter(time==max(DFall$time)) %>% 
  group_by(theta) %>% 
  summarize(pol_opt=mean(fstar)) %>% 
  mutate(ins='Insurance')

temp<-DFall %>% 
  filter(time==max(DFall$time)) %>% 
  group_by(b) %>% 
  summarize(pol_opt=weighted.mean(fstar,theta_prop)) %>% 
  mutate(ins='Insurance')



b_df<-rbind(no_pol_conv_b,temp)
t_df<-rbind(no_pol_conv,temp_ins_t)



ggplot(b_df,aes(x=b,y=pol_opt,color=ins))+
  geom_line()

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


