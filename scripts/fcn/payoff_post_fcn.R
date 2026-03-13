payoff_post_fcn<-function(q,b,shock,V,pay,prem_rate,gamma,mod,c_cost,w_probs,r,K,biomass_grid,w_grid,beta,trigger){


  if(pay==0){
    pi<-q-c_cost*q/((1+shock)*b)
    
    e_ut<-utility(pi,gamma=gamma,boost=0,mod=mod)
  } else{
    
    # Fishers buy insurance based on expected value
    fish_obj<-optim(par=b/10,
                    fn=fisher_ut,
                    b=b,
                    q=q,
                    prem_rate=prem_rate,
                    shock=w_grid,
                    gamma=gamma,
                    c_cost=c_cost,
                    w_probs=w_probs,
                    mod=mod,
                    method='Brent',
                    trigger=trigger,
                    lower=0,
                    upper=K)
    #the actual the recieve with the known shock to the manager
    payout<-insurance_payout(shock,fish_obj$par,trigger)
    premium<-fish_obj$par*prem_rate
    
    pi=q-c_cost*q/((1+shock)*b)+payout-premium
    
    e_ut<-utility(pi,
                  gamma = gamma,
                  boost=0,
                  mod=mod)
  }
  

  
  esc<-max((1+shock)*b-q,K/1000)
  
  
  xnext<-growth(esc,r=r,K=K)
  
  sp_store<-vector(mode="numeric",length=ncol(V))
  
  
  
  for(n in 1:ncol(V)){
    
    temp<-spline(x=biomass_grid,y=V[,n],xout=xnext,method="natural")
    sp_store[n]<-temp$y
    
    
  }
  
  Vfuture=beta*sum(sp_store*w_probs)
  
  return(-(e_ut+beta*Vfuture))
  
}
