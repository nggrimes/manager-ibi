payoff_fcn<-function(q,b,shock,V,pay,prem_rate,gamma,mod,c_cost,w_probs,r,K,biomass_grid,w_grid,beta,trigger){
  
  
  if(pay==0){
    pi<-map_dbl(shock,~q-c_cost*q/((1+.x)*b))
    
    e_ut<-sum(map_dbl(pi,~utility(.x,gamma=gamma,boost=0,mod=mod))*w_probs)
  } else{
    fish_obj<-optim(par=b/10,
                    fn=fisher_ut,
                    b=b,
                    q=q,
                    prem_rate=prem_rate,
                    shock=shock,
                    gamma=gamma,
                    c_cost=c_cost,
                    w_probs=w_probs,
                    mod=mod,
                    method='Brent',
                    trigger=trigger,
                    lower=0,
                    upper=K)
    
    e_ut=-fish_obj$value
  }
  
  
  
  esc<-map2_dbl(.x=q,.y=w_grid,~max((1+.y)*b-.x,K/1000))
  
  
  xnext<-map_dbl(esc,~growth(.x,r=r,K=K))
  
  temp<-map_dbl(xnext,~spline(x=biomass_grid,y=V,xout=.x,method='natural')$y)
  
  Vfuture=sum(temp*w_probs)
  
  return(-(e_ut+beta*Vfuture))
  
}
