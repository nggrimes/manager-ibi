fisher_ut<-function(cov,b,q,prem_rate,shock,gamma,mod='cara',c_cost,w_probs,trigger){
  
  payout<-map_dbl(.x=shock,
                  .f=~insurance_payout(shock=.x,
                                       coverage=cov,
                                       trigger=trigger))
  premium<-cov*prem_rate
  
  pi<-map2_dbl(.x=shock,.y=payout,~q-c_cost*q/((1+.x)*b)+.y-premium)
  
  e_ut<-map_dbl(pi,~utility(.x,gamma,boost=0,mod=mod))*w_probs
  
  return(-sum(e_ut))
}
