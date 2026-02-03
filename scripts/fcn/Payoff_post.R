Payoff_post = function(f,b,theta,V,r,K,delta,theta_prop,ra,ut_mod,xgrid,p,c,cshape,trigger,premium,gamma){
  
  xnext = growth(b,f,theta,r,K)
  
  sp_store<-vector(mode="numeric",length=ncol(V))
  
  
  
  for(n in 1:ncol(V)){
    
    temp<-spline(x=xgrid,y=V[,n],xout=xnext,method="natural")
    sp_store[n]<-temp$y
    
    
  }
  
  Vfuture=delta*sp_store*theta_prop
  #Discounts and add shock probability weights for all splined future Vs
  temp_prof=pi(f,b,theta,p,c,cshape,trigger,premium,gamma)
  
  
  utility=-(ut_post(temp_prof,ra,ut_mod)+sum(Vfuture))
  
  
  #browser()
  return(utility)
}
