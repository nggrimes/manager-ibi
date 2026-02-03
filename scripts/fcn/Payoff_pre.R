Payoff_pre<-function(f,b,V,r,K,delta,ra,ut_mod,xgrid,p,c,cshape,trigger,premium,gamma,sim){

  xnext=map_dbl(sim,~growth(b=b,f=f,theta=.x,r=r,K=K))
  
  temp<-map_dbl(xnext,~spline(x=xgrid,y=V,xout=.x,method='natural')$y)
  
  pi_df<-map_dbl(sim,~pi(f=f,b=b,theta=.x,p=p,c=c,cshape=cshape,trigger=trigger,gamma=gamma,premium=premium))

  df<-data.frame(pi_out=pi_df,ut_mod=ut_mod,ra=ra)
  

  utility<-pmap_dbl(df,ut_post)
  
  eu=mean(utility)+mean(temp*delta,na.rm=TRUE)
  
  
  
  
  return(-eu)
}
