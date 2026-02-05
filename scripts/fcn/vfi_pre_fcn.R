vfi_pre_fcn<-function(sizex,n_z,cshape,delta,
                      r,K,p,c,a,
                      mean_theta,sigma_theta,
                      gamma,trigger,m,
                      ut_mod){
  
  

  
  small=K/100
  
  set.seed(42)
  sim<-rnorm(1000,mean=mean_theta,sd=sigma_theta)
  
  sim<-sim[which(sim>=-1)]

  premium_rate=mean(map_dbl(trigger-sim,~max(0,.x)))*(1+m) 
  
  
  
  
  
  #### Set up state grid ####
  xgrid = seq(small,K,length.out=sizex)
  
  DFall = data.frame()
  Vnext = vector(mode='numeric',length=sizex)
  V = rep(0,length.out=sizex)
  pol_converged=vector(mode='numeric',length=sizex) #Harvest Policy function
  polnow=vector(mode='numeric',length=sizex)
  
  step=1
  error="False"
  tol=0.0001
  
  #options=list("algorithm"="NLOPT_LN_COBYLA",xtol_rel=1e-06)
  
  
  while(error=="False"){
    t=step
    #print(t)
    #Loop over shocks space
   
    for(i in 1:length(xgrid)){
      
      b=xgrid[i]
      guess=0.1
      
      # Initial guess
      if (t == 1) {
          guess <- 0.1
        
      } else {
        # Use previous iteration's solution as warm start

          guess <- polnow[i]
      }
      
      
      payout=p*b*gamma
      premium=payout*payout
      #output = nloptr(x0=guess,eval_f=Payoff_pre,lb=0.0001,ub=1,opts=options,b=b,V=V,sim=sim,r=r,K=K,delta=delta,ra=a,ut_mod=ut_mod,xgrid=xgrid,p=p,c=c,cshape=cshape,trigger=trigger,premium=premium,gamma=gamma)
      
      output = optim(par=guess,fn=Payoff_pre,lower=0.00001,upper=1,b=b,V=V,sim=sim,r=r,K=K,delta=delta,ra=a,ut_mod=ut_mod,xgrid=xgrid,p=p,c=c,cshape=cshape,trigger=trigger,premium=premium,gamma=payout,method='L-BFGS-B')
      fstar = output$par
      Vstar = -output$val
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
  
  conv<-DFall |> 
    filter(time==max(DFall$time)) |> 
    select(-time)
  
  parameters<-data.frame(
    r=r,
    K=K,
    small=small,
    p=p,
    c=c,
    a=a,
    delta=delta,
    cshape=cshape,
    mean_theta=mean_theta,
    sigma_theta=sigma_theta,
    gamma=gamma,
    trigger=trigger,
    premium=premium,
    trigger=trigger,
    ut_mod=ut_mod
  )
  
  policy<-list(conv=conv,parameters=parameters)
 
  return(policy)
}
