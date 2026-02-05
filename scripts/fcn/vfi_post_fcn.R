vfi_post_fcn<-function(sizex,n_z,cshape,delta,
                       r,K,p,c,a,
                       mean_theta,sigma_theta,
                       gamma,trigger,m,
                       ut_mod){

  
  
  small=K/100

  #### Set up state grid ####
  xgrid = seq(small,K,length.out=sizex)
  

  set.seed(42)
  
  
  shock_process <- discretize_normal_shock(n_z, mean_theta, sigma_theta)
  w_grid <- shock_process$w_grid
  w_probs <- shock_process$probs
  
  premium_rate<-sum(map_dbl(trigger-w_grid,~max(0,.x))*w_probs)*(1+m) 
  
  #### Set up maximization settings and storage ####
  options=list("algorithm"="NLOPT_LN_COBYLA",xtol_rel=1e-06)
  DFall = data.frame()
  Vnext = matrix(ncol = length(w_probs),nrow =sizex)
  V = matrix(0,ncol=length(w_probs),nrow=sizex)
  pol_converged=matrix(0,nrow=sizex,ncol = length(w_probs)) #Harvest Policy function
  polnow=matrix(0,nrow=sizex,ncol = length(w_probs))
  
  # Define errors
  step=1
  error="False"
  tol=0.001
  
  # Run backwards induction
  while(error=="False"){
    t=step
    print(t)

    #Loop over shocks space
    for(j in 1:length(w_probs)){
      
      temptheta=w_grid[j]
      
      
      for(i in 1:length(xgrid)){
        
        #print(paste0('State: ',i,' of ',length(xgrid),', Shock: ',j,' of ',length(theta_prop),', Time: ',t))
        b=xgrid[i]
        # Initial guess
        if (t == 1) {
          guess <- 0.1
          
        } else {
          # Use previous iteration's solution as warm start
          
          guess <- polnow[i,j]
        }
        
        
        payout=p*b*gamma
        premium=payout*premium_rate
      
        output = nloptr(x0=guess,eval_f=Payoff_post,lb=0.0001,ub=1,opts=options,b=b,V=V,theta=temptheta,r=r,K=K,delta=delta,theta_prop=w_probs,ra=a,ut_mod=ut_mod,xgrid=xgrid,p=p,c=c,cshape=cshape,trigger=trigger,premium=premium,gamma=payout)
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
    if(step==100){
      print('Maximum number of iterations (100) reached.')
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
