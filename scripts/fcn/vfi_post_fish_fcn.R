vfi_post_fish_fcn<-function(n_shocks,mu_w,sigma_w,m,trigger,n_biomass,K,r,c_cost,pay,beta,gamma,ut_mod){

  # Biomass grid
  biomass_min <- K/100       # Minimum viable biomass
  biomass_max <- K       # Carrying capacity
  biomass_grid <- seq(biomass_min, biomass_max, length.out = n_biomass)
  max_coverage <- K  
  
  shock_process <- discretize_normal_shock(n_shocks, mu_w, sigma_w)
  w_grid <- shock_process$w_grid
  w_probs <- shock_process$probs
  
  
  premium_rate<-sum(pmax(trigger - w_grid, 0) * w_probs) * (1 + m)
  DFall = data.frame()
  Vnext = matrix(ncol = length(w_probs),nrow =n_biomass)
  V = matrix(0,ncol=length(w_probs),nrow=n_biomass)
  pol_converged=matrix(0,nrow=n_biomass,ncol = length(w_probs)) #Harvest Policy function
  polnow=matrix(0,nrow=n_biomass,ncol = length(w_probs))
  

  
  step=1
  error="False"
  tol=0.01
  
  #options=list("algorithm"="NLOPT_LN_COBYLA",xtol_rel=1e-06)
  
  
  while(error=="False"){
    t=step
    
    #Loop over shocks space
    
    for(i_w in 1:n_shocks){
      shock=w_grid[i_w]
      
      for(i in 1:n_biomass){
        
        b=biomass_grid[i]
        guess=b/10
        
        output = optim(par=guess,
                       fn=payoff_post_fcn,
                       lower=0,
                       upper=b*(1+shock),
                       b=b,
                       shock=shock,
                       w_grid=w_grid,
                       V=V,
                       gamma=gamma,
                       mod=ut_mod,
                       pay=pay,
                       prem_rate=premium_rate,
                       c_cost=c_cost,
                       w_probs=w_probs,
                       r=r,
                       K=K,biomass_grid=biomass_grid,
                       trigger=trigger,
                       beta=beta,
                       method='Brent') 
        
        
        ins_obj<-optim(par=b/10,
                       fn=fisher_ut,
                       b=b,
                       q=output$par,
                       c_cost=c_cost,
                       w_probs=w_probs,
                       shock=w_grid,
                       mod=ut_mod,
                       gamma=gamma,
                       prem_rate=premium_rate,
                       trigger=trigger,
                       method='Brent',
                       lower=0,
                       upper=max_coverage)
        
        cstar<-ins_obj$par
        
        
        
        fstar = output$par
        Vstar = -output$val
        Vnext[i,i_w]=Vstar
        DFnow = data.frame(time=t,b=b,w=shock,fstar=fstar,cstar=cstar,Vstar=Vstar)
        DFall = bind_rows(DFall,DFnow)
        polnow[i,i_w]=fstar
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
    gamma=gamma,
    beta=beta,
    c_cost=c_cost,
    mu_w=mu_w,
    sigma_w=sigma_w,
    trigger=trigger,
    m=m,
    ut_mod=ut_mod,
    pay=pay,
    n_biomass=n_biomass,
    n_shocks=n_shocks
  )
  
  policy<-list(conv=conv,parameters=parameters)
  
  return(policy)
}
