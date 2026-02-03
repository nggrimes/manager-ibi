vfi_post_fcn<-function(sizex,n_z,cshape,delta,
                       r,K,p,c,a,
                       mean_theta,sigma_theta,
                       gamma,coverage,
                       ut_mod){

  
  
  small=K/100

  #### Set up state grid ####
  xgrid = seq(small,K,length.out=sizex)
  
  #### Set up shock grid ####
  n_sims=10000
  
  set.seed(42)
  
  #Get distribution and bins of
  theta_raw<-rnorm(n_sims,mean=mean_theta,sd=sigma_theta)
  theta_dist<-theta_raw[which(theta_raw>-1)]
  
  thetahist<-hist(theta_dist,breaks=seq(min(theta_dist),max(theta_dist),l=n_z+1),plot=FALSE)
  ## Make bins classification
  bins=theta_dist
  for(k in 2:length(thetahist$breaks)){
    index=which(theta_dist>thetahist$breaks[k-1]&theta_dist<=thetahist$breaks[k])
    bins[index]=k-1
  }
  
  theta_df<-data.frame(theta_dist=theta_dist,bins=bins)
  
  means<- theta_df |> 
    group_by(bins) |> 
    summarize(mean=mean(theta_dist))
  
  #elminate first observation
  theta_prop=thetahist$counts/n_sims
  theta_val=means$mean
  
  # Use distribution to get insurance parameters
  trigger_index=findNearest(coverage,theta_val)
  trigger=theta_val[trigger_index]
  trigger_prop=sum(theta_prop[1:trigger_index])
  
  premium=gamma*trigger_prop
  
  #### Set up maximization settings and storage ####
  options=list("algorithm"="NLOPT_LN_COBYLA",xtol_rel=1e-06)
  DFall = data.frame()
  Vnext = matrix(ncol = length(theta_prop),nrow =sizex)
  V = matrix(0,ncol=length(theta_prop),nrow=sizex)
  pol_converged=matrix(0,nrow=sizex,ncol = length(theta_prop)) #Harvest Policy function
  polnow=matrix(0,nrow=sizex,ncol = length(theta_prop))
  
  # Define errors
  step=1
  error="False"
  tol=0.001
  
  # Run backwards induction
  while(error=="False"){
    t=step
    print(t)

    #Loop over shocks space
    for(j in 1:length(theta_prop)){
      
      temptheta=theta_val[j]
      
      
      for(i in 1:length(xgrid)){
        
        #print(paste0('State: ',i,' of ',length(xgrid),', Shock: ',j,' of ',length(theta_prop),', Time: ',t))
        b=xgrid[i]
        guess=0.1
        output = nloptr(x0=guess,eval_f=Payoff_post,lb=0.0001,ub=1,opts=options,b=b,V=V,theta=temptheta,r=r,K=K,delta=delta,theta_prop=theta_prop,ra=a,ut_mod=ut_mod,xgrid=xgrid,p=p,c=c,cshape=cshape,trigger=trigger,premium=premium,gamma=gamma)
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
    coverage=coverage,
    premium=premium,
    trigger=trigger,
    ut_mod=ut_mod
  )
  
  policy<-list(conv=conv,parameters=parameters)
  
return(policy)
}
