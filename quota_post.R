#### Post-Shock With Quota ####

library(ggplot2)
library(reshape2)
library(gridExtra)
library(fields)
library(tidyverse)

# ==============================================================================
# 1. MODEL PARAMETERS
# ==============================================================================

set.seed(123)

# Economic parameters
beta <- 0.95           # Discount factor
gamma <- 0.1          # Risk aversion (CRRA utility)


# Harvest and profit parameters
q <- 0.2               # Catchability coefficient
p <- 1                # Price per unit harvest
c_cost <- 2.0          # Cost coefficient for effort (c*E^2)


# Biological parameters (logistic growth)
r_growth <- 0.8        # Intrinsic growth rate
K <- 100               # Carrying capacity

# Shock parameters (multiplicative on biomass)
mu_w <- 0              # Mean shock (centered at zero)
sigma_w <- 0.2         # Shock volatility

# Insurance parameters
trigger <- -0.11       # Shock trigger (payout when w < trigger)
m<-0.2               # Premium as % of coverage (actuarially unfair)

# State space - keep moderate size for speed
n_biomass <- 25        # Grid points for biomass
n_shocks <- 11         # Number of shock realizations (use odd number for symmetry)

# Biomass grid
biomass_min <- K/100       # Minimum viable biomass
biomass_max <- K       # Carrying capacity
biomass_grid <- seq(biomass_min, biomass_max, length.out = n_biomass)

# Choice bounds
effort_max <- 5       # Maximum effort
max_coverage <- q*effort_max*K    # Maximum insurance coverage
cov_rate<-0.5
# Iteration parameters
max_iter <- 300
tol <- 1e-4
alpha=1 # Weight on harvest vs utility (1=harvest only, 0=utility only)


# ==============================================================================
# 2. UTILITY FUNCTION (CRRA)
# ==============================================================================

utility <- function(c, gamma,boost=100,mod='power') {
  c <- c+boost  # Ensure positive
  
  if(mod=='cara'){
    ut<-(1-exp(-gamma*c))/gamma
    return(ut)
  }
  else if (gamma == 1) {
    return(log(c))
  } else {
    return((c^(1 - gamma) - 1) / (1 - gamma))
  }
}

# ==============================================================================
# 3. SHOCK PROCESS DISCRETIZATION
# ==============================================================================

discretize_normal_shock <- function(n, mu, sigma, m = 3) {
  w_max <- mu + m * sigma
  w_min <- mu - m * sigma
  
  w_grid <- seq(w_min, w_max, length.out = n)
  step <- (w_max - w_min) / (n - 1)
  
  probs <- numeric(n)
  
  for (j in 1:n) {
    if (j == 1) {
      probs[j] <- pnorm((w_grid[1] - mu + step/2) / sigma)
    } else if (j == n) {
      probs[j] <- 1 - pnorm((w_grid[n] - mu - step/2) / sigma)
    } else {
      upper <- (w_grid[j] - mu + step/2) / sigma
      lower <- (w_grid[j] - mu - step/2) / sigma
      probs[j] <- pnorm(upper) - pnorm(lower)
    }
  }
  
  probs <- probs / sum(probs)
  
  return(list(w_grid = w_grid, probs = probs))
}

shock_process <- discretize_normal_shock(n_shocks, mu_w, sigma_w)
w_grid <- shock_process$w_grid
w_probs <- shock_process$probs


premium_rate<-sum(pmax(trigger - w_grid, 0) * w_probs) * (1 + m)

# ==============================================================================
# 4. BIOLOGICAL AND ECONOMIC FUNCTIONS
# ==============================================================================

# Logistic growth
growth <- function(escapement) {
  growth <- escapement + r_growth * escapement * (1 - escapement / K)
  
  growth<-min(growth,K)
  return(max(growth, biomass_min))
}

# Harvest function: h = q * E * B * (1 + w)
harvest_function <- function(effort, biomass, shock,quota) {
  
  available_biomass <- biomass * (1 + shock)
  available_biomass <- max(available_biomass, 0)
  harvest <- q * effort * available_biomass
  harvest<-min(harvest,quota)
  return(harvest)
}

# Profit from harvest (before consumption and insurance)
profit <- function(harvest, effort) {
  revenue <- p * harvest
  cost <- c_cost * effort^2
  return(revenue - cost)
}

# Insurance payout
insurance_payout <- function(shock, coverage) {
  if (shock < trigger) {
    gap <- (trigger - shock)
    payout <- coverage * gap
    return(payout)
  } else {
    return(0)
  }
}

## Insurance parameters ##
pay<-map_dbl(.x=w_grid,
             .f=~insurance_payout(shock=.x,
                                  coverage=cov_rate*max_coverage))

premium<-premium_rate*cov_rate*max_coverage

# Max functions
ut_max_post<-function(effort,biomass,w,quota,pay,premium){
  
  e_harvest<-harvest_function(effort=effort,
                              biomass=biomass,
                              shock=w,
                              quota=quota)
  
  
  pi<-e_harvest-c_cost*effort^2+pay-premium
  
  e_ut<-map_dbl(.x=pi,
                .f=~utility(c=.x,gamma=gamma,boost=0,mod='cara'))
  
  return(-sum(e_ut))
}

objective_fcn_post<-function(choice,b,w,V_now,pay_vec,premium){
  
  quota<-choice[1]
  
  obj<-optim(par=0.1,
             fn=ut_max_post,
             biomass=b,
             quota=quota,
             w=w,
             pay=pay_vec,
             premium=premium,
             method="L-BFGS-B",
             lower=0,
             upper=effort_max)
  
  effort_choice<-obj$par
  
  e_harvest<-harvest_function(effort=effort_choice,
                             biomass=b,
                             shock=w,
                             quota=quota)
  
  esc<-min(b*(1+w)-e_harvest,biomass_min)
    
  
  xnext<-growth(esc)
  
  sp_store<-vector(mode="numeric",length=ncol(V_now))
  
  for(n in 1:ncol(V_now)){
    
    temp<-spline(x=biomass_grid,y=V_now[,n],xout=xnext,method="natural")
    sp_store[n]<-temp$y
    
    
  }
 
  
  Vfuture=beta*sum(sp_store*w_probs)
  
  
  pay<-pay_vec
  
  pi<-e_harvest-c_cost*effort_choice^2+pay-premium
  
  e_ut<-utility(c=pi,gamma=gamma,boost=0,mod='cara')
  
  if(alpha==1){
    
    return(-(sum(e_harvest)+Vfuture))
    
  } else if(alpha==0){
    
    return(-((1-alpha)*sum(e_ut)+Vfuture))
    
  } else{
    
    scale=max_ut/max_harvest # Make sure we have a more equal comparison between harvest and ut
    
    return(-(alpha*sum(e_harvest)*scale+(1-alpha)*sum(e_ut)+Vfuture))
  }
  
  
}


DFall = data.frame()
# Initialize value function and policy functions
V <- matrix(0, n_biomass, n_shocks)
Vnext <- matrix(0, n_biomass, n_shocks)

 #Harvest Policy function
polnow=matrix(0, n_biomass, n_shocks)

step=1
error="False"
tol=0.01

#options=list("algorithm"="NLOPT_LN_COBYLA",xtol_rel=1e-06)


while(error=="False"){
  t=step
  
  #Loop over shocks space
  for(i_w in 1:n_shocks){
    w=w_grid[i_w]
    for(i in 1:n_biomass){
      
      b=biomass_grid[i]
      guess=b/10
      output = optim(par=guess,
                     fn=objective_fcn_post,
                     lower=0.00001,
                     upper=b*(1+w),
                     b=b,
                     V_now=V,
                     w=w,
                     pay_vec=pay[i_w],
                     premium=premium,
                     method='Brent') 
      fstar = output$par
      Vstar = -output$val
      Vnext[i,i_w]=Vstar
      DFnow = data.frame(time=t,b=b,w=w,fstar=fstar,Vstar=Vstar)
      DFall = bind_rows(DFall,DFnow)
      polnow[i,i_w]=fstar
    }
  
  }
  
  comp<-abs(polnow-pol_converged)
  error_vec<-colMeans(comp)
  print(step)
  #Check if errors reach tolerance for policy function convergence
  if(sum(error_vec<tol)==length(error_vec) & step>10){
    error="True"
  }else{
    error="False"
  }
  
  step=step+1
  #browser()
  V=Vnext
  pol_converged=polnow
  
  
  #Check to make sure we stop eventually
  if(step==200){
    print('Maximum number of iterations (400) reached.')
    break
  }
}

conv<-DFall |> 
  filter(time==max(DFall$time)) |> 
  select(-time)

parameters<-data.frame(
  r=r_growth,
  K=K,
  small=biomass_min,
  p=p,
  a=ra,
  beta=beta,
  cshape=2,
  mu_w=mu_w,
  sigma_w=sigma_w,
  gamma=gamma,
  premium=premium,
  trigger=trigger,
  cov_rate=cov_rate,
  m=m,
  q=q,
  alpha=alpha
)

policy<-list(conv=conv,parameters=parameters)

save(policy,file=here::here("data","vfi_pre_quota_noi_alpha1.Rdata"))
