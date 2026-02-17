#### Weighted optimization VFI ####

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
gamma <- 0.15          # Risk aversion (CRRA utility)


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
cov_rate<-0
# Iteration parameters
max_iter <- 300
tol <- 1e-4
alpha=0 # Weight on harvest vs utility (1=harvest only, 0=utility only)


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
ut_max<-function(effort,biomass,quota,pay,premium){
 
  e_harvest<-map_dbl(.x=w_grid,
                  .f=~harvest_function(effort=effort,
                                       biomass=biomass,
                                       quota=quota,
                                       shock=.x))
  
  
  pi<-e_harvest-c_cost*effort^2+pay-premium
  
  e_ut<-map_dbl(.x=pi,
             .f=~utility(c=.x,gamma=gamma,boost=0,mod='cara'))*w_probs
  
  return(-sum(e_ut*w_probs))
}

objective_fcn<-function(choice,b,V_now,pay_vec,premium){

  quota<-choice[1]
  
  obj<-optim(par=0.1,
        fn=ut_max,
        biomass=b,
        quota=quota,
        pay=pay_vec,
        premium=premium,
        method="L-BFGS-B",
        lower=0,
        upper=effort_max)

  effort_choice<-obj$par

  e_harvest<-map_dbl(.x=w_grid,
                  .f=~harvest_function(effort=effort_choice,
                                       biomass=b,
                                       quota=quota,
                                       shock=.x))
  
  esc<-map2_dbl(.x=e_harvest,.y=w_grid,~max((1+.y)*b-.x,biomass_min))
  
  
  xnext<-map_dbl(esc,~growth(.x))
  
  temp<-map_dbl(xnext,~spline(x=biomass_grid,y=V_now,xout=.x,method='natural')$y)
  
  pay<-pay_vec
  
  pi<-e_harvest-c_cost*effort_choice^2+pay-premium
  
  e_ut<-map_dbl(.x=pi,
                .f=~utility(c=.x,gamma=gamma,boost=0,mod='cara'))*w_probs
  
  if(alpha==1){
    
    return(-(sum(e_harvest*w_probs)+beta*sum(temp*w_probs)))
  
    } else if(alpha==0){
    
    return(-((1-alpha)*sum(e_ut)+beta*sum(temp*w_probs)))
  
    } else{
    
    scale=max_ut/max_harvest # Make sure we have a more equal comparison between harvest and ut
    
    return(-(alpha*sum(e_harvest*w_probs)*scale+(1-alpha)*sum(e_ut)+beta*sum(temp*w_probs)))
  }
  
  
}


DFall = data.frame()
Vnext = vector(mode='numeric',length=n_biomass)
V = rep(0,length.out=n_biomass)
pol_converged=vector(mode='numeric',length=n_biomass) #Harvest Policy function
polnow=vector(mode='numeric',length=n_biomass)

step=1
error="False"
tol=0.01

#options=list("algorithm"="NLOPT_LN_COBYLA",xtol_rel=1e-06)


while(error=="False"){
  t=step
 
  #Loop over shocks space
  
  for(i in 1:n_biomass){

    b=biomass_grid[i]
    guess=b/10

    
    #output = nloptr(x0=guess,eval_f=Payoff_pre,lb=0.0001,ub=1,opts=options,b=b,V=V,sim=sim,r=r,K=K,delta=delta,ra=a,ut_mod=ut_mod,xgrid=xgrid,p=p,c=c,cshape=cshape,trigger=trigger,premium=premium,gamma=gamma)
    output = optim(par=guess,
                   fn=objective_fcn,
                   lower=0.00001,
                   upper=b,
                   b=b,
                   V_now=V,
                   pay_vec=pay,
                   premium=premium,
                   method='Brent') 
    #output = optim(par=guess,fn=objective_fcn,lower=0.00001,upper=K,b=b,V_now=V,pay_vec=pay,premium=premium,method='L-BFGS-B')
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
  
  print(error_vec)
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
