




set.seed(123)

# Economic parameters
beta <- 0.95           # Discount factor
gamma <- 0.5          # Risk aversion (CRRA utility)


# Harvest and profit parameters

p <- 1                # Price per unit harvest
c_cost <- 15          # Cost coefficient for effort (c*E^2)


# Biological parameters (logistic growth)
r_growth <- 0.8        # Intrinsic growth rate
K <- 100               # Carrying capacity

# Shock parameters (multiplicative on biomass)
mu_w <- 0              # Mean shock (centered at zero)
sigma_w <- 0.2         # Shock volatility

# Insurance parameters
trigger <- -0.01       # Shock trigger (payout when w < trigger)
m<-0               # Premium as % of coverage (actuarially unfair)

# State space - keep moderate size for speed
n_biomass <- 25        # Grid points for biomass
n_shocks <- 13         # Number of shock realizations (use odd number for symmetry)

# Biomass grid
biomass_min <- K/100       # Minimum viable biomass
biomass_max <- K       # Carrying capacity
biomass_grid <- seq(biomass_min, biomass_max, length.out = n_biomass)

# Choice bounds
effort_max <- 5       # Maximum effort
max_coverage <- K    # Maximum insurance coverage

# Iteration parameters
max_iter <- 300
tol <- 1e-4

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

pay<-0 # If we want to endogenize insurance choice into the quota decision or not

# ==============================================================================
# 4. BIOLOGICAL AND ECONOMIC FUNCTIONS
# ==============================================================================

# Logistic growth
growth <- function(escapement) {
  growth <- escapement + r_growth * escapement * (1 - escapement / K)
  
  growth<-min(growth,K)
  return(max(growth, biomass_min))
}



## Insurance parameters ##



fisher_ut<-function(cov,b,q,prem_rate,shock,gamma,mod='cara'){

  payout<-map_dbl(.x=w_grid,
                       .f=~insurance_payout(shock=.x,
                                            coverage=cov))
  premium<-cov*prem_rate
  
  pi<-map2_dbl(.x=shock,.y=payout,~q-c_cost*q/((1+.x)*b)+.y-premium)
  
  e_ut<-map_dbl(pi,~utility(.x,gamma,boost=0,mod=mod))*w_probs
  
  return(-sum(e_ut))
}

payoff_fcn<-function(q,b,shock,V,pay,prem_rate,gamma,mod){


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
                    mod=mod,
                    method='Brent',
                    lower=0,
                    upper=max_coverage)
    
    e_ut=-fish_obj$value
  }
  
  
  
  esc<-map2_dbl(.x=q,.y=w_grid,~max((1+.y)*b-.x,biomass_min))
  
  
  xnext<-map_dbl(esc,~growth(.x))
  
  temp<-map_dbl(xnext,~spline(x=biomass_grid,y=V,xout=.x,method='natural')$y)
  
  Vfuture=sum(temp*w_probs)
  
  return(-(e_ut+beta*Vfuture))
  
}

DFall = data.frame()
Vnext = vector(mode='numeric',length=n_biomass)
V = rep(0,length.out=n_biomass)
pol_converged=vector(mode='numeric',length=n_biomass) #Harvest Policy function
polnow=vector(mode='numeric',length=n_biomass)

step=1
error="False"
tol=0.001

#options=list("algorithm"="NLOPT_LN_COBYLA",xtol_rel=1e-06)


while(error=="False"){
  t=step
  
  #Loop over shocks space
  
  for(i in 1:n_biomass){
    
    b=biomass_grid[i]
    guess=b/10
    
    output = optim(par=guess,
                   fn=payoff_fcn,
                   lower=0,
                   upper=b,
                   b=b,
                   shock=w_grid,
                   V=V,
                   gamma=gamma,
                   mod=ut_mod,
                   pay=pay,
                   prem_rate=premium_rate,
                   method='Brent') 
    

    ins_obj<-optim(par=b/10,
                   fn=fisher_ut,
                                 b=b,
                                 q=output$par,
                                 prem_rate=premium_rate,
                                 shock=w_grid,
                                 gamma=gamma,
                                 mod=ut_mod,
                                 method='Brent',
                                 lower=0,
                                 upper=max_coverage)
    
    cstar<-ins_obj$par
    


    fstar = output$par
    Vstar = -output$val
    Vnext[i]=Vstar
    DFnow = data.frame(time=t,b=b,fstar=fstar,Vstar=Vstar,cstar=cstar)
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


## find fisher utility at each

noi_ut<-map2_dbl(.x=conv$fstar,.y=conv$b,~fisher_ut(0,.y,.x,0,w_grid))

i_ut<-map2(.x=conv$fstar,.y=conv$b,~optim(.y/10,fn=fisher_ut,b=.y,q=.x,prem_rate=premium_rate,shock=w_grid,method='Brent',lower=0,upper=max_coverage))
see<-bind_rows(i_ut)

