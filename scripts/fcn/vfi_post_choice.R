

vfi_post_choice<-function(sizex,n_z,mean_theta,sigma_theta,p,c,cshape,r,K,ra,delta,m,trigger){
  #### Parameters ####
  set.seed(123)
  
  # Economic parameters
  beta <- delta          # Discount factor
  gamma <- ra          # Risk aversion (CRRA utility)
  
  # Harvest/profit parameters
  p <- p              # Price per unit harvest
  c_cost <- c        # Cost coefficient for harvest effort (c*f^2)
  c_shape<-cshape
  # Biological parameters (logistic growth)
  r <- r             # Intrinsic growth rate
  K <- K               # Carrying capacity
  
  # Shock parameters
  mu_w <- mean_theta           # Mean shock (centered at zero)
  sigma_w <- sigma_theta        # Shock volatility (proportional)
  
  # Insurance parameters
  m <- m   # Premium as % of coverage (actuarially unfair)
  trigger <- trigger       # Shock trigger (payout when w < trigger)
  
  
  # State space
  n_biomass <- sizex        # Grid points for biomass
  n_shocks <- n_z         # Number of shock realizations
  biomass_min <- K*0.001       # Minimum biomass on grid
  biomass_max <- K       # Maximum biomass (carrying capacity)
  biomass_grid <- seq(biomass_min, biomass_max, length.out = n_biomass)
  
  # Optimization bounds
  effort_max <- .999  # Maximum harvest effort
  max_coverage <- p * K * 0.5  # Maximum insurance coverage
  
  # Iteration parameters
  max_iter <- 500
  tol <- 1e-5
  
  #### Functions #####
  
  utility <- function(profit, gamma) {
    if (gamma == 1) {
      return(log(pmax(profit+100, 1e-10)))  # Log utility
    } else {
      profit_safe <- pmax(profit+100, 1e-10)
      return((profit_safe^(1 - gamma) - 1) / (1 - gamma))
    }
  }
  
  # Discretize i.i.d. normal shock process
  discretize_normal_shock <- function(n, mu, sigma, m = 3) {
    w_max <- mu + m * sigma
    w_min <- mu - m * sigma
    
    w_grid <- seq(w_min, w_max, length.out = n)
    step <- (w_max - w_min) / (n - 1)
    
    # Probabilities (i.i.d., so same for all current states)
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
  
  logistic_growth <- function(escapement) {
    growth <- escapement + r * escapement * (1 - escapement / K)
    return(pmax(growth, 0))
  }
  
  insurance_payout <- function(shock, coverage) {
    if (shock < trigger) {
      gap <- (trigger - shock)
      payout <- coverage *gap
      return(payout)
    } else {
      return(0)
    }
  }
  
  
  ### Calculate shocks and premium ####
  shock_process <- discretize_normal_shock(n_shocks, mu_w, sigma_w)
  w_grid <- shock_process$w_grid
  w_probs <- shock_process$probs
  
  premium_rate<-sum(map_dbl(trigger-w_grid,~max(0,.x))*w_probs)*(1+m) # add the insurance markup
  
  solve_vfi <- function(with_insurance = TRUE) {
    
    
    # Initialize value function and policy functions
    V <- matrix(0, n_biomass, n_shocks)
    V_new <- matrix(0, n_biomass, n_shocks)
    
    policy_effort <- matrix(0, n_biomass, n_shocks)

    
    if (with_insurance) {
      policy_coverage <- matrix(0, n_biomass, n_shocks)
      policy_premium <- matrix(0, n_biomass, n_shocks)
    }
    
    
    # Objective function to MINIMIZE (negative of value function)
    objective_fn <- function(choices, current_biomass, w_shock, V_current) {
      
      if (with_insurance) {
        effort <- choices[1]
        coverage <- choices[2]
      } else {
        effort <- choices[1]
        coverage <- 0
      }
      
      # Bounds checking (optim should handle this, but just in case)
      if (effort < 0 || effort > effort_max) return(1e10)
      if (coverage < 0 || coverage > max_coverage) return(1e10)
      
      # Premium payment
      premium <- premium_rate * coverage
      
      # Available biomass and harvest
      available_biomass <- pmax(current_biomass * (1 + w_shock), 1e-6)
      harvest <- effort * available_biomass
      
      # Escapement
      escapement <- pmax(available_biomass - harvest, 1e-6)
      
      # Profit this period
      profit_base <- p * harvest - c_cost * effort^c_shape
      payout <- insurance_payout(w_shock, coverage)
      profit_total <- profit_base + payout - premium
      
      # Next period biomass
      biomass_next <- logistic_growth(escapement)
  
      sp_store<-vector(mode="numeric",length=ncol(V_current))
      
      #Find the expected future value via spline interpolation for each shock state
      for(n in 1:ncol(V_current)){
        
        temp<-spline(x=biomass_grid,y=V_current[,n],xout=biomass_next,method="natural")
        sp_store[n]<-temp$y
        
        
      }
      
      Vfuture=delta*sp_store*w_probs
      
      # Total value
      total_value <- utility(profit_total,gamma) + sum(Vfuture)
      
      # Return negative for minimization
      return(-total_value)
    }
    
    # Iteration
    for (iter in 1:max_iter) {
      for (i_w in 1:n_shocks){
        w_shock <- w_grid[i_w]
        # Loop over current biomass state
        for (i_b in 1:n_biomass) {
          
          current_biomass <- biomass_grid[i_b]
          
          # Initial guess
          if (iter == 1) {
            if (with_insurance) {
              init_guess <- c(0.3, max_coverage * 0.3)
            } else {
              init_guess <- c(0.3)
            }
          } else {
            # Use previous iteration's solution as warm start
            if (with_insurance) {
              init_guess <- c(policy_effort[i_b, i_w], policy_coverage[i_b, i_w])
            } else {
              init_guess <- c(policy_effort[i_b, i_w])
            }
          }
          
          # Optimize using L-BFGS-B (handles box constraints efficiently)
          if (with_insurance) {
            opt_result <- optim(
              par = init_guess,
              fn = objective_fn,
              current_biomass = current_biomass,
              V_current = V,
              w_shock=w_shock,
              method = "L-BFGS-B",
              lower = c(0, 0),
              upper = c(effort_max, max_coverage)
            )
            
            best_effort <- opt_result$par[1]
            best_coverage <- opt_result$par[2]
          } else {
            opt_result <- optim(
              par = init_guess,
              fn = objective_fn,
              current_biomass = current_biomass,
              w_shock=w_shock,
              V_current = V,
              method = "L-BFGS-B",
              lower = c(0),
              upper = c(effort_max)
            )
            
            best_effort <- opt_result$par[1]
            best_coverage <- 0
          }
          
          max_value <- -opt_result$value  # Convert back to positive
          
          # Calculate metrics for storage
          premium <- premium_rate * best_coverage
          
          # Store results
          V_new[i_b, i_w] <- max_value
          policy_effort[i_b, i_w] <- best_effort
          
          if (with_insurance) {
            policy_coverage[i_b, i_w] <- best_coverage
            policy_premium[i_b, i_w] <- premium
          }
        }
      }
      # Check convergence
      diff <- colMeans(abs(V_new - V))
    
      if (iter %% 50 == 0) {
        cat(sprintf("Iteration %d: max diff = %.6f\n", iter, sum(diff)))
      }
      
      if (sum(diff<tol)==length(diff)) {
        cat(sprintf("Converged in %d iterations!\n", iter))
        break
      }
      
      
      V <- V_new
    }
    
    # Return results
    result <- list(
      V = V,
      fstar = policy_effort,
      b = biomass_grid
    )
    
    if (with_insurance) {
      result$coverage = policy_coverage
      result$premium = policy_premium
    }
    
    return(result)
  }
  

  
  solution_no_ins <- solve_vfi(with_insurance = FALSE)
  
  solution_with_ins <- solve_vfi(with_insurance = TRUE)
  
  parameters<-data.frame(
    r=r,
    K=K,
    p=p,
    c=c,
    a=ra,
    delta=delta,
    cshape=cshape,
    mean_theta=mean_theta,
    sigma_theta=sigma_theta,
    trigger=trigger,
    ut_mod=ut_mod,
    m=m,
    sizex=sizex,
    n_z=n_z,
    trigger=trigger
  )
  
  policy<-list(ins=solution_with_ins,noi=solution_no_ins,parameters=parameters)
  
  return(policy)
  
}
