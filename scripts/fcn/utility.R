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
