growth <- function(escapement,r_growth,K) {
  biomass_min=K/1000
  growth <- escapement + r_growth * escapement * (1 - escapement / K)
  
  growth<-min(growth,K)
  return(max(growth, biomass_min))
}
