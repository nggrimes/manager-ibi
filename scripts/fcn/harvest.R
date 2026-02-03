harvest <- function(f, biomass, theta){
  
  h = f * biomass * (1+theta)
  
  return(h)
}
