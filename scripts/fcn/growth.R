growth = function(b,f,theta,r,K)
{
  
  escapement=b*(1+theta)-harvest(f=f,biomass=b,theta = theta)
  

  if(is.na(escapement)) {
    print(paste0("b: ",b," f: ",f," theta: ",theta))
  }
  
  if(escapement>K){
    bnext=K
  }else if(escapement<0){
    bnext=K*0.001
  }  else{
    
    bnext = escapement + r*escapement*(1-escapement/K)
  }
  
  return(bnext)
}
