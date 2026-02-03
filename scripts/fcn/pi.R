pi = function(f,b,theta,p,c,cshape,trigger,premium,gamma)
{
 
  #Make Calculations with insurance, insurance defined globally in param space

    if(is.na(trigger)){
      payout=0
      premium=0
    } else if(theta<=trigger){ #check if we payout, note: constant payment for now so premium and payout easy to calculate
      payout=gamma*(trigger-theta)
    }else{
      payout=0
    }
    
    
    profit=(p*harvest(f=f,biomass=b,theta = theta))-(c*(f^cshape))+payout-premium
  
  return(profit)
}
