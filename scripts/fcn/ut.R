ut<-function(f,b,theta,trigger,pi,payout,premium,ut_mod='cara',ra){
  
  
  if(theta<trigger){
    pi_out<-pi+payout-premium
  } else{
    pi_out<-pi-premium
  }
  
  
  if(ut_mod=='cara'){
    out=1-exp(-ra*(pi_out))
  }else if(ut_mod=='log'){
    out=log(pi_out+1000)  #Plus 1000 is to make sure we dont have negative values
  }else if (ut_mod=='rn'){
    out=pi_out
  } else if(ut_mod=='power'){
    out=((pi_out+1000)^(1-ra))/(1-ra)  #Plus 1000 is to make sure we dont have negative values
  } else{
    print('Not a valid utility function. Please use cara: for constant absolute risk aversion, log: for log utility, or rn: for risk neutral')
  }
  
  return(out)
}
