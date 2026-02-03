ut_post<-function(pi_out,ra,ut_mod='cara'){

  pi_out<-pmax(pi_out,1e-10)
  
  if(ut_mod=='cara'){
    out=1-exp(-ra*(pi_out))
  }else if(ut_mod=='log'){
    out=log(pi_out)  
  }else if (ut_mod=='rn'){
    out=pi_out
  } else if(ut_mod=='power'){
    out=pi_out^(1-ra)/(1-ra)  
  } else{
    print('Not a valid utility function. Please use cara: for constant absolute risk aversion, log: for log utility, or rn: for risk neutral')
  }
  
  return(out)
}
