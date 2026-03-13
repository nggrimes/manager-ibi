# Insurance payout
insurance_payout <- function(shock, coverage,trigger) {
  if (shock < trigger) {
    gap <- (trigger - shock)
    payout <- coverage * gap
    return(payout)
  } else {
    return(0)
  }
}
