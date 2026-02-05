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