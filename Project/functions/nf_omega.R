nf_omega <- nimbleFunction(
  run = function(N = double(2), midbins = double(1)) {
    # Declare variables
    nRows <- dim(N)[1]  # Get the number of rows in N
    N_sum <- numeric(nRows)  # Initialize a vector for the summed values
    
    # Loop to compute the row-wise sum
    for (i in 1:nRows) {
      N_sum[i] <- sum(N[i, ])
    }
    
    total_N <- sum(N_sum)  # Total sum of N
    
    # Compute the weighted mean and standard deviation
    mu <- sum(midbins * N_sum) / total_N
    sd <- sqrt(sum(N_sum * (midbins - mu)^2) / total_N)
               
    omega <- nimC(mu, sd)
               
    returnType(double(1))
    return(omega)
  }
)