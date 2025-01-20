# Definitions of weighted average and normalization functions

  # Weighted contribution of datas function and transformation into nimbleFunction
  # Weighted contribution of datas function : f <- function(x, mean, sd) {x * dnorm(x, mean, sd)}
nf <- nimbleFunction(
  run = function(x = double(1), omega = double(1)) {
    
    numerat <- x * dnorm(x, omega[1], omega[2]) # Calculate the numerator for weighted average (mean and sd)
    
    returnType(double(1))
    return(numerat)
  }
)

  # Normalization function and transformation into nimbleFunction
  # Normalization function : f2 <- function(x, mean, sd) {dnorm(x, mean, sd)}
nf2 <- nimbleFunction(
  run = function(x = double(1), omega = double(1)) {
    
    denomin <- dnorm(x, omega[1], omega[2]) # Calculate the denominator for normalization (mean and sd)
    
    returnType(double(1))
    return(denomin)
  }
)

# Midpoint interval function
nf_midbins <- nimbleFunction(
  run = function(mu = double(0), sd = double(0), nb_LC = double(0), Lmin = double(0), Lmax = double(0)) {
    # Definition of omega
    omega <- nimC(mu, sd)
    
    # Definition of min/max bounds
    borne_min <- 0
    borne_max <- 5
    
    # Calculate breakpoints based on the value range
    breaks_base <- seq(Lmin, Lmax, length.out = nb_LC - 1)
    breaks_names <- c(borne_min, breaks_base, borne_max)
    
    # Initialization of vectors
    int_num <- numeric(length = nb_LC)
    int_denom <- numeric(length = nb_LC)
    
    # Loop to calculate integrals over each interval
    for (cl in 1:nb_LC) {
      # Calculate the weighted average for the interval
      int_num[cl] <- nimIntegrate(nf, lower = breaks_names[cl], upper = breaks_names[cl + 1], param = omega)[1] # Index 1 to retrieve only the value
      
      # Calculate the normalization for the interval
      int_denom[cl] <- nimIntegrate(nf2, lower = breaks_names[cl], upper = breaks_names[cl + 1], param = omega)[1] # Index 1 to retrieve only the value
    }
    
    # Finally, calculate mid_breaks as the ratio of weighted average and normalization
    midbins <- int_num / int_denom
    
    returnType(double(1))
    return(midbins)
  }
)
