# Créer la fonction Nimble
nf_prop <- nimbleFunction(
  run = function(mu = double(0), sd = double(0), nb_LC = double(0), Lmin = double(0), Lmax = double(0)) {
    # Definition of omega
    omega <- nimC(mu, sd)
    
    # Définition des bornes min/max
    borne_min <- 0
    borne_max <- 5
    
    # Calculer les points de rupture en fonction de la plage de valeurs
    breaks_base <- seq(Lmin, Lmax, length.out = (nb_LC - 1))
    breaks_names <- c(borne_min, breaks_base, borne_max)
    
    # Calculer les probabilités cumulatives pour chaque intervalle
    prop_intervals <- pnorm(breaks_names, mean = omega[1], sd = omega[2])
    
    # Calculer la différence entre les probabilités cumulatives consécutives
    prop <- prop_intervals[2:(nb_LC + 1)] - prop_intervals[1:(nb_LC)]
    
    returnType(double(1))
    return(prop)
  }
)