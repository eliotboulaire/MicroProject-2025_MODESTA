model_code <- nimbleCode({
  # -------------------------------------------------------------------------
  ##             ATLANTIC SALMON SIZE-STRUCTURED LIFE CYLE MODEL
  ##                        FOCUS ON THE MARINE PHASE
  ##                           (SALMSIZE-MAR MODEL)
  ## 
  ## @ Eliot BOULAIRE
  ## Supervised by Etienne RIVOT & Marie NEVOUX
  ## Version 17/12/2024
  # -------------------------------------------------------------------------
  
  
  # -------------------------------------------------------------------------
  #                  GLOSSARY OF PRINCIPAL VARIABLE NAMES
  # -------------------------------------------------------------------------
  
  # Life stages
  # ----------------
  ## 3 | Smolts
  ## 4 | Post-smolts
  ### 41 | Post-smolts with size structure of migrating smolt 
  ### 42 | Post-smolts with size structure of post-smolt at the end of the first summer
  ## 5 | Maturing post-smolts
  ## 8 | Non-maturing post-smolts
  ## 6 | 1SW adult returns
  ## 9 | 2SW adult returns
  
  # Indices
  # ----------------
  ## "(stage)" can be one of the life stages aforementioned
  ##
  ## L(stage)[t] = 22>l>60 | Number of individual scale in the population sample (by cohorts)
  ## T = 24 | Number of cohorts of smolts migrating downstream (1996-2019)
  ## S = 2 | Number of sexes (females and males)
  ## LC = 50 | Number of population length classes created
  
  # Parameters
  # ----------------
  ## "(stage)" can be one of the life stages aforementioned
  ##
  ## Abundances
  ### N(stage)[lc, t, s] | Number of individuals (by length classes, cohorts & sex)
  ### Tot_N(stage)[t, s] | Total number of individuals (by cohorts & sex)
  ### Mu_N(stage)[t] | Mean number of individuals connected to the lognormal likelihood (by cohorts)
  ##
  ## Sex ratio
  ### Psex_SR(stage)[t, s] | Sex proportion in the population (by cohorts & sex)
  ### Pfem_SR(stage)[t] | Female proportion in the population connected to the binomial likelihood (by cohorts)
  ##
  ## Scale size structure
  ### Mu_L(stage)[t] | Mean scale size of the population size structure connected to the Gaussian likelihood (by cohorts)
  ### Sd_L(stage)[t] | Standard deviation of the population size structure connected to the Gaussian likelihood (by cohorts)
  ##
  ## Scale size histogram
  ### Midbins_LC(stage)[lc, t] | Mean length of each length classes (by length classes & cohorts)
  ### Prop_LC(stage)[lc, t] | Proportion of individuals in each length classes (by length classes & cohorts)
  
  # Demographic transitions
  # ----------------
  ## Survival rate (first months at sea)
  ### Theta3[lc, t] | Survival rate in natural scale (by length class & cohorts)
  ### Beta3 and Alpha3 are the parameters of the logistic regression of survival we want to estimate
  ### Min_Theta3 = 0 / Max_Theta3 = 1 | minimum and maximum survival authorized by logistic regression
  ##
  ## Maturation rate (end of the first summer at sea)
  ### Theta4[lc, t] | Maturation rate in natural scale (by length class, cohorts)
  ### Beta4 and Alpha4 are the parameters of the logistic regression of maturation we want to estimate
  ##
  ## Fixed survival rate between maturation decision and adult returns
  ### Theta5 = exp(-(M * Delta5)) | Fixed survival rate between maturation and return as 1SW
  ### M = 0.03 | Sea mortality rate by months
  ### Delta5 = 9 | Number of months (between november y and july y+1)
  ###
  ### Theta6 = exp(-(M * Delta6))| Fixed survival rate between maturation and return as 2SW
  ### M = 0.09 | Sea mortality rate by months
  ### Delta6 = 17 | Number of months (between november y and may y+2)
  
  # Likelihood
  # ----------------
  ## "(stage)" can be one of the life stages aforementioned
  ## 
  ## Eff_N(stage) | Number of individuals (by cohorts) estimated by a CMR model
  ##
  ## Scales_L(stage) | Each scale size (by individuals & cohorts)
  ##
  ## Nfem_SR(stage) | Number of females (by cohorts)
  
  # Constants
  # ----------------
  ## Min and Max
  ### Min_L3 = 0.4 / Max_L3 = 1.9 | minimum and maximum scale size of stage 3 for all cohorts
  ### Min_L42 = 1.7 / Max_L42 = 3.9 | minimum and maximum scale size of stage 42 for all cohorts
  ### Min_L6 = 3.1 / Max_L6 = 5.4 | minimum and maximum scale size of stage 6 for all cohorts
  ### Min_L9 = 3.6 / Max_L9 = 6.2 | minimum and maximum scale size of stage 9 for all cohorts
  ##
  ## Mean
  ### Mean_L3 = 1.053 | mean scale size of stage 3 for all cohorts
  ### Mean_L42 = 2.657 | mean scale size of stage 42 for all cohorts
  
  
  # -------------------------------------------------------------------------
  ##                            Stage 3 -> Stage 4
  # -------------------------------------------------------------------------
  
  # --------------------------------------------------------
  ##                          PRIOR
  # --------------------------------------------------------
  
  # Parameters of Stage 3
  # -----------------------------
  ## Abundance parameter
  for (t in 1:C) {
    Log_Mu_N3[t] <- log(Mu_N3[t]) - ((1/2) * Log_Sd_N3[t])
    Mu_N3[t] ~ dunif(min = 0, max = 50000)
  }
  
  ## size structure parameters
  for (t in 1:C) {
    ### Scale size structure normal parameters
    Mu_L3[t] ~ dunif(min = 0, max = 10)
    Sd_L3[t] ~ T(dt(0, 1/(0.5^2), 1), 0, )
    
    ### Scale size structure histogram parameters
    Midbins_LC3[1:LC, t] <- nf_midbins(mu = Mu_L3[t], sd = Sd_L3[t], nb_LC = LC, Lmin = Min_L3, Lmax = Max_L3)
    Prop_LC3[1:LC, t] <- nf_prop(mu = Mu_L3[t], sd = Sd_L3[t], nb_LC = LC, Lmin = Min_L3, Lmax = Max_L3)
  }
  
  ## Sex ratio parameters
  for (t in 1:C) {
    ### Probability of females
    Pfem_SR3[t] ~ dbeta(shape1 = 1, shape2 = 1)
    
    ### Probability of each sex
    Psex_SR3[t, 1:2] <- nimC(Pfem_SR3[t], (1 - Pfem_SR3[t]))
  }
  
  # Survival rate
  # -----------------------------
  ## Regression hyperparameters
  ### Slope
  Beta3 ~ dnorm(mean = 0, sd = 5)
  ### Intercept
  Alpha3 ~ dnorm(mean = 0, sd = 5)
  
  ## Survival rate logistic regression
  for (t in 1:C) {
    Logit_Theta3[1:LC, t] <- Alpha3 + Beta3 * (Midbins_LC3[1:LC, t] - Mean_L3)
    Theta3[1:LC, t] <- Min_Theta3 + ((Max_Theta3 - Min_Theta3) / (1 + exp(-Logit_Theta3[1:LC, t])))
  }
  
  # Population dynamics
  # -----------------------------
  ## Application of the survival rate
  for (t in 1:C) {
    for (s in 1:S) {
      N3[1:LC, t, s] <- Mu_N3[t] * Psex_SR3[t, s] * Prop_LC3[1:LC, t]
      N41[1:LC, t, s] <- N3[1:LC, t, s] * Theta3[1:LC, t]
      
      Tot_N41[t, s] <- sum(N41[1:LC, t, s])
    }
  }
  
  # Parameters of Stage 41
  # -----------------------------
  for (t in 1:C) {
    ## Scale size parameters
    Mu_L41[t] <- nf_omega(N = N41[1:LC, t, 1:S], midbins = Midbins_LC3[1:LC, t])[1]
    Sd_L41[t] <- nf_omega(N = N41[1:LC, t, 1:S], midbins = Midbins_LC3[1:LC, t])[2]
  }
  
  # --------------------------------------------------------
  ##                      LIKELIHOOD
  # --------------------------------------------------------
  
  # Scale size likelihood
  # -----------------------------
  for (t in 1:C) {
    ## For Stage 3
    for (l in 1:L3[t]) {
      Scales_L3[l, t] ~ dnorm(mean = Mu_L3[t], sd = Sd_L3[t])
    }
    
    ## For Stage 41
    for (l in 1:L41[t]) {
      Scales_L41[l, t] ~ dnorm(mean = Mu_L41[t], sd = Sd_L41[t])
    }
  }
  
  # Abundance pseudolikelihood
  # -----------------------------
  ## For Stage 3
  for (t in 1:C) {
    Eff_N3[t] ~ dlnorm(meanlog = Log_Mu_N3[t], sdlog = Log_Sd_N3[t])
  }
  
  # Sex ratio likelihood
  # -----------------------------
  ## For Stage 3
  for (t in 1:C) {
    Nsex_SR3[t, 1:2] ~ dmulti(prob = Psex_SR3[t, 1:2], size = Tot_SR3[t])
  }
  
  
  # -------------------------------------------------------------------------
  ##                          Stage 4 -> Stage 5 & 8
  # -------------------------------------------------------------------------
  
  # --------------------------------------------------------
  ##                          PRIOR
  # --------------------------------------------------------
  
  # Parameters of Stage 42
  # -----------------------------
  ## size structure parameters
  for (t in 1:C) {
    ### Scale size structure normal parameters
    Mu_L42[t] ~ dunif(min = 0, max = 10)
    Sd_L42[t] ~ T(dt(0, 1/(0.5^2), 1), 0, )
    
    ### Scale size structure histogram parameters
    Midbins_LC42[1:LC, t] <- nf_midbins(mu = Mu_L42[t], sd = Sd_L42[t], nb_LC = LC, Lmin = Min_L42, Lmax = Max_L42)
    Prop_LC42[1:LC, t] <- nf_prop(mu = Mu_L42[t], sd = Sd_L42[t], nb_LC = LC, Lmin = Min_L42, Lmax = Max_L42)
  }
  
  # Maturation probability
  # -----------------------------
  ## Regression hyperparameters
  ### Slope
  Beta4 ~ dnorm(mean = 0, sd = 5)
  ### Intercept
  Alpha4 ~ dnorm(mean = 0, sd = 5)
  
  ## Maturation rate logistic regression
  for (t in 1:C) {
    Logit_Theta4[1:LC, t] <- Alpha4 + Beta4 * (Midbins_LC42[1:LC, t] - Mean_L42)
    Theta4[1:LC, t] <- 1/(1+exp(-Logit_Theta4[1:LC, t]))
  }
  
  # Population dynamics
  # -----------------------------
  ## Application of the maturation rate
  for (t in 1:C) {
    for (s in 1:S) {
      N42[1:LC, t, s] <- Tot_N41[t, s] * Prop_LC42[1:LC, t]
      
      N5[1:LC, t, s] <- N42[1:LC, t, s] * Theta4[1:LC, t]
      N8[1:LC, t, s] <- N42[1:LC, t, s] * (1-Theta4[1:LC, t])
      
      Tot_N5[t, s] <- sum(N5[1:LC, t, s])
      Tot_N8[t, s] <- sum(N8[1:LC, t, s]) 
    }
  }
  
  # Parameters of Stage 5 and 8
  # -----------------------------
  ## Parameters of stage 5
  for (t in 1:C) {
    ### Scale size parameters
    Mu_L5[t] <- nf_omega(N = N5[1:LC, t, 1:S], midbins = Midbins_LC42[1:LC, t])[1]
    Sd_L5[t] <- nf_omega(N = N5[1:LC, t, 1:S], midbins = Midbins_LC42[1:LC, t])[2]
  }
  
  ## Parameters of stage 8
  for (t in 1:C) {
    ### Scale size parameters
    Mu_L8[t] <- nf_omega(N = N8[1:LC, t, 1:S], midbins = Midbins_LC42[1:LC, t])[1]
    Sd_L8[t] <- nf_omega(N = N8[1:LC, t, 1:S], midbins = Midbins_LC42[1:LC, t])[2]
  }
  
  # --------------------------------------------------------
  ##                      LIKELIHOOD
  # --------------------------------------------------------
  
  # Scale size likelihood
  # -----------------------------
  for (t in 1:C) {
    ## For Stage 42
    for (l in 1:L42[t]) {
      Scales_L42[l,t] ~ dnorm(mean = Mu_L42[t], sd = Sd_L42[t])
    }
    
    ## For Stage 5
    for (l in 1:L5[t]) {
      Scales_L5[l,t] ~ dnorm(mean = Mu_L5[t], sd = Sd_L5[t])
    }
    
    ## For Stage 8
    for (l in 1:L8[t]) {
      Scales_L8[l,t] ~ dnorm(mean = Mu_L8[t], sd = Sd_L8[t])
    }
  }
  
  
  # -------------------------------------------------------------------------
  ##                          Stage 5 & 8 -> Stage 6 & 9
  # -------------------------------------------------------------------------
  
  # --------------------------------------------------------
  ##                          PRIOR
  # --------------------------------------------------------
  
  # Parameters of Stage 6
  # -----------------------------
  ## size structure parameters
  for (t in 1:C) {
    ### Scale size structure normal parameters
    Mu_L6[t] ~ dunif(min = 0, max = 10)
    Sd_L6[t] ~ T(dt(0, 1/(0.5^2), 1), 0, )
    
    ### Scale size structure histogram parameters
    Prop_LC6[1:LC, t] <- nf_prop(mu = Mu_L6[t], sd = Sd_L6[t], nb_LC = LC, Lmin = Min_L6, Lmax = Max_L6)
  }
  
  # Parameters of Stage 9
  # -----------------------------
  ## size structure parameters
  for (t in 1:C) {
    ### Scale size structure normal parameters
    Mu_L9[t] ~ dunif(min = 0, max = 10)
    Sd_L9[t] ~ T(dt(0, 1/(0.5^2), 1), 0, )
    
    ### Scale size structure histogram parameters
    Prop_LC9[1:LC, t] <- nf_prop(mu = Mu_L9[t], sd = Sd_L9[t], nb_LC = LC, Lmin = Min_L9, Lmax = Max_L9)
  }
  
  
  # Population dynamics
  # -----------------------------
  ## Application of the fixed survival rate
  for (s in 1:S) {
    for (t in 1:C) {
      N6[1:LC, t, s] <- Tot_N5[t, s] * Prop_LC6[1:LC, t] * Theta5
      N9[1:LC, t, s] <- Tot_N8[t, s] * Prop_LC9[1:LC, t] * Theta8
      
      Tot_N6[t, s] <- sum(N6[1:LC, t, s])
      Tot_N9[t, s] <- sum(N9[1:LC, t, s])
    }
  }
  
  # Parameters 2 of Stage 6
  # -----------------------------
  for (t in 1:C) {
    ## Abundance parameters
    Mu_N6[t] <- Tot_N6[t, 1] + Tot_N6[t, 2]
    Log_Mu_N6[t] <- log(Mu_N6[t]) - ((1/2) * Log_Sd_N6[t])
    
    ## Sex ratio parameters
    Psex_SR6[t, 1] <- Tot_N6[t, 1]/Mu_N6[t]
    Psex_SR6[t, 2] <- Tot_N6[t, 2]/Mu_N6[t]
  }
  
  # Parameters 2 of Stage 9
  # -----------------------------
  for (t in 1:C) {
    ## Abundance parameters
    Mu_N9[t] <- Tot_N9[t, 1] + Tot_N9[t, 2]
    Log_Mu_N9[t] <- log(Mu_N9[t]) - ((1/2) * Log_Sd_N9[t])
    
    ## Sex ratio parameters
    Psex_SR9[t, 1] <- Tot_N9[t, 1]/Mu_N9[t]
    Psex_SR9[t, 2] <- Tot_N9[t, 2]/Mu_N9[t]
  }
  
  # --------------------------------------------------------
  ##                      LIKELIHOOD
  # --------------------------------------------------------
  
  # Scale size likelihood
  # -----------------------------
  for (t in 1:C) {
    ## For stage 6
    for (l in 1:L6[t]) {
      Scales_L6[l,t] ~ dnorm(mean = Mu_L6[t], sd = Sd_L6[t])
    }
    
    ## For stage 9
    for (l in 1:L9[t]) {
      Scales_L9[l,t] ~ dnorm(mean = Mu_L9[t], sd = Sd_L9[t])
    }
  }
  
  # Sex ratio likelihood
  # -----------------------------
  for (t in 1:C) {
    ## For stage 6
    Nsex_SR6[t, 1:2] ~ dmulti(prob = Psex_SR6[t, 1:2], size = Tot_SR6[t])
    
    ## For stage 9
    Nsex_SR9[t, 1:2] ~ dmulti(prob = Psex_SR9[t, 1:2], size = Tot_SR9[t])
  }
  
  # Abundance likelihood
  # -----------------------------
  for (t in 1:C) {
    ## For stage 6
    Eff_N6[t] ~ dlnorm(meanlog = Log_Mu_N6[t], sdlog = Log_Sd_N6[t])
    
    ## For stage 9
    Eff_N9[t] ~ dlnorm(meanlog = Log_Mu_N9[t], sdlog = Log_Sd_N9[t])
  }
  
  
  # --------------------------------------------------------
  ##                   END OF THE MODEL
  # --------------------------------------------------------
})