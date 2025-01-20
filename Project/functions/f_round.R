f_round <- function(x, digit = 1, method = 1) {
  ratio <- 10^digit
  
  if (method == 1) {
    return(floor(x * ratio) / ratio)  # Arrondit vers le bas
  } else if (method == 2) {
    return(ceiling(x * ratio) / ratio)  # Arrondit vers le haut
  } else {
    stop("Invalid method. Use 1 for floor or 2 for ceiling.")  # Message d'erreur si méthode invalide
  }
}