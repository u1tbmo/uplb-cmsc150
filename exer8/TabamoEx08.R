# Tabamo, Euan Jed S. - CMSC 150 AB-3L
# Exercise 8 - Muller's Method
# Description: Solves problems involving root finding using Muller's Method
# Date Created: November 13, 2024

# Validates the input of the MullerMethod function
InputValidator <- function(f = NA,
                           x0 = NA,
                           x1 = NA,
                           x2 = NA) {
  if (missing(f) || missing(x0) || missing(x1) || missing(x2)) {
    return(FALSE)
  }
  if (!is.function(f)) {
    return(FALSE)
  }
  if (!(is.numeric(x0) && is.numeric(x1) && is.numeric(x2))) {
    return(FALSE)
  }
  TRUE
}

# Returns a list :D
MullerMethod <- function(f = NA,
                         x0 = NA,
                         x1 = NA,
                         x2 = NA,
                         macheps = 0.00001,
                         maxiter = 1000,
                         verbose = TRUE) {
  if (!InputValidator(f, x0, x1, x2)) {
    if (!is.na(f)) {
      warning("Missing arguments!")
      return (list(
        f = NA,
        coefficients = NA,
        root = NA,
        iterations = NA,
        ea = NA
      ))
    } else {
      warning("Missing arguments!")
      return (list(
        f = f,
        coefficients = NA,
        root = NA,
        iterations = NA,
        ea = NA
      ))
    }
  }
  
  ea <- Inf
  iter <- 1L
  while (ea > macheps || iter < maxiter) {
    # Compute the function value y for each x
    y0 <- f(x0)
    y1 <- f(x1)
    y2 <- f(x2)
    
    # Compute for h0, d0, h1, and d1
    h0 <- x1 - x0
    h1 <- x2 - x1
    d0 <- f(x1) - f(x0) / h0
    d1 <- f(x2) - f(x1) / h1
    
    # Compute the coefficients and constant of the derived quadratic polynomial
    A <- (d1 - d0) / (h1 + h0)
    B <- A * h1 + d1
    C <- f(x2)
    
    # Check the sign of the denominator of the alternative quadratic formula
    if ()
  }
}
