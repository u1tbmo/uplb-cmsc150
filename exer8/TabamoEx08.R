# Tabamo, Euan Jed S. - CMSC 150 AB-3L
# Exercise 8 - Muller's Method
# Description: Solves problems involving root finding using Muller's Method
# Date Created: November 13, 2024

# Validates the input of the MullerMethod function
InputValidator <- function(f = NA, x0 = NA, x1 = NA, x2 = NA) {
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

# Converts a complex number to a real number if possible
ComplexToReal <- function(z) {
  # If z is complex and the imaginary part is 0, return the real part
  # Otherwise, return z
  if (is.complex(z) && Im(z) == 0) Re(z) else z
}

# Returns a list :D
MullerMethod <- function(f = NA, x0 = NA, x1 = NA, x2 = NA, macheps = 0.00001, maxiter = 1000, verbose = TRUE) {
  if (!InputValidator(f, x0, x1, x2)) {
    if (!is.na(f)) {
      warning("Missing arguments!")
      return(list(
        f = NA,
        coefficients = NA,
        root = NA,
        iterations = NA,
        ea = NA
      ))
    } else {
      warning("Missing arguments!")
      return(list(
        f = f,
        coefficients = NA,
        root = NA,
        iterations = NA,
        ea = NA
      ))
    }
  }

  table <- data.frame(
    x0 = numeric(),
    x1 = numeric(),
    x2 = numeric(),
    f_x0 = numeric(),
    f_x1 = numeric(),
    f_x2 = numeric(),
    A = numeric(),
    B = numeric(),
    C = numeric(),
    x3 = numeric(),
    f_x3 = numeric(),
    Error = numeric()
  )

  ea <- Inf
  iter <- 0L
  while (ea > macheps && iter < maxiter) {
    # Compute the function value y for each x
    y0 <- f(x0)
    y1 <- f(x1)
    y2 <- f(x2)

    # Compute for h0, d0, h1, and d1
    h0 <- x1 - x0
    h1 <- x2 - x1
    d0 <- (f(x1) - f(x0)) / h0
    d1 <- (f(x2) - f(x1)) / h1

    # Compute the coefficients and constant of the derived quadratic polynomial
    A <- (d1 - d0) / (h1 + h0)
    B <- A * h1 + d1
    C <- f(x2)

    # Check the sign of the denominator of the alternative quadratic formula
    # Conversion to complex before taking the root allows taking the root of negative numbers
    denominator <- NA
    denominatorPositive <- B + sqrt(as.complex(B^2 - 4 * A * C))
    denominatorNegative <- B - sqrt(as.complex(B^2 - 4 * A * C))

    # Choose the denominator with the larger magnitude
    if (abs(denominatorNegative) < abs(denominatorPositive)) {
      denominator <- denominatorPositive
    } else {
      denominator <- denominatorNegative
    }

    # Compute for x3 and the approximate relative error
    x3 <- x2 - 2 * C / denominator
    ea <- abs((x3 - x2) / x3)

    # Append the values to the table
    table <- rbind(
      table,
      # Convert complex numbers to real if possible
      data.frame(
        x0 = ComplexToReal(x0),
        x1 = ComplexToReal(x1),
        x2 = ComplexToReal(x2),
        f_x0 = ComplexToReal(y0),
        f_x1 = ComplexToReal(y1),
        f_x2 = ComplexToReal(y2),
        A = ComplexToReal(A),
        B = ComplexToReal(B),
        C = ComplexToReal(C),
        x3 = ComplexToReal(x3),
        f_x3 = ComplexToReal(f(x3)),
        Error = ComplexToReal(ea)
      )
    )

    # Update the values of x0, x1, x, and iter for the next iteration
    x0 <- x1
    x1 <- x2
    x2 <- x3
    iter <- iter + 1L
  }

  # Print the table if verbose is TRUE
  if (verbose) print(table)

  # Return the results as a list
  list(
    f = f,
    coefficients = sapply(c(A, B, C), ComplexToReal),
    root = ComplexToReal(x3),
    iterations = iter,
    ea = ComplexToReal(ea)
  )
}

# Real roots sample
fxn1 <- function(x) x^3 - 13 * x - 12
result1 <- MullerMethod(fxn1, 4.5, 5.5, 5, 10^-09)

# Complex roots sample
fxn2 <- function(x) x^3 - x^2 + 2 * x - 2
result2 <- MullerMethod(fxn2, 0, 3, 5, 1e-09, 1000, TRUE)

# Cosine function sample
