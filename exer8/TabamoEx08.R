# Tabamo, Euan Jed S. - CMSC 150 AB-3L
# Exercise 8 - Muller's Method
# Description: Solves problems involving root finding using Muller's Method
# Date Created: November 13, 2024

# Prints a named list's elements
ListPrint <- function(lst, varName = "List") {
  cat(varName)
  cat(" =========================\n\n")
  for (name in names(lst)) {
    # Print the name of the element and the element itself
    cat(paste0(name, ":\n"))
    if (is.matrix(lst[[name]]) || is.function(lst[[name]])) {
      # If the element is a matrix or a function, print it
      print(lst[[name]])
    } else if (is.numeric(lst[[name]]) || is.complex(lst[[name]])) {
      # If the element is numeric or complex, use cat and respect the digits option
      cat(format(lst[[name]], digits = getOption("digits")))
      cat("\n")
    } else {
      # Otherwise, use cat and collapse the elements with spaces
      cat(paste(lst[[name]], collapse = "    "))
      cat("\n")
    }
    cat("\n")
  }
}

# Validates the input of the MullerMethod function
InputValidator <- function(f, x0, x1, x2, macheps = 0.00001, maxiter = 1000, verbose = TRUE) {
  if (missing(f) || missing(x0) || missing(x1) || missing(x2)) {
    warning(" : Missing arguments!")
    return(FALSE)
  }
  if (!is.function(f)) {
    warning("Argument 'f' is not a function!")
    return(FALSE)
  }
  if (!(is.numeric(x0) && is.numeric(x1) && is.numeric(x2))) {
    warning("Arguments 'x0', 'x1', and 'x2' must be numeric!")
    return(FALSE)
  }
  if (!is.numeric(macheps) || macheps <= 0) {
    warning("Argument 'macheps' must be a positive numeric value!")
    return(FALSE)
  }
  if (!is.numeric(maxiter) || maxiter <= 0 || maxiter != as.integer(maxiter)) {
    warning("Argument 'maxiter' must be a positive integer!")
    return(FALSE)
  }
  if (!is.logical(verbose)) {
    warning("Argument 'verbose' must be a logical value!")
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

# Returns a list containing the results of Muller's Method
MullerMethod <- function(f, x0, x1, x2, macheps = 0.00001, maxiter = 1000, verbose = TRUE) {
  # Validate the input
  if (!InputValidator(f, x0, x1, x2, macheps, maxiter, verbose)) {
    return(
      list(
        f = if (!is.function(f)) NA else f,
        coefficients = NA,
        root = NA,
        iterations = NA,
        ea = NA
      )
    )
  }

  # Initialize the variables for the loop
  table <- data.frame()
  ea <- Inf
  iter <- 0L
  while (ea == "undefined" || (ea > macheps && iter < maxiter)) {
    # Compute the function value y for each x
    y0 <- f(x0)
    y1 <- f(x1)
    y2 <- f(x2)

    # Compute for h0, h1, which are the step sizes between x0, x1, and x2
    h0 <- x1 - x0
    h1 <- x2 - x1
    # Compute d0, d1, which are the divided differences that approximate the slope
    d0 <- (f(x1) - f(x0)) / h0
    d1 <- (f(x2) - f(x1)) / h1

    # Compute the coefficients and constant of the derived quadratic polynomial
    A <- (d1 - d0) / (h1 + h0) # Rate of change of the slope between intervals
    B <- A * h1 + d1 # Slope at x2
    C <- y2 # Function value at x2

    # Calculate the discriminant, allowing for complex numbers
    discriminant <- ComplexToReal(sqrt(as.complex(B^2 - 4 * A * C)))

    # Determine which denominator to use based on magnitude and compute x3
    if (abs(B + discriminant) > abs(B - discriminant)) {
      x3 <- x2 - (2 * C) / (B + discriminant)
    } else {
      x3 <- x2 - (2 * C) / (B - discriminant)
    }

    # Calculate the approximate relative error and account for x3 being 0
    ea <- abs((x3 - x2) / x3)

    # Append the values to the table
    table <- rbind(
      table,
      # Convert complex numbers to real if possible before storing in the data frame
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
        Error = ea
      )
    )

    # Update the values of x0, x1, x, and iter for the next iteration
    x0 <- x1
    x1 <- x2
    x2 <- x3
    iter <- iter + 1L
  }

  # Print the table if verbose is TRUE
  if (verbose) {
    cat("Iterations:\n\n")
    print(table)
    cat("\n")
  }

  # Return the results as a list, zapsmall any very small numbers
  list(
    f = f,
    coefficients = sapply(c(A, B, C), function(x) zapsmall(ComplexToReal(x))), # Apply ComplexToReal and zapsmall to the coefficients
    root = zapsmall(ComplexToReal(x3)),
    iterations = iter,
    ea = zapsmall(ea)
  )
}

# Digit option
options(digits = 4)

# Real roots sample
fxn_real <- function(x) x^3 - 13 * x - 12
ListPrint(MullerMethod(fxn_real, 4.5, 5.5, 5, 10^-09, 1000, TRUE), "Real Roots Results")

# Complex roots sample
fxn_complex <- function(x) x^3 - x^2 + 2 * x - 2
ListPrint(MullerMethod(fxn_complex, 0, 3, 5, 1e-09, 1000, TRUE), "Complex Roots Results")

# Cosine function sample
ListPrint(MullerMethod(function(x) cos(x), 0, 2, 4, 1e-05, 1000, TRUE), "Cosine Function Results")
