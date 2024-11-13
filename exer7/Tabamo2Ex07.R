# Tabamo, Euan Jed S. - CMSC 150 AB-3L
# Exercise 7 - Polynomial Interpolation
# Description: Solves problems using newton divided difference interpolation
# Date Created: November 2, 2024

# Prints a named list's elements
ListPrint <- function(lst, varName = "List") {
  cat(varName)
  cat(" =========================\n\n")
  for (name in names(lst)) {
    # Print the name of the element and the element itself
    cat(paste0(name, ":\n"))
    if (is.matrix(lst[[name]]) || is.function(lst[[name]])) {
      print(lst[[name]])
    } else {
      cat(paste(lst[[name]], collapse = "    "))
      cat("\n")
    }
    cat("\n")
  }
}

# Function to validate the input for Lagrange the interpolation functions
validate_input <- function(x, data) {
  if (!is.numeric(x)) {
    warning(paste0("Error: ", x, " is not numeric"))
    return(FALSE)
  }
  if (!is.list(data)) {
    warning(paste0("Error: ", data, " is not a list"))
    return(FALSE)
  }
  if (length(data) != 2) {
    warning(paste0("Error: ", data, " does not contain exactly two elements"))
    return(FALSE)
  }
  if (!is.numeric(data[[1]]) || !is.numeric(data[[2]])) {
    warning(paste0("Error: ", data, " does not contain numeric vectors"))
    return(FALSE)
  }
  if (length(data[[1]]) != length(data[[2]])) {
    warning(paste0("Vectors in data are not of the same length."))
    return(FALSE)
  }

  TRUE
}

# Predicts the value of the function given by the data at x.
# Returns the coefficients of the polynomial and the computed value f(x) = y.
NDD <- function(x, data) {
  # Check if data is invalid
  if (!validate_input(x, data)) {
    return(list(coefficients = NA, y = NA))
  }

  # Obtain x and y data and their length
  x_vector <- data[[1]]
  y_vector <- data[[2]]
  n <- length(x_vector)

  # Create an n×n matrix
  table <- matrix(NA, nrow = n, ncol = n)

  # Fill the first column with values from y = f(x)
  table[, 1] <- y_vector

  # Calculate the divided differences
  for (j in 2:n) { # Iterate over column 2 to the end
    for (i in 1:(n - j + 1)) { # Iterate over rows, decreasing per column
      table[i, j] <- (table[i + 1, j - 1] - table[i, j - 1]) /
        (x_vector[i + j - 1] - x_vector[i])
    }
  }

  # Get coefficients (first row contains all needed coefficients)
  coefficients <- table[1, 1:n]

  # Evaluate the polynomial at x
  # f(x) = b0 + b1(x - x0) + b2(x - x0)(x - x1) + ...
  y <- coefficients[1] # The first coefficient is the constant term
  product_term <- 1 # We can hold the previous product term and multiply it by (x - xi) for each iteration
  for (i in 2:n) {
    product_term <- product_term * (x - x_vector[i - 1])
    y <- y + coefficients[i] * product_term
  }

  # Return the list of coefficients and the computed value
  list(coefficients = coefficients, y = y)
}

# Sample Data
sample_data <- list(
  c(0, pi / 4, pi / 2, 3 * pi / 4), # 0, π/4, π/2, 3π/4
  round(sin(c(0, pi / 4, pi / 2, 3 * pi / 4)), 3) # 0, 0.707, 1, 0.707
)
sample_result <- NDD(0.5, sample_data)

# Turtle Data from Part 2
turtle_data <- list(
  c(1, 2, 3, 4, 5, 6),
  c(280, 310, 295, 700, 750, 720)
)
turtle_result <- NDD(3.5, turtle_data)

ListPrint(sample_result, "sin(0.5) approximation")
ListPrint(turtle_result, "Expected turtle population at 3.5 months (April 15)")