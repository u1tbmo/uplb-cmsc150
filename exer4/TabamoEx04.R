# Tabamo, Euan Jed S. - CMSC 150 AB-3L
# Exercise 4 - Systems of Linear Equations
# Description: Solves systems of linear equations using
#              Gaussian and Gauss-Jordan Elimination.
# Date Created: September 18, 2024

# Library
library(MASS)

# Prints a named list's elements
listPrint <- function(lst) {
  for (name in names(lst)) {
    cat("-- ")
    cat(name)
    cat(" --------------------\n")
    if (is.matrix(lst[[name]])) {
      print(lst[[name]])
    } else {
      cat(paste(lst[[name]], collapse = "    "), "\n")
    }
    cat("\n\n")
  }
}

# Prints a named list's elements with fractional values
listPrintFractional <- function(lst) {
  for (name in names(lst)) {
    cat("-- ")
    cat(name)
    cat(" --------------------\n")
    if (is.matrix(lst[[name]])) {
      print(fractions(lst[[name]]))
    } else {
      cat(paste(lst[[name]], collapse = "    "), "\n")
    }
    cat("\n\n")
  }
}

# Uses the Gaussian Elimination method to solve a system of linear equations
# given the augmented coefficient matrix and the variables
GaussianMethod <- function(input) {
  # Check if the input is a list and has the required elements
  if (!is.list(input)) {
    cat("Invalid argument: input must be a list.\n")
    return(NA)
  }
  for (name in c("augcoeffmatrix", "variables")) {
    if (!name %in% names(input)) {
      cat("Missing argument: ", name, "\n")
      return(NA)
    }
  }

  # Get the input variables
  augcoeffmatrix <- input$augcoeffmatrix
  variables <- input$variables

  # Get the row count of the augmented coefficient matrix
  n <- nrow(augcoeffmatrix)

  # Check if arguments are valid
  if (length(variables) != nrow(augcoeffmatrix)) {
    stop(
      "The number of variables do not match the number of rows in the matrix."
    )
  }
  if (nrow(augcoeffmatrix) != ncol(augcoeffmatrix) - 1) {
    stop(
      "The augmented coefficient matrix size is invalid."
    )
  }

  # Initialize needed variables
  u <- augcoeffmatrix
  x <- numeric(n)

  # Forward Elimination

  # Iterate over every row except the last (pivot rows)
  for (i in 1:(n - 1)) {
    # Get the pivot row
    pivot_row <- which.max((abs(u[i:n, i]))) + (i - 1)
    # If the pivot element is 0, there is no solution
    if (u[pivot_row, i] == 0) {
      return(NA)
    }

    # Swap rows
    temp <- u[pivot_row, ]
    u[pivot_row, ] <- u[i, ]
    u[i, ] <- temp

    # Iterate over every row below the pivot row (eval rows)
    for (j in (i + 1):n) {
      # Find the multiplier and multiply by the pivot row
      # to get the normalized row, then subtract the normalized row
      # from the eval row and assign to the eval row.
      multiplier <- u[j, i] / u[i, i]
      normalized_row <- multiplier * u[i, ]
      u[j, ] <- u[j, ] - normalized_row
    }
  }


  # Backward Elimination
  x[n] <- u[n, n + 1] / u[n, n]
  for (i in (n - 1):1) {
    x[i] <- (u[i, n + 1] - sum(u[i, (i + 1):n] * x[(i + 1):n])) / u[i, i]
  }

  # Return the variables, augmented coefficient matrix, and solutions.
  list(
    variables = variables,
    augcoeffmatrix = u,
    solution = x
  )
}

GaussJordanMethod <- function(input) {
  # Check if the input is a list and has the required elements
  if (!is.list(input)) {
    cat("Invalid argument: input must be a list.\n")
    return(NA)
  }
  for (name in c("augcoeffmatrix", "variables")) {
    if (!name %in% names(input)) {
      cat("Missing argument: ", name, "\n")
      return(NA)
    }
  }

  # Get the input variables
  augcoeffmatrix <- input$augcoeffmatrix
  variables <- input$variables

  # Get the row count of the augmented coefficient matrix
  n <- nrow(augcoeffmatrix)

  # Check if arguments are valid
  if (length(variables) != nrow(augcoeffmatrix)) {
    stop(
      "The number of variables do not match the number of rows in the matrix."
    )
  }
  if (nrow(augcoeffmatrix) != ncol(augcoeffmatrix) - 1) {
    stop(
      "The augmented coefficient matrix size is invalid."
    )
  }

  # Initialize needed variables
  a <- augcoeffmatrix
  x <- numeric(n)

  # Iterate over every row
  for (i in 1:n) {
    if (i != n) {
      # Get the pivot row
      pivot_row <- which.max((abs(a[i:n, i]))) + (i - 1)
      # If the pivot element is 0, there is no solution
      if (a[pivot_row, i] == 0) {
        return(NA)
      }

      # Swap rows
      temp <- a[pivot_row, ]
      a[pivot_row, ] <- a[i, ]
      a[i, ] <- temp
    }
    # Scale the pivot row so that the pivot element is 1
    a[i, ] <- a[i, ] / a[i, i]

    # Iterate over every row other than the pivot row and normalize
    for (j in 1:n) {
      if (i == j) {
        next
      }
      normalized_row <- a[j, i] * a[i, ]
      a[j, ] <- a[j, ] - normalized_row
    }
  }

  # Get the solutions
  for (i in 1:n) {
    x[i] <- a[i, n + 1]
  }

  # Return the variables, augmented coefficient matrix, and solutions.
  list(
    variables = variables,
    augcoeffmatrix = a,
    solution = x
  )
}

# Handout Sample
handout_matrix <- matrix(
  c(8, 0.6, 20, -3, -0.3, 0, 0.9, 0.1, 1, 12, -45, 46),
  nrow = 3,
  ncol = 4,
  dimnames = list(
    c("r1", "r2", "r3"), c("c1", "c2", "c3", "b")
  )
)

result1a <- GaussianMethod(
  list(
    augcoeffmatrix = handout_matrix,
    variables = c("x1", "x2", "x3")
  )
)

result1b <- GaussJordanMethod(
  list(
    augcoeffmatrix = handout_matrix,
    variables = c("x1", "x2", "x3")
  )
)

# Print results
listPrintFractional(result1a)
listPrintFractional(result1b)

# Exercise Part 2 Item 1
part2a_matrix <- matrix(
  c(
    8000, 4500, 4000, 3000, 2000, 1000, 900, 250, 143145000,
    7800, 6500, 5800, 0, 3100, 1600, 1000, 300, 158870000,
    10000, 0, 3100, 0, 2600, 1300, 850, 150, 108440000,
    5200, 3700, 3100, 2700, 2400, 1800, 1200, 450, 143805000,
    7700, 7100, 0, 5700, 5100, 1300, 950, 95, 181390500,
    9300, 8700, 6100, 5100, 4000, 1000, 700, 70, 209273000,
    6000, 0, 5000, 4300, 3000, 1900, 1400, 920, 174388000,
    8500, 3700, 4200, 3900, 3500, 2400, 1000, 250, 183065000
  ),
  nrow = 8,
  ncol = 9,
  dimnames = list(
    c("A", "B", "C", "D", "E", "F", "G", "H"),
    c("Gen Ad B", "Gen Ad A", "Upper Box B", "Upper Box A", "Lower Box B", "Lower Box A", "VIP", "Royalty", "Total Profit")
  ),
  byrow = TRUE
)
result2a1 <- GaussianMethod(
  list(
    augcoeffmatrix = part2a_matrix,
    variables = c("Gen Ad B", "Gen Ad A", "Upper Box B", "Upper Box A", "Lower Box B", "Lower Box A", "VIP", "Royalty")
  )
)
result2a2 <- GaussJordanMethod(
  list(
    augcoeffmatrix = part2a_matrix,
    variables = c("Gen Ad B", "Gen Ad A", "Upper Box B", "Upper Box A", "Lower Box B", "Lower Box A", "VIP", "Royalty")
  )
)

# Print results
listPrint(result2a1)
listPrint(result2a2)

# Exercise Part 2 Item 2
part2b_matrix <- matrix(
  c(
    4, -1, 0, -1, 0, 0, 0, 0, 0, 80,
    -1, 4, -1, 0, -1, 0, 0, 0, 0, 30,
    0, -1, 4, 0, 0, -1, 0, 0, 0, 80,
    -1, 0, 0, 4, -1, 0, -1, 0, 0, 50,
    0, -1, 0, -1, 4, -1, 0, -1, 0, 0,
    0, 0, -1, 0, -1, 4, 0, 0, -1, 50,
    0, 0, 0, -1, 0, 0, 4, -1, 0, 120,
    0, 0, 0, 0, -1, 0, -1, 4, -1, 70,
    0, 0, 0, 0, 0, -1, 0, -1, 4, 120
  ),
  nrow = 9,
  ncol = 10,
  dimnames = list(
    c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
    c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "b")
  ),
  byrow = TRUE
)
result2b1 <- GaussianMethod(
  list(
    augcoeffmatrix = part2b_matrix,
    variables = c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9")
  )
)
result2b2 <- GaussJordanMethod(
  list(
    augcoeffmatrix = part2b_matrix,
    variables = c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9")
  )
)

# Print results
listPrint(result2b1)
listPrint(result2b2)
