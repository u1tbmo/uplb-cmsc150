# Tabamo, Euan Jed S. - CMSC 150 AB-3L
# Exercise 3 - LU Decomposition
# Description: Finds the LU Decomposition of matrices.
# Date Created: September 18, 2024

# Library
library(MASS)

# Returns a logical type based on whether the input matrix is a square matrix or not.
SquareMatrix <- function(mat) {
  # Ensure that the size of the rows is equal to the cols
  nrow(mat) == ncol(mat)
}

# Returns the list containing the resulting upper and lower triangular matrices
# of an LU Decomposition
LUDecomposition <- function(mat) {
  # Ensure the matrix is square
  if (!SquareMatrix(mat)) {
    return(NA)
  }

  # Dimensions of the matrix
  n <- nrow(mat)

  # Set up the L and U matrices
  matrix_l <- diag(x = 1, nrow = n)
  matrix_u <- mat

  # Iterate over the pivot rows
  for (i in 1:(n - 1)) {
    # If the matrix is singular, LU decomposition cannot proceed
    if (matrix_u[i, i] == 0) {
      stop("Pivot element is 0! LU Decomposition is not possible.")
    }
    # Iterate over the eval rows
    for (j in (i + 1):n) {
      # Find the multiplier and put it into the L matrix
      multiplier <- matrix_u[j, i] / matrix_u[i, i]
      matrix_l[j, i] <- multiplier
      # Perform the row operations for the U matrix
      matrix_u[j, ] <- matrix_u[j, ] - multiplier * matrix_u[i, ]
    }
  }

  # Return the list containing the original matrix, the resulting
  # L and U matrices, and the result of the LU matrix multiplication
  # A = LU is always true
  list(
    A = fractions(mat),
    L = fractions(matrix_l),
    U = fractions(matrix_u),
    LU = fractions(matrix_l %*% matrix_u)
  )
}

# Prints a named list's elements
listPrint <- function(lst) {
  for (name in names(lst)) {
    cat("-- ")
    cat(name)
    cat(" --------------------\n")
    print(lst[[name]])
    cat("\n\n")
  }
}

# Example
matrix_a <- matrix(
  c(2, 2, -1, 3 , 7, 2, 2, 3, 8),
  nrow = 3,
  ncol = 3,
  byrow = TRUE
)
listPrint(LUDecomposition(matrix_a))

# Handout
matrix_b <- matrix(
  c(3, -0.1, -0.2, 0.1, 7, -0.3, 0.3, -0.2, 10),
  nrow = 3,
  ncol = 3,
  byrow = TRUE
)
listPrint(LUDecomposition(matrix_b))

# Part 2
matrix_c <- matrix(
  c(9, 3, 2, 2, 7, 3, -1, 2, 8),
  nrow = 3,
  ncol = 3,
  byrow = FALSE
)
listPrint(LUDecomposition(matrix_c))
