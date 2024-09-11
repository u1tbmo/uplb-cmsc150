lu_decomposition <- function(matrix_a) {
  # Get the dimensions of matrix_a
  n <- nrow(matrix_a)

  # Initialize l as the identity matrix and u as a copy of a
  matrix_l <- diag(x = 1, nrow = n, ncol = n)
  matrix_u <- matrix_a

  # Perform forward elimination
  for (i in 1:(n - 1)) {
    if (matrix_u[i, i] == 0) {
      stop("Pivot element is 0. LU decomposition cannot proceed.")
    }
    # Loop over rows below the current row i
    for (j in (i + 1):n) {
      # Find the multiplier
      multiplier <- matrix_u[j, i] / matrix_u[i, i]

      # Update matrix_l matrix
      matrix_l[j, i] <- multiplier

      # Update matrix_u matrix
      matrix_u[j, ] <- matrix_u[j, ] - multiplier * matrix_u[i, ]
    }
  }

  list(L = matrix_l, U = matrix_u)
}

forward_elim_determinant <- function(matrix_u) {
  # Get the dimensions of matrix_a
  n <- nrow(matrix_u)

  # Calculate determinant
  determinant <- 1
  for (i in 1:n) {
    determinant <- determinant * matrix_u[i, i]
  }

  determinant
}

# Example

matrix_a <- matrix(
  c(3, 0.1, 0.3, -0.1, 7, -0.2, -0.2, -0.3, 10),
  nrow = 3,
  ncol = 3
)

cat("LU Decomposition: \n")
triangles <- lu_decomposition(matrix_a)
matrix_l <- triangles$L
matrix_u <- triangles$U
print(matrix_l)
print(matrix_u)
cat("\nDeterminant: ")
print(forward_elim_determinant(matrix_u))
cat("\n")
