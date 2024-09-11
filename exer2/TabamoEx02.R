# Tabamo, Euan Jed S. - CMSC 150 AB-3L
# Exercise 2 - Matrices
# Description: Solves for the minors, cofactors, adjoint, and inverse of matrices.
# Date Created: September 11, 2024

# Returns a logical type based on whether the input matrix is a square matrix or not.
SquareMatrix <- function(mat) {
  # Ensure that the size of the rows is equal to the cols
  nrow(mat) == ncol(mat)
}

# Returns the determinant of a square matrix using recursion
# Returns NA if the matrix is not square
ComputeDeterminant <- function(mat) {
  # Check if the matrix is square
  if (!SquareMatrix(mat)) {
    return(NA)
  }

  n <- ncol(mat)

  sum <- 0
  # Base Case
  if (n == 1) {
    return(mat[1, 1])
  }
  # Recursive Case
  else {
    for (j in 1:n) {
      sum <- sum + (-1)^(1 + j) * mat[1, j] * ComputeDeterminant(as.matrix(mat[-1, -j]))
    }
  }

  # Return the determinant
  sum
}

# Returns the minor of the matrix given the position of an element in the matrix
# Returns NA if parameters are invalid
MatrixMinor <- function(mat, i, j) {
  # Check if the i and j parameters are in-bounds
  if (i < 1 || i > nrow(mat) || j < 1 || j > ncol(mat)) {
    return(NA)
  }

  # Return the determinant of the sub-matrix removing row i and j
  ComputeDeterminant(as.matrix(mat[-i, -j]))
}

# Returns the cofactor of the matrix given the position of an element in the matrix
# Returns NA if the parameters are invalid
MatrixCofactor <- function(mat, i, j) {
  # Check if the i and j parameters are in-bounds
  minor <- MatrixMinor(mat, i, j)
  if (is.na(minor)) {
    return(NA)
  }

  # Return the cofactor of the matrix at i and j
  (-1)^(i + j) * minor
}

# Returns the adjoint of the square matrix mat
# Returns NA if the matrix is nonsquare
MatrixAdjoint <- function(mat) {
  # Ensure that the matrix is square
  if (!SquareMatrix(mat)) {
    return(NA)
  }

  # Create a matrix of cofactors
  rows <- nrow(mat)
  cols <- ncol(mat)
  adjoint <- matrix(
    data = NA,
    nrow = rows,
    ncol = cols
  )
  for (i in 1:rows) {
    for (j in 1:cols) {
      adjoint[i, j] <- MatrixCofactor(mat, i, j)
    }
  }

  # Transpose the matrix of cofactors and return
  t(adjoint)
}

# Returns the inverse of a square matrix
# Returns NA if the matrix is not square or if the matrix is singular
MatrixInverse <- function(mat) {
  matrixDet <- ComputeDeterminant(mat)

  if (!SquareMatrix(A) || matrixDet == 0) {
    return(NA)
  }

  # Solve the inverse matrix and return
  (1 / matrixDet) * MatrixAdjoint(mat)
}

# Handout Example
cat("\n")
cat("HANDOUT EXAMPLE ---------------------------------------\n")
A <- matrix(c(1, -1, 1, -1, 2, 1, -1, 3, 4), nrow = 3, ncol = 3)
aSize <- nrow(A)

cat("Minors\n")
for (i in 1:aSize) {
  for (j in 1:aSize) {
    print(MatrixMinor(A, i, j))
  }
}

cat("Cofactors\n")
for (i in 1:aSize) {
  for (j in 1:aSize) {
    print(MatrixCofactor(A, i, j))
  }
}

cat("Adjoint\n")
print(MatrixAdjoint(A))

cat("Inverse\n")
print(MatrixInverse(A))

cat("Checking\n")
print(A %*% MatrixInverse(A))

# Exercise Matrix
cat("\n")
cat("EXERCISE PART B ---------------------------------------\n")
B <- matrix(c(2, 1, 4, -1, -3, -2, 0, -1, 5), nrow = 3, ncol = 3)
bSize <- nrow(B)
cat("Minors\n")
for (i in 1:bSize) {
  for (j in 1:bSize) {
    print(MatrixMinor(B, i, j))
  }
}

cat("Cofactors\n")
for (i in 1:bSize) {
  for (j in 1:bSize) {
    print(MatrixCofactor(B, i, j))
  }
}

cat("Adjoint\n")
print(MatrixAdjoint(B))

cat("Inverse\n")
print(MatrixInverse(B))

cat("Checking\n")
print(B %*% MatrixInverse(B))


# 2 x 2 Matrix
cat("\n")
cat("2 x 2 EXAMPLE -----------------------------------------\n")
C <- matrix(c(3, -1, 4, -2), nrow = 2, ncol = 2)
cSize <- nrow(C)

cat("Minors\n")
for (i in 1:cSize) {
  for (j in 1:cSize) {
    print(MatrixMinor(C, i, j))
  }
}

cat("Cofactors\n")
for (i in 1:cSize) {
  for (j in 1:cSize) {
    print(MatrixCofactor(C, i, j))
  }
}

cat("Adjoint\n")
print(MatrixAdjoint(C))

cat("Inverse\n")
print(MatrixInverse(C))

cat("Checking\n")
print(C %*% MatrixInverse(C))
