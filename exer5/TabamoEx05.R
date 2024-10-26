# Tabamo, Euan Jed S. - CMSC 150 AB-3L
# Exercise 5 - Simplex Method
# Description: Solves linear programming questions using the Simplex Method
# Date Created: October 9, 2024

# Prints a named list's elements
listPrint <- function(lst, varName = "List") {
  cat("==============================\n\n")
  for (name in names(lst)) {
    cat("-- ")
    cat(name)
    cat(" --\n")
    if (is.matrix(lst[[name]])) {
      print(lst[[name]])
    } else {
      cat(paste(lst[[name]], collapse = "    "))
      cat("\n")
    }
    cat("\n")
  }
}

# Calculates the magnitude of a vector
# If the vector's magnitude is 1, then the vector is clear
Magnitude <- function(vector) {
  # Check if vector is a vector
  if (!is.vector(vector)) {
    cat("Invalid argument: vector must be a vector!\n")
    return(NA)
  }

  sum <- 0L
  # Iterate over the vector
  for (num in vector) {
    sum <- sum + num^2
  }

  sqrt(sum)
}

# Returns a matrix of the basic solution to the linear programming problem
BasicSolution <- function(finalTableau, isMax = TRUE) {
  # Check if isMax exists
  if (!exists("isMax")) {
    cat("Missing argument: isMax\n")
    return(NA)
  }

  # Check if finalTableau is a matrix
  if (!is.matrix(finalTableau)) {
    cat("Invalid argument: finalTableau must be a matrix!")
    return(NA)
  }

  # Get the size of finalTableau
  m <- nrow(finalTableau)
  n <- ncol(finalTableau)

  # Get only n-1 colnames of finalTableau
  col_names <- head(colnames(finalTableau), n - 1)

  # Create a matrix of the basic solution
  soln <- matrix(data = 0, nrow = 1, ncol = n - 1, dimnames = list("", col_names))

  if (isMax) {
    # Case 1: Maximization
    # Iterate over every column and check if the column is clear or not
    for (i in 1:(n - 1)) {
      # If the magnitude is 1, then the column is clear
      if (Magnitude(finalTableau[, i]) == 1) {
        for (j in 1:m) {
          if (finalTableau[j, i] == 1) {
            soln[1, i] <- finalTableau[j, n]
          }
        }
      }
    }
  } else {
    # Case 2: Minimization
    # For minimization, the bottom row except for the minimized value is the bottom row
    for (i in 1:(n - 1)) {
      if (i != n - 1) {
        soln[1, i] <- finalTableau[m, i]
      } else {
        soln[1, i] <- finalTableau[m, n]
      }
    }
  }

  # Return the soln matrix
  soln
}

# Returns a labelled list with the final tableau, basic solution,
# and the maximum or minimum value derived from the simplex method
Simplex <- function(tableau, isMax = TRUE) {
  # Check if isMax exists
  if (!is.logical(isMax)) {
    cat("Invalid argument: isMax must be a logical value!\n")
    return(list(finalTableau = NA, basicSolution = NA, Z = NA))
  }

  # Check if tableau is a matrix
  if (!is.matrix(tableau)) {
    cat("Invalid argument: tableau must be a matrix!\n")
    return(list(finalTableau = NA, basicSolution = NA, Z = NA))
  }


  # Get the number of rows and cols in the tableau
  m <- nrow(tableau)
  n <- ncol(tableau)

  # Logic variable for loop condition
  soln_found <- FALSE

  # Repeat until solution is found
  while (!soln_found) {
    # Obtain the indices of the negative numbers in the objective function
    p_col_idx <- which.min(tableau[m, 1:(n - 1)])

    # The solution is found if the smallest number is nonnegative
    if (tableau[m, p_col_idx] >= 0) {
      soln_found <- TRUE
      break
    }

    # Find the smallest positive test ratio in the column
    p_row_idx <- NA
    smallest_tr <- Inf

    # Iterate through every row except the objective function row
    for (i in 1:(m - 1)) {
      # If the divisor for the test ratio is nonpositive, then ignore the row
      if (tableau[i, p_col_idx] < 0) {
        next
      }

      # Calculate the test ratio for the current row
      tr <- tableau[i, n] / tableau[i, p_col_idx]

      # Check if the test ratio is the smallest found so far
      if (tr < smallest_tr) {
        # Update the smallest test ratio and index of the pivot row
        smallest_tr <- tr
        p_row_idx <- i
      }
    }

    # If no valid pivot row is found, the problem is unbounded
    if (is.na(p_row_idx)) {
      return(list(finalTableau = tableau, basicSolution = NA, Z = NA))
    }

    # For maximization, indicate the entering variable
    # This is not necessary for minimization since the last element of the column itself is the solution
    if (isMax) {
      rownames(tableau)[p_row_idx] <- colnames(tableau)[p_col_idx]
    }

    # Gauss-Jordan Elimination
    # Normalize the row
    tableau[p_row_idx, ] <- tableau[p_row_idx, ] / tableau[p_row_idx, p_col_idx]
    for (i in 1:m) {
      # Row Elimination
      if (i != p_row_idx) {
        tableau[i, ] <- tableau[i, ] - tableau[i, p_col_idx] * tableau[p_row_idx, ]
      }
    }
  }

  # Get the solution from the tableau
  soln <- BasicSolution(tableau, isMax)

  # Return the named list
  list(finalTableau = tableau, basicSolution = soln, Z = soln[1, n - 1])
}

# Exercise Part 1
exer_matrix <- matrix(
  c(
    1, 1, 1, 0, 0, 0, 95,
    3, 7, 0, 1, 0, 0, 400,
    1, 2, 0, 0, 1, 0, 120,
    -4000, -5000, 0, 0, 0, 1, 0
  ),
  byrow = TRUE,
  nrow = 4,
  dimnames = list(
    c("s1", "s2", "s3", "Z"),
    c("x1", "x2", "s1", "s2", "s3", "P", "Solution")
  )
)

# Maximization Example
max_sample_matrix <- matrix(
  c(
    1, 0, 1, 0, 0, 6, # x + s1 = 6
    3, 1, 0, 1, 0, 9, # 3x + y + s2 = 9
    -3, -4, 0, 0, 1, 0 # -3x -4y + Z = 0
  ),
  byrow = TRUE,
  nrow = 3,
  dimnames = list(
    c("s1", "s2", "Z"),
    c("x", "x2", "s1", "s2", "Z", "Solution")
  )
)

# Minimization Example
#   Tableau must be the dual problem tableau
min_sample_matrix <- matrix(
  c(
    1, 7, 1, 0, 0, 14, # s1 + 7s2 + x1 = 14
    2, 6, 0, 1, 0, 20, # 2s1 + 6s2 + x2 = 20
    -4, -20, 0, 0, 1, 0 # -4s1 -20s2 + Z = 0
  ),
  byrow = TRUE,
  nrow = 3,
  dimnames = list(
    rep("", 3),
    c("s1", "s2", "x1", "x2", "Z", "Solution")
  )
)

# Unbounded Example
# Occurs when:
#   Optimal solution still not found
#   A pivot row cannot be selected (all test ratios are nonpositive)
unbounded_matrix <- matrix(
  c(
    1, -1, 1, 0, 0, 10, # x1 - x2 + s1 = 10
    2, -1, 0, 1, 0, 40, # 2x1 - x2 + s2 = 40
    -2, -1, 0, 0, 1, 0 # -2x1 - x2 + Z = 0
  ),
  byrow = TRUE,
  nrow = 3,
  dimnames = list(
    rep("", 3),
    c("x1", "x2", "s1", "s2", "Z", "Solution")
  )
)

# Degeneracy Example
# Occurs when:
#   A tie occurs when computing for the test ratio
#   The RHS has a 0 in the BFS
degeneracy_matrix <- matrix(
  c(
    4, 3, 1, 0, 0, 0, 12,
    4, 1, 0, 1, 0, 0, 8,
    4, 2, 0, 0, 1, 0, 8,
    -2, -1, 0, 0, 0, 1, 0
  ),
  byrow = TRUE,
  nrow = 4,
  dimnames = list(
    c("s1", "s2", "s3", "Z"),
    c("x1", "x2", "s1", "s2", "s3", "Z", "Solution")
  )
)

degeneracy_matrix_2 <- matrix(
  c(
    2, 1, 1, 0, 0, 0, 4,   # 2x1 + x2 + s1 = 4
    1, 2, 0, 1, 0, 0, 4,   # x1 + 2x2 + s2 = 4
    1, 1, 0, 0, 1, 0, 2,   # x1 + x2 + s3 = 2
    -3, -2, 0, 0, 0, 1, 0  # -3x1 - 2x2 + Z = 0
  ),
  byrow = TRUE,
  nrow = 4,
  dimnames = list(
    c("s1", "s2", "s3", "Z"),
    c("x1", "x2", "s1", "s2", "s3", "Z", "Solution")
  )
)

result_exer <- Simplex(exer_matrix, TRUE)
result_max <- Simplex(max_sample_matrix, TRUE)
result_min <- Simplex(min_sample_matrix, FALSE)
result_unbounded <- Simplex(unbounded_matrix, TRUE)
result_degeneracy <- Simplex(degeneracy_matrix, TRUE)
result_degeneracy_2 <- Simplex(degeneracy_matrix_2, TRUE)
listPrint(result_exer)
listPrint(result_max)
listPrint(result_min)
listPrint(result_unbounded)
listPrint(result_degeneracy)
listPrint(result_degeneracy_2)
