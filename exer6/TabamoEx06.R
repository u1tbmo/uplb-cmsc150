# Tabamo, Euan Jed S. - CMSC 150 AB-3L
# Exercise 6 - Linear and Polynomial Regression
# Description: Solves problems using linear and polynomial regression
# Date Created: October 23, 2024

# Prints a named list's elements
ListPrint <- function(lst, varName = "List") {
  cat("==============================\n\n")
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

# Returns the solution vector of an input augmented coefficient matrix
GaussJordanMethod <- function(augcoeffmatrix) {
  # Get the row count of the augmented coefficient matrix
  n <- nrow(augcoeffmatrix)

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

  # Return the solution
  x
}

# Returns a list containing the augmented coefficient matrix, coefficients,
# polynomial string, and polynomial function of the polynomial regression
PolynomialRegression <- function(order, matrix) {
  # Check if matrix is a matrix and has the correct dimensions
  if (!is.matrix(matrix) || ncol(matrix) != 2 || nrow(matrix) < order + 1) {
    warning("Invalid argument: matrix must be a matrix with 2 columns and at least order + 1 rows\n")
    return(
      list(
        augcoeffmatrix = ifelse(is.matrix(matrix), matrix, NA),
        coefficients = NA,
        polynomial_string = NA,
        polynomial_function = NA
      )
    )
  }

  # Check if order is a positive integer but also
  # allow numeric type for order as long as it is an integer
  if (!is.numeric(order) || as.integer(order) != order || order < 1) {
    warning("Invalid argument: order must be a positive integer.\n")
    return(
      list(
        augcoeffmatrix = NA,
        coefficients = NA,
        polynomial_string = NA,
        polynomial_function = NA
      )
    )
  }

  # Get references to the x and y vectors (more readable)
  x_vector <- matrix[, 1]
  y_vector <- matrix[, 2]

  # Create an augmented coefficient matrix to solve for the coefficients
  augcoeffmatrix <- NULL
  for (i in 1:(order + 1)) {
    # Create a new row and compute for each column
    row <- NULL

    # Coefficient part of the matrix
    for (j in 1:(order + 1)) {
      row <- c(row, sum(x_vector^(i + j - 2)))
    }

    # Augmented part of the matrix (RHS)
    row <- c(row, sum(x_vector^(i - 1) * y_vector))

    # Bind the rows to a matrix
    augcoeffmatrix <- rbind(augcoeffmatrix, row)
  }

  # Add dimnames to the matrix
  colnames(augcoeffmatrix) <- c(paste0("x^", 0:order), "RHS")
  rownames(augcoeffmatrix) <- paste0("Order ", 0:order)

  # Solve the augmented coefficient matrix
  coefficients <- GaussJordanMethod(augcoeffmatrix)
  if (any(is.na(coefficients))) {
    warning("No solution found for the augmented coefficient matrix.\n")
    return(
      list(
        augcoeffmatrix = augcoeffmatrix,
        coefficients = NA,
        polynomial_string = NA,
        polynomial_function = NA
      )
    )
  }

  # Create another coefficient vector that is fully positive
  coefficient_magnitudes <- c()
  for (c in coefficients) {
    if (c < 0) {
      coefficient_magnitudes <- c(coefficient_magnitudes, -c)
    } else {
      coefficient_magnitudes <- c(coefficient_magnitudes, c)
    }
  }

  # Create the polynomial string then iterate over the coefficients/magnitudes
  # to add the terms to the polynomial string
  polynomial_string <- ""
  for (i in 0:order) {
    if (i == 0) { # Constant term, since x^0 = 1
      polynomial_string <- paste0(polynomial_string, coefficients[i + 1]) # We care about the sign of this term
    } else if (i == 1) { # Linear term, since x^1 = x
      polynomial_string <- paste0(polynomial_string, coefficient_magnitudes[i + 1], "*x")
    } else { # Polynomial terms, x^i
      polynomial_string <- paste0(polynomial_string, coefficient_magnitudes[i + 1], "*x^", i)
    }

    # Add a plus or negative sign to separate terms
    if (i != order) {
      if (coefficients[i + 2] >= 0) {
        polynomial_string <- paste0(polynomial_string, " + ")
      } else {
        polynomial_string <- paste0(polynomial_string, " - ")
      }
    }
  }

  # Create a function that evaluates the polynomial string
  polynomial_function <- eval(parse(text = paste("function(x)", polynomial_string)))

  # Prepend f(x) = to the polynomial string
  polynomial_string <- paste0("f(x) = ", polynomial_string)

  # Return the results of the polynomial regression
  list(
    augcoeffmatrix = augcoeffmatrix,
    coefficients = coefficients,
    polynomial_string = polynomial_string,
    polynomial_function = polynomial_function
  )
}

# A wrapper for the PolynomialRegression function that plots the regression line
# Produces a png file with the plot
RegressionPlotter <- function(order, matrix, plot_name = "Polynomial Regression") {
  # Get the output of the polynomial regression
  regression <- PolynomialRegression(order, matrix)

  # Check if the regression was successful
  if (any(is.na(regression$coefficients))) {
    warning("Regression failed. No plot generated.\n")
    return(regression)
  }

  # Get the directory of the current script
  current_dir <- dirname(sys.frame(1)$ofile)

  # Create a png file with the plot in the current directory
  png(file.path(current_dir, paste0(plot_name, ".png")), width = 1000, height = 1000, pointsize = 24)

  # Plot the points of the data
  plot(
    x = matrix[, 1],
    y = matrix[, 2],
    main = plot_name,
    xlab = colnames(matrix)[1],
    ylab = colnames(matrix)[2],
    pch = 19,
    col = "#508a9e"
  )

  # Get the polynomial function from the regression result
  func <- regression$polynomial_function

  # Plot the polynomial regression line/curve
  curve(
    func,
    add = TRUE,
    col = "#006800",
    lwd = 2
  )

  # Close the png file
  dev.off()

  # Return the result of the polynomial regression
  regression
}

# Quadratic regression taken from video
quadratic_matrix <- matrix(
  c(
    0, 1, 2, 3, 4, 5,
    2.1, 7.7, 13.6, 27.2, 40.9, 61.1
  ),
  ncol = 2,
  dimnames = list(
    rep("", 6),
    c("x", "y")
  )
)

# Plants data from exercise part 2
plants_matrix <- matrix(
  c(
    5, 10, 15, 20, 25, 30, 35, 40,
    2.1, 4.5, 9.3, 15.8, 24.1, 34.7, 48.3, 60.2
  ),
  ncol = 2,
  dimnames = list(
    1:8,
    c("Days after Planting", "Plant Height (cm)")
  )
)

# Unsolvable case due to repeated x values
unsolvable_matrix <- matrix(
  c(
    1, 1, 1, # repeated x values
    2.0, 4.0, 8.0 # different y values
  ),
  ncol = 2,
  dimnames = list(
    rep("", 3),
    c("x", "y")
  )
)

quadratic_result <- PolynomialRegression(2, quadratic_matrix)
ListPrint(quadratic_result)

plants_result_linear <- RegressionPlotter(1, plants_matrix, "Days after Planting vs. Plant Height (Linear)")
ListPrint(plants_result_linear)

plants_result_quadratic <- RegressionPlotter(2, plants_matrix, "Days after Planting vs. Plant Height (Quadratic)")
ListPrint(plants_result_quadratic)

plants_result_cubic <- RegressionPlotter(3, plants_matrix, "Days after Planting vs. Plant Height (Cubic)")
ListPrint(plants_result_cubic)

unsolvable_result <- RegressionPlotter(1, unsolvable_matrix)
ListPrint(unsolvable_result)
