# Tabamo, Euan Jed S. - CMSC 150 AB-3L
# Exercise 1 - Introduction to R Programming
# Description: Creates and prints a matrix that lists the frequencies of objects in a vector
# Date Created: August 28, 2024

# Counts the number of occurrences of an element in a vector
countElement <- function(vector, elementToCount) {
  # Initialize the count
  count <- 0L
  
  # Increment the count if we encounter the element in the vector
  for (element in vector) {
    if (element == elementToCount) {
      count <- count + 1L
    }
  }
  
  # Return the count number
  count
}

# Counts the number of times every element in setVector occurs in the originalVector
# Returns an integer vector with the number of occurrences of each element in setVector
countElements <- function(originalVector, setVector) {
  # Initialize a counts vector
  counts <-  c()
  
  # Append the number of times the element appears in the originalVector
  # To the counts vector
  for (element in setVector) {
    counts <- c(counts, countElement(originalVector, element))
  }
  
  # Return the counts vector
  counts
}

# Searches a vector for an element.
# Returns TRUE if the vector contains the element and FALSE if not.
vectorContains <- function(vector, elementToSearch) {
  # Return TRUE if we can find the element in the vector
  for (element in vector) {
    if (elementToSearch == element) {
      return(TRUE)
    }
  }
  # Return FALSE if we cannot find the element in the vector
  FALSE
}

# Removes duplicate items from a set, creating a new set vector.
# Returns the new set vector.
createSetVector <-  function(vector) {
  # Initialize a setVector
  setVector <- c()
  
  # Add all elements to the setVector if they are not yet in the setVector
  for (element in vector) {
    if (!vectorContains(setVector, element)) {
      setVector <- c(setVector, element)
    }
  }
  
  # Return the new setVector
  setVector
}

# Creates a character vector from a vector
# Returns the character vector
normalizeVector <- function(vector) {
  # Initialize the normalizedVector
  normalizedVector <- c()
  
  # Add every element to the normalizedVector as a character
  for (element in vector) {
    normalizedVector <- c(normalizedVector, as.character(element))
  }
  
  # Return the normalizedVector
  normalizedVector
}

# Creates a two-column matrix from two vectors.
# Returns the matrix with vectors as columns.
createMatrix <- function(setVector, countVector) {
  # Create normalized vectors for setVector and countVector
  normalizedSetVector <- normalizeVector(setVector)
  normalizedCountVector <- normalizeVector(countVector)
  
  # Get the length of the column to use in the matrix dimnames
  columnLength <-  length(setVector)
  
  # Return the matrix
  newMatrix <- matrix(
    c(normalizedSetVector, normalizedCountVector),
    nrow <- columnLength,
    ncol <-  2,
    byrow <-  FALSE,
    dimnames <- list(1:columnLength, c("Characters", "Frequency"))
  )
}

# Creates a frequency matrix from an input character vector
# Returns the frequency matrix
createFrequencyMatrix <- function(inputVector) {
  # Get the set of unique elements in a vector
  setVector <- createSetVector(inputVector)
  # Count the occurrences of each vector in the setVector
  countVector <- countElements(inputVector, setVector)
  # Create a matrix of the occurrences of each character
  finalMatrix <- createMatrix(setVector, countVector)
  
  # Return the matrix
  finalMatrix
}

# Exercise Sample
sample <- c("a", "b", "a", "c", "b", "a", "d")
sampleMatrix <- createFrequencyMatrix(sample)
print(sampleMatrix)

# Test Cases
test1 <-  c("E",
            "u",
            "a",
            "n",
            " ",
            "J",
            "e",
            "d",
            " ",
            "S",
            ".",
            " ",
            "T",
            "a",
            "b",
            "a",
            "m",
            "o")
matrix1 <- createFrequencyMatrix(test1)
print(matrix1)

test2 <- c("H", "a", "l", "l", "o", "W", "o", "r", "l", "d")
matrix2 <- createFrequencyMatrix(test2)
print(matrix2)

test3 <- c("I", "c", "h", " ", "h", "e", "i", "ß", "e", " ", "E", "u", "a", "n")
matrix3 <- createFrequencyMatrix(test3)
print(matrix3)
