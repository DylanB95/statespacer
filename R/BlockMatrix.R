#' Combine Matrices into a Block Diagonal Matrix
#'
#' Creates a block diagonal matrix with its arguments as the blocks.
#'
#' @param A The first block matrix to be put on the diagonal.
#' @param ... Other block matrices that should be put on the diagonal.
#'
#' @details
#' `BlockMatrix()` tries to coerce its arguments to a matrix,
#' using \code{\link[base:matrix]{as.matrix}}.
#'
#' @return Block diagonal matrix consisting of the specified matrices.
#'
#' @examples
#' BlockMatrix(diag(ceiling(50 * stats::runif(5))), matrix(1:8, 4, 2), c(14,8))
#'
#' @export
BlockMatrix <- function(A = NULL, ...) {

  # Retrieve the other arguments and put them in a list
  dots <- list(...)

  # Initialise first matrix
  new <- A
  if (!is.null(new) & !is.matrix(new)) {
    new <- as.matrix(new)
  }

  # If no other block matrices are specified, return the first matrix
  if (length(dots) == 0) {
    return(new)
  }

  # Add block matrices on the diagonal of the new matrix
  for (block in dots) {

    # If current block is NULL, go to the next block
    if (is.null(block)) {
      next
    }
    # If current block has a 0 dimension, go to the next block
    if (any(dim(block) == 0)) {
      next
    }
    if (!is.matrix(block)) {
      block <- as.matrix(block)
    }

    # Store the matrix before adding a block to it
    old <- new

    # Store the dimensions before adding the block
    dim_old <- dim(old)

    # If Null, dimensions equal 0
    if (is.null(dim_old)) {
      dim_old <- c(0, 0)
    }

    # How many rows and columns will be added?
    dim_block <- dim(block)

    # Dimensions after adding the block to the matrix
    dim_new <- dim_old + dim_block

    # Initialise the new matrix
    new <- matrix(0, dim_new[1], dim_new[2])

    # Store the old matrix into the new matrix, using the proper indices
    if (dim_old[1] > 0 & dim_old[2] > 0) {
      new[1:dim_old[1], 1:dim_old[2]] <- old
    }

    # Store the block into the new matrix, using the proper indices
    new[(dim_old[1] + 1) : dim_new[1], (dim_old[2] + 1) : dim_new[2]] <- block
  }

  # Return the final constructed matrix
  return(new)
}
