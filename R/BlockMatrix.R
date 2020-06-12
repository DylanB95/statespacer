#' Combine Matrices into a Block Diagonal Matrix
#'
#' Creates a block diagonal matrix with its arguments as the blocks.
#'
#' @param ... Matrices that should be put on the diagonal.
#'
#' @details
#' `BlockMatrix()` tries to coerce its arguments to a matrix,
#' using \code{\link[base:matrix]{as.matrix}}.
#'
#' @return Block diagonal matrix having the specified matrices on its diagonal.
#'
#' @author Dylan Beijers, \email{dylanbeijers@@gmail.com}
#'
#' @examples
#' BlockMatrix(diag(ceiling(9 * stats::runif(5))), matrix(1:8, 4, 2), c(14, 8))
#' @export
BlockMatrix <- function(...) {

  # Retrieve the arguments and put them in a list
  dots <- list(...)

  # Eliminate NULLs
  dots[vapply(dots, is.null, logical(1))] <- NULL

  # If no block matrices are specified, return a 0 x 0 matrix
  if (length(dots) == 0) {
    return(matrix(nrow = 0, ncol = 0))
  }

  # Coerce arguments to matrices
  dots <- lapply(dots, as.matrix)

  # Dimensions of blocks
  dim_blocks <- vapply(dots, dim, integer(2))

  # Dimensions of result matrix
  dim_result <- rowSums(dim_blocks)

  # Initialise result matrix
  result <- matrix(0, nrow = dim_result[[1]], ncol = dim_result[[2]])

  # Indices of blocks in result matrix
  row_2 <- cumsum(dim_blocks[1, ])
  col_2 <- cumsum(dim_blocks[2, ])
  row_1 <- c(1, 1 + row_2)
  col_1 <- c(1, 1 + col_2)

  # Add block matrices on the diagonal of the result matrix
  for (i in seq_along(dots)) {

    # If current block has a 0 dimension, go to the next block
    if (any(dim_blocks[, i] == 0)) {
      next
    }

    # Indices of the block in the result
    row_ind <- row_1[[i]]:row_2[[i]]
    col_ind <- col_1[[i]]:col_2[[i]]

    # Add the block into the result matrix
    result[row_ind, col_ind] <- dots[[i]]
  }

  # Return the final constructed matrix
  return(result)
}
