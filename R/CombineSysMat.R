#' Combine the Z System Matrices of two State Space Components
#'
#' Combines the Z system matrices of two State Space components.
#'
#' @param Z1 Z system matrix of the first component.
#' @param Z2 Z system matrix of the second component.
#'
#' @return
#' Returns a matrix if both `Z1` and `Z2` are matrices.
#' Returns an array if at least one of them is an array.
#'
#' @noRd
CombineZ <- function(Z1, Z2) {

  # Check if Z1 and Z2 are matrices
  is_mat_Z1 <- is.matrix(Z1)
  is_mat_Z2 <- is.matrix(Z2)

  # Dimensions of Z matrices
  dimZ1 <- dim(Z1)
  dimZ2 <- dim(Z2)

  # Check if rows match
  if (dimZ1[1] != dimZ2[1]) {
    stop("Number of rows of Z matrices must match!")
  }

  # Compute result
  if (is_mat_Z1 & is_mat_Z2) {
    result <- cbind(Z1, Z2)
  } else if (is_mat_Z1 & !is_mat_Z2) {
    result <- array(
      apply(
        Z2, 3,
        function(x) {
          cbind(
            Z1,
            matrix(x, dimZ2[1], dimZ2[2])
          )
        }
      ),
      dim = c(dimZ1[1], sum(dimZ1[2], dimZ2[2]), dimZ2[3])
    )
  } else if (!is_mat_Z1 & is_mat_Z2) {
    result <- array(
      apply(
        Z1, 3,
        function(x) {
          cbind(
            matrix(x, dimZ1[1], dimZ1[2]),
            Z2
          )
        }
      ),
      dim = c(dimZ2[1], sum(dimZ2[2], dimZ1[2]), dimZ1[3])
    )
  } else {

    # Check if 3rd dimensions match
    if (dimZ1[3] != dimZ2[3]) {
      stop("3rd dimensions of Z matrices must match!")
    }

    result <- sapply(
      1:dimZ1[3],
      FUN = function(i) {
        cbind(
          matrix(Z1[,,i], dimZ1[1], dimZ1[2]),
          matrix(Z2[,,i], dimZ2[1], dimZ2[2])
        )
      },
      simplify = "array"
    )
  }

  return(result)
}

#' Combine non-Z System Matrices of two State Space Components
#'
#' Combines the T, R, or Q system matrices of two State Space components.
#'
#' @param S1 system matrix of the first component.
#' @param S2 system matrix of the second component.
#'
#' @return
#' Returns a matrix if both `S1` and `S2` are matrices.
#' Returns an array if at least one of them is an array.
#'
#' @noRd
CombineTRQ <- function(S1, S2) {

  # Check if S1 and S2 are matrices
  is_mat_S1 <- is.matrix(S1)
  is_mat_S2 <- is.matrix(S2)

  # Dimensions of S matrices
  dimS1 <- dim(S1)
  dimS2 <- dim(S2)

  # Compute result
  if (is_mat_S1 & is_mat_S2) {
    result <- BlockMatrix(S1, S2)
  } else if (is_mat_S1 & !is_mat_S2) {
    result <- array(
      apply(
        S2, 3,
        function(x) BlockMatrix(S1, as.matrix(x))
      ),
      dim = c(sum(dimS1[1], dimS2[1]), sum(dimS1[2], dimS2[2]), dimS2[3])
    )
  } else if (!is_mat_S1 & is_mat_S2) {
    result <- array(
      apply(
        S1, 3,
        function(x) BlockMatrix(as.matrix(x), S2)
      ),
      dim = c(sum(dimS1[1], dimS2[1]), sum(dimS1[2], dimS2[2]), dimS1[3])
    )
  } else {

    # Check if 3rd dimensions match
    if (dimS1[3] != dimS2[3]) {
      stop("3rd dimensions of S matrices must match!")
    }

    result <- sapply(
      1:dimS1[3],
      FUN = function(i) BlockMatrix(as.matrix(S1[,,i]), as.matrix(S2[,,i])),
      simplify = "array"
    )
  }

  return(result)
}
