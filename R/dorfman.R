#' @name dorfman
#' @aliases dorfman
#' @title Dorfman Group Testing Decoder
#' @description Implements two-stage Dorfman group testing algorithm with sensitivity and specificity considerations.
#' @details The function performs:
#' \enumerate{
#'   \item Group testing: Samples are pooled into groups of size \code{c} and tested
#'   \item Individual testing: If a group tests positive, all members are tested individually
#' }
#' The probability models used:
#' \deqn{P(Z_g = 1|Y_g) = \begin{cases}
#' Se[1] & \text{if } \sum Y_g > 0 \\
#' 1-Sp[1] & \text{otherwise}
#' \end{cases}}
#' \deqn{P(Z_i = 1|Y_i) = \begin{cases}
#' Se[2] & \text{if } Y_i = 1 \\
#' 1-Sp[2] & \text{otherwise}
#' \end{cases}}
#' where \eqn{Z_g} is group test result, \eqn{Z_i} is individual test result, and \eqn{Y_g}, \eqn{Y_i} are true statuses.
#'
#' @param Y.true Numeric vector (n x 1) of true binary statuses (1=positive, 0=negative)
#' @param Se Numeric vector of length 2: (group sensitivity, individual sensitivity)
#' @param Sp Numeric vector of length 2: (group specificity, individual specificity)
#' @param c Integer group size (must divide length(Y.true))
#'
#' @return List with:
#' \describe{
#'   \item{Z}{Test results matrix (J x (c+3)) where columns are:
#'             [Test result, Group size, Test type (1=group, 2=individual), Sample indices]}
#'   \item{Y}{Individual results matrix (n x 4) where columns are:
#'             [Test result, Test type, Group test ID, Individual test ID]}
#' }
#'
#' @importFrom stats plogis rbinom
#'
#' @examples
#' set.seed(123)
#' Y.true <- rbinom(10, 1, 0.3)
#' result <- dorfman(Y.true, c(0.95, 0.98), c(0.98, 0.99), 2)
#' result
#' @export
dorfman <- function(Y.true, Se, Sp, c) {
  # Validate inputs
  if(length(Se) != 2 || length(Sp) != 2)
    stop("Se and Sp must be length 2 vectors")
  if(length(Y.true) %% c != 0)
    stop("Group size c must divide length(Y.true)")

  N <- length(Y.true)            # Total number of samples
  Jmax <- N + N/c                # Maximum possible tests (all groups positive)
  J <- 1                         # Test counter

  # Initialize output matrices
  Y <- matrix(-99, nrow = N, ncol = 4)  # Individual tracking matrix
  Z <- matrix(-99, nrow = Jmax, ncol = c + 3)  # Test results matrix

  # Process each group
  for(j in 1:(N/c)) {
    group_indices <- ((j-1)*c+1):(j*c)  # Indices for current group

    # --- Group Testing Stage ---
    group_pos <- sum(Y.true[group_indices]) > 0  # True group status
    prob <- ifelse(group_pos, Se[1], 1-Sp[1])   # Test probability
    Z[J, 1] <- rbinom(1, 1, prob)               # Simulate group test result
    Z[J, 2] <- c                                # Record group size
    Z[J, 3] <- 1                                # Mark as group test (type 1)
    Z[J, 4:(c+3)] <- group_indices              # Store sample indices

    # Initialize individual records
    Y[group_indices, 1] <- 0       # Temporary status
    Y[group_indices, 2] <- 1       # Test count (1 group test)
    Y[group_indices, 3] <- J       # Group test ID
    J <- J + 1                      # Increment test counter

    # --- Individual Testing Stage (if group positive) ---
    if(Z[J-1, 1] == 1) {
      for(k in group_indices) {
        ind_pos <- Y.true[k] > 0              # True individual status
        prob <- ifelse(ind_pos, Se[2], 1-Sp[2]) # Individual test probability
        Z[J, 1] <- rbinom(1, 1, prob)         # Simulate individual test
        Z[J, 2] <- 1                          # Group size = 1 (individual)
        Z[J, 3] <- 2                          # Mark as individual test (type 2)
        Z[J, 4] <- k                          # Sample index

        # Update individual record
        Y[k, 1] <- Z[J, 1]    # Individual test result
        Y[k, 2] <- 2          # Test count (group + individual)
        Y[k, 4] <- J          # Individual test ID
        J <- J + 1            # Increment test counter
      }
    }
  }

  # Trim unused test rows
  Z <- Z[1:(J-1), , drop = FALSE]

  return(list(Z = Z, Y = Y))
}
