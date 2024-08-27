#' A Tree-guided Compositional Variable Selection method (TCVS)
#'
#' This function allows you to
#' identify outcome-associated components (OTUs) within high-dimensional microbial compositional data with
#' hierarchical taxonomic tree information inherent among microbial taxa.
#'
#' @param X the compositional covariate matrix.
#' @param Z the CLR transforamtion matrix (design matrix).
#' @param y the response vector.
#' @param P the taxonomic structure of OTUs.
#' @param method the method used for selecting the tuning parameter. The default is "BIC".
#' @param maxlam the maximum lambda value for the tuning parameter lambda candidates.
#' @param minlam the minimum lambda value for the tuning parameter lambda candidates.
#' @param nlam the number of lambda value candidates.
#' @param fdr the nominal FDR level. The default is 0.05.
#' @param seed an integer value used to set the seed of the random number generator for ensuring reproducibility.
#' The default is 2023.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item TCVS.beta.hat - The estimated coefficients from the TCVS method.
#'   \item S.TCVS - The selected variables from the TCVS method.
#' }
#'
#' @examples
#' # Example usage
#' \dontrun{
#' n = 200
#' p = 60
#' sim.setting = 1
#' fdr.normial = 0.05
#' method = "BIC"
#' loop_start_index = 1
#' seed = 2023
#' type = "Dirmult" # type to generate compositional data
#' normalizeMethod = "Rowsum"
#' beta.slack.factor = NULL
#' pseudocount = 0.5
#' maxlam = 0.1
#' minlam = 1e-7
#' nlam = 2
#' X <- get.all.OTU(n, p, trans = 0, type, normalizeMethod, pseudocount, seed = seed)
#' # clr transformation
#' Z <- get.all.OTU(n, p, trans = 2, type, normalizeMethod, pseudocount, seed = seed)
#' y <- get.all.y(Z, setting = sim.setting, beta.slack.factor, seed = seed)
#' result.TCVS <- TCVS(
#'   X = X,
#'   Z = Z,
#'   y = y,
#'   P = P,
#'   method = method,
#'   maxlam = maxlam,
#'   minlam = minlam,
#'   nlam = nlam,
#'   fdr,
#'   seed = seed
#' )
#' }
#' @export
#' @importFrom knockoff knockoff.threshold
TCVS <- function(X,Z,y,
                 P,
                 method = "BIC",
                 maxlam,
                 minlam,
                 nlam,
                 fdr = 0.05,
                 seed = 2023) {
   set.seed(seed)
   Result <- list()
   lambdas.TCVS <- vector(mode = "numeric")
   n = dim(X)[1]
   p = dim(X)[2]

   current_date <- format(Sys.Date(), format = "%m%d")

   A0 <- matrix(c(rep(1, p), rep(0, p), rep(0, p), rep(1, p)), nrow = 2, ncol = 2 * p, byrow = TRUE)
   Z.tilde <- get.all.knockoff(Z, seed) # generate knockoff copies
   Z.comb <- list()
   cvxrResult.TCVS <- list()
   lambdas.TCVS.all.results <- list()

   W_TCVS <- matrix(nrow = 1, ncol = p)
   TCVS.beta.hat <- matrix(nrow = 1, ncol = 2 * p)
   threshold.TCVS <- vector(mode = "numeric")
   S.TCVS <- matrix(0, nrow = 1, ncol = p)
   S.TCVS.name <- paste(current_date,
    "n", n,
    "p", p,
    "maxlam", maxlam,
    "minlam", minlam,
    "nlam", nlam,
    "S.TCVS.txt",
    sep = "_"
   )
    Z.comb <- cbind(Z, Z.tilde)

    y.vector <- unlist(y)
    y.mean <- mean(y.vector)
    y.vector <- y.vector - y.mean

    Z.comb.matrix <- matrix(unlist(Z.comb), nrow = nrow(Z.comb))
    Z.comb.mean <- colMeans(Z.comb.matrix)
    Z.comb.matrix <-
      scale(Z.comb.matrix, center = Z.comb.mean, scale = F) # centering before scaling

    lambdas.TCVS.all.results<- BIC_TCVS(
      Z = Z.comb.matrix,
      Y = y.vector,
      p,
      group = P,
      maxlam,
      minlam,
      nlam,
      lambda_seq = NULL,
      A0 = A0,
      CoefNormalization = T
    )

    lambdas.TCVS <- lambdas.TCVS.all.results$best_lambda_BIC # select optimal lambda by BIC
    bic_values <- sapply(lambdas.TCVS.all.results[[1]], function(x) x$BIC)

    min_bic_index <- which.min(bic_values)
    cvxrResult.TCVS <- lambdas.TCVS.all.results[[1]][[min_bic_index]]$betahat
    TCVS.beta.hat[1, ] <- as.numeric(unlist(cvxrResult.TCVS))
    TCVS.beta.hat[1, which(abs(TCVS.beta.hat[1, ]) < 1e-08)] <- 0
    tmp2 <- TCVS.beta.hat[1, ]
    for (j in 1:p) {
      W_TCVS[1, j] <- abs(tmp2[j]) - abs(tmp2[j + p])
    }

    threshold.TCVS[1] <- knockoff.threshold(W = W_TCVS[1, ], fdr = fdr, offset = 0)
    S.TCVS[1, which(W_TCVS[1, ] >= threshold.TCVS[1])] <- 1

    # save variables
    #Result$Z.tilde <- Z.tilde
    Result$TCVS.beta.hat <- TCVS.beta.hat
    Result$S.TCVS <- S.TCVS
    return(Result)
}

#' Generate compositional matrix X and CLR transformation matrix Z
#'
#' This function generates the X and Z matrices based on the specified type and parameters.
#'
#' @param n the number of samples for the simulated compositional data.
#' @param p the number of features the for the simulated compositional data.
#' @param type the type of the simulated data to be generated. The options are "Lognormal_previous" and "Dirmult".
#' @param normalizeMethod the method of normalization when type = "Dirmult". The default is "Rowsum".
#' @param pseudocount  a pseudo-count value added to the count matrix to prevent zero values in the read count.
#'
#' @return A list containing the following matrices:
#' \itemize{
#'   \item X - The compositional matrix.
#'   \item log_X - The logarithm of the compositional matrix.
#'   \item Z - The CLR transformation matrix.
#' }
#'
#' @examples
#' # Example usage for Dirmult type
#' result <- TCVS_Generate_X_Z(n = 200, p = 60, type="Lognormal_previous",
#' normalizeMethod, pseudocount = 0.5)
#'
#' @export
#' @importFrom cluster pam
#' @importFrom dirmult rdirichlet
TCVS_Generate_X_Z <- function(
    n,
    p,
    type = c("Lognormal_previous", "Dirmult"),
    normalizeMethod = "Rowsum",
    pseudocount) {
  if (type == "Lognormal_previous") {
    mu <- c(rep(log(p / 2), 5), rep(0, p - 5))
    rho <- 0.5
    Sigma <- toeplitz(rho^(0:(p - 1)))

    W <- matrix(rnorm(n * p), n) %*% chol(Sigma)
    X <- matrix(nrow = n, ncol = p)
    log_X <- matrix(nrow = n, ncol = p)
    Z <- matrix(nrow = n, ncol = p)
    X <- exp(W) / rowSums(exp(W))
    log_X <- log(X)
    Z <- log_X - rowMeans(log_X)
  } else if (type == "Dirmult") {
    m <- 50
    k <- 4 ## 4 time pts within each subject
    n <- m * k ## total sample size
    data(throat.tree)
    data(throat.otu.tab)
    data(DirMultOutput)
    nClus <- 20
    depth <- 10000
    tree <- throat.tree # tree = midpoint(tree)
    tree.dist <- cophenetic(tree)
    obj <- pam(tree.dist, nClus)
    clustering <- obj$clustering
    otu.ids <- tree$tip.label
    p.est <- dd$pi
    names(p.est) <- names(dd$pi)
    theta <- dd$theta
    gplus <- (1 - theta) / theta
    p.est <- p.est[otu.ids]
    g.est <- p.est * gplus
    p.clus <- sort(tapply(p.est, clustering, sum), decreasing = T)
    comm <-
      matrix(0, n, length(g.est)) ## OTU matrix n \times micro_dim: all 20 clusters
    rownames(comm) <- 1:nrow(comm)
    colnames(comm) <- names(g.est)
    comm.p <- comm ## comm.p hold the underlying proportions
    dim(comm.p) ## 200 samples, 856 OTUs

    nSeq <- rnbinom(n, mu = depth, size = 25)
    for (i in 1:n) {
      comm.p[i, ] <- rdirichlet(1, g.est)[1, ]
      comm[i, ] <- rmultinom(1, nSeq[i], prob = comm.p[i, ])[, 1]
    }
    comm <- comm + pseudocount
    comm <- comm[, 1:p]
    if (normalizeMethod == "Rowsum") {
      X <- comm / rowSums(comm) # normalization by Rowsum
      log_X <- log(X)
      Z <- log_X - rowMeans(log_X)
    } else {
      stop("normalization method must be Rowsum")
    }
  } else {
    stop("type must be either Lognormal_previous or dirmult")
  }
  return(list(X, log_X, Z))
}

#' Generate OTU Matrix
#'
#' This function generates the OTU matrix.
#'
#' @param n the number of samples for the simulated compositional data.
#' @param p the number of features the for the simulated compositional data.
#' @param trans an integer specifying the type of transformation applied to the data.
#' The options are 0 (no transformation), 1 (log transformation), and 2 (clr transformation).
#' @param type the type of the simulated data to be generated. The options are "Lognormal_previous" and "Dirmult".
#' @param normalizeMethod the normalization method.
#' @param pseudocount the pseudocount value.
#' @param seed an integer value used to set the seed of the random number generator for ensuring reproducibility.
#' The default is 2023.
#' @return The generated OTU matrix.
#' @export
get.all.OTU <- function(
    n,
    p,
    trans = c(0, 1, 2),
    type,
    normalizeMethod,
    pseudocount,
    seed = 2023) {
  OTU <- NULL
  set.seed(seed)
    if (trans == 0) {
      OTU <- TCVS_Generate_X_Z(n, p,
        type,
        normalizeMethod = "Rowsum",
        pseudocount
      )[[1]]
    } else if (trans == 1) {
      OTU <- TCVS_Generate_X_Z(n,
        p,
        type,
        normalizeMethod = "Rowsum",
        pseudocount
      )[[2]]
    } else if (trans == 2) {
      OTU <- TCVS_Generate_X_Z(n,
        p,
        type,
        normalizeMethod = "Rowsum",
        pseudocount
      )[[3]]
    }
  OTU
  return(OTU) # matrix
}

#' Generate Response Variable Based on the CLR matrix Z
#'
#' This function generates a response variable `y` as a linear combination of predictors
#' specified in matrix `Z`, influenced by a vector of coefficients `beta`. The coefficients
#' are determined based on predefined settings and optionally modified by a slack factor.
#' Random noise is added to the response to simulate more realistic scenarios.
#'
#' @param Z the CLR transforamtion matrix (design matrix).
#' @param setting a numeric vector indicating the setting for coefficient generation.
#'                Setting 1 generates a fixed set of coefficients with specific non-zero values.
#'                Setting 2 generates coefficients by sampling from a uniform distribution.
#' @param beta.slack.factor an optional numeric factor to adjust the coefficients by a given factor.
#'                          If not provided, coefficients are used as defined by the `setting`.
#'                          The default is 1.
#' @param seed an integer value to set the seed for random number generation to ensure reproducibility.
#'
#' @return A numeric vector representing the generated response variable `y`.
#'
#' @examples
#'   Z <- matrix(rnorm(100 * 20), ncol = 20)
#'   y <- get.all.y(Z, setting = 1, beta.slack.factor, seed)
#' @export
get.all.y <- function(Z, setting = c(1, 2), beta.slack.factor = 1, seed = 2023) {
  n <- dim(Z)[1]
  p <- dim(Z)[2]

  beta <- rep(0, p)
  if (is.null(beta.slack.factor)) {
    if (setting == 1) {
      beta[1] <- 1
      beta[2] <- -1
      beta[4] <- 0.8
      beta[5] <- -0.8
      beta[8] <- -1.5
      beta[9] <- -0.5
      beta[10] <- 2
      beta[11] <- 1.2
      beta[12] <- -1.2
      beta[18] <- 0.7
      beta[19] <- 0.8
      beta[20] <- -1.5
    }
    if (setting == 2) {
      set.seed(2023)
      y <- list()
      case3 <- sort(sample(seq(1, 60), size = 12))
      case <- case3

      set.seed(2023)
      beta_1 <- sample(seq(-2, 2, by = 0.1), size = 12, replace = TRUE)
      beta <- rep(0, p)

      beta[case3] <- beta_1
      beta[case3[12]] <- (-sum(beta[case3[1:11]]))
    }
  } else {
    if (setting == 1) {
      beta[1] <- 1
      beta[2] <- -1
      beta[4] <- 0.8
      beta[5] <- -0.8
      beta[8] <- -1.5
      beta[9] <- -0.5
      beta[10] <- 2
      beta[11] <- 1.2
      beta[12] <- -1.2
      beta[18] <- 0.7
      beta[19] <- 0.8
      beta[20] <- -1.5
      beta <- beta.slack.factor * beta
    }
    if (setting == 2) {
      set.seed(2023)
      case3 <- sort(sample(seq(1, 60), size = 12))
      case <- case3

      set.seed(2023)
      beta_1 <- sample(seq(-2, 2, by = 0.1), size = 12, replace = TRUE)
      beta <- rep(0, p)

      beta[case3] <- beta_1
      beta[case3[12]] <- (-sum(beta[case3[1:11]]))
      beta <- beta.slack.factor * beta
    }
  }
    set.seed(seed = seed + 1)
    y<- as.numeric(Z %*% beta) + rnorm(n)
  return(y)
}

#' Generate Knockoff Copy Matrix
#'
#' This function generates the knockoff copy matrix.
#'
#' @param Z the original matrix.
#' @param seed the seed for random number generation. The default is 2023.
#'
#' @return The knockoff copy matrix.
#' @export
#' @importFrom knockoff create.second_order
get.all.knockoff <- function(Z, seed = 2023) {
    set.seed(seed + 2)
    KF <- create.second_order(Z, method = "sdp")
  return(KF)
}

#' Group Lasso Penalty with L2 Norm (Without Square)
#'
#' This function computes the penalty term in the augmented problem incorporating auxiliary knockoff copies and tree structure with
#' L2 norm (without square) for a given set of coefficients, groups, and a regularization parameter.
#'
#' @param beta  a numeric vector of coefficients.
#' @param group the taxonomic structure of OTUs.
#' @param lambda the regularization parameter.
#' @param p the number of features the for the simulated compositional data.
#' @param CoefNormalization Logical. Whether to normalize the coefficients within each group. The default is T.
#'
#' @return The computed penalty.
#' @export
#' @importFrom CVXR p_norm
#' @importFrom CVXR multiply
Group.lasso <- function(beta, group, lambda, p, CoefNormalization = T) {
  N <- nrow(group)
  if (N > (2 * p)) {
    group_start <- 2 * p + 1

    if (CoefNormalization == T) {
      add_terms <- apply(group[group_start:N, ], 1, function(row) {
        p_norm(multiply(row, beta), 2) / length(which(row == 1)) # calculate the elementwise product of the inputs.
      })
    } else if (CoefNormalization == F) {
      add_terms <- apply(group[group_start:N, ], 1, function(row) {
        p_norm(multiply(row, beta), 2) # calculate the elementwise product of the inputs.
      })
    } else {
      stop("CoefNormalization  must be T or F")
    }
    total_sum_part2 <- Reduce(`+`, add_terms)
  } else {
    total_sum_part2 <- 0
  }
  total_sum_part1 <- p_norm(beta, 1) # L1 norm
  total_penalty <- lambda * (total_sum_part1 + total_sum_part2)
  return(total_penalty)
}

#' CVXR Optimization Result
#'
#' This function performs an optimization using the CVXR package to solve an augmented regression problem with the penalty.
#'
#' @param y the response vector.
#' @param Z the design matrix matrix.
#' @param A0 The constraint matrix, which enforces the condition that the coefficients of the original features and their
#' knockoff counterparts sum to zero, respectively.
#' @param lambda the regularization parameter.
#' @param p the number of features the for the simulated compositional data.
#' @param group the taxonomic structure of OTUs.
#' @param CoefNormalization Logical. Whether to normalize the coefficients within each group. The default is T.
#'
#' @return The estimated coefficients from the optimization.
#' @export
#' @importFrom CVXR Variable
#' @importFrom CVXR sum_squares
#' @importFrom CVXR Problem
#' @importFrom CVXR Minimize
#' @importFrom CVXR solve
cvxr.result <- function(y, Z, A0, lambda, p, group, CoefNormalization = T) {
  beta <- Variable(2 * p) # optimization variable
  constr <- list(matrix(A0[1, ], nrow = 1) %*% beta == 0, matrix(A0[2, ], nrow = 1) %*% beta == 0)
  n <- length(y)
  loss <- sum_squares(y - Z %*% beta) / (2 * n)
  obj <- loss + Group.lasso(beta, group, lambda = lambda, p,  CoefNormalization = T)
  prob <- Problem(Minimize(obj), constr)
  result <- solve(prob)
  return(betahat = result$getValue(beta))
}

#' BIC Calculation for CVXR Optimization
#'
#' This function calculates the Bayesian Information Criterion (BIC) for a given set of parameters using the CVXR package.
#'
#' @param Z the design matrix matrix.
#' @param Y the response vector.
#' @param p the number of features the for the simulated data.
#' @param group the taxonomic structure of OTUs.
#' @param maxlam the maximum lambda value for the tuning parameter lambda candidates.
#' @param minlam the minimum lambda value for the tuning parameter lambda candidates.
#' @param nlam the number of lambda value candidates.
#' @param lambda_seq the sequence of lambda values. If NULL, a sequence will be generated.
#' @param A0 the constraint matrix.
#' @param CoefNormalization Logical. Whether to normalize the coefficients within each group. The default is T.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item total.results - A list of results for each lambda value.
#'   \item best_lambda_BIC - The lambda value with the minimum BIC.
#' }
#' @export
BIC_TCVS <- function(Z, Y, p, group, maxlam, minlam, nlam, lambda_seq, A0, CoefNormalization = T) {
  N <- nrow(group)
  n <- length(Y)
  m <- length(which(group != 0))
  lambda_results <- list()
  if (is.null(lambda_seq)) {
    if (nlam > 1) {
      # 2*max(abs(crossprod(Z.centered,Y.centered)/n))
      lambda_seq <- exp(seq(from = log(maxlam), to = log(minlam), length.out = nlam))
    } else {
      lambda_seq <- 0.1
    }
  }
  lambda_results <- lapply(lambda_seq, function(lambda) {
    result <- cvxr.result(Y, Z, A0, lambda, p, group, CoefNormalization)
    result[which(abs(result[, ]) < 1e-08), ] <- 0
    k <- sum(result != 0)
    loss <- sum((Y - Z %*% result)^2) / n
    BIC_values <- log(loss) + k * log(n) / n
    return(list(
      betahat = result,
      lambda = lambda,
      BIC_values = BIC_values
    ))
  })
  best_lambda_BIC <- lambda_seq[which.min(sapply(lambda_results, function(x) x$BIC_values))]
  return(list(
    total.results <- lambda_results,
    best_lambda_BIC = best_lambda_BIC
  ))
}

## simulation of multi-normal
#' @importFrom MASS mvrnorm
simunorm=function(Sigma,n,p){
  out=mvrnorm(n, rep(0,p), Sigma, tol = 1e-6)
  return(out)
}

### transfer distance matrix to kernel matrix
D2K = function(D){
  n = NROW(D)
  centerM = diag(n ) - 1/n
  K = -0.5*centerM %*% (D*D) %*% centerM
  eK= eigen(K, symmetric = T)
  K = eK$vector %*% diag(abs(eK$values)) %*% t(eK$vector)
  return(K)
}

