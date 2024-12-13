## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS PLNfixedVsupfit
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' An R6 Class to represent a PLNfit in a standard, general framework, with fixed eigenvectors of Covariance matrix
#'
#' @param responses the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param covariates design matrix (called X in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param offsets offset matrix (called O in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which PLN is called.
#' @param weights an optional vector of observation weights to be used in the fitting process.
#' @param nullModel null model used for approximate R2 computations. Defaults to a GLM model with same design matrix but not latent variable.
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param control a list for controlling the optimization. See details.
#' @param config part of the \code{control} argument which configures the optimizer
#'
#' @rdname PLNfixedVsupfit
#' @importFrom R6 R6Class

PLNfixedVsupfit <- R6Class(
  classname = "PLNfixedVsupfit",
  inherit = PLNfit_fixedv,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public  = list(
    #' @description Initialize a [`PLNfixedVsupfit`] model
    initialize = function(grouping, responses, covariates, offsets, weights, formula, control) {
      private$grouping <- grouping
      super$initialize(responses, cbind(covariates, grouping), offsets, weights, formula, control)
      
    },
    #' @description Update a [`PLNfixedVsupfit`] object
    update = function(B=NA, Sigma=NA, Omega=NA, L=NA, M=NA, S=NA, Z=NA, A=NA, Ji=NA, R2=NA, monitoring=NA) {
      super$update(B = B, Sigma = Sigma, Omega = Omega, M = M, S = S, Z = Z, A = A, Ji = Ji, R2 = R2, monitoring = monitoring)
      if (!anyNA(L)) private$L <- L
    },
    
    #' @description Call to the NLopt or TORCH optimizer and update of the relevant fields
    optimize = function(grouping, responses, covariates, offsets, weights, config) {
      super$optimize(responses, cbind(covariates, grouping), offsets, weights, config)
    },
  
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## Prediction methods --------------------
  #' @description Predict group of new samples
  # @inheritParams predict.PLNLDAfit
  #' @param newdata A data frame in which to look for variables, offsets and counts  with which to predict.
  #' @param type The type of prediction required. The default are posterior probabilities for each group (in either unnormalized log-scale or natural probabilities, see "scale" for details), "response" is the group with maximal posterior probability and "scores" is the average score along each separation axis in the latent space, with weights equal to the posterior probabilities.
  #' @param scale The scale used for the posterior probability. Either log-scale ("log", default) or natural probabilities summing up to 1 ("prob").
  #' @param prior User-specified prior group probabilities in the new data. If NULL (default), prior probabilities are computed from the learning set.
  #' @param control a list for controlling the optimization. See [PLN()] for details.
  #' @param envir Environment in which the prediction is evaluated
  predict = function(newdata,
                     type = c("posterior", "group", "latent"),
                     scale = c("log", "prob"),
                     prior = NULL,
                     control = PLN_param(backend="nlopt"), envir = parent.frame()) {
    
    type <- match.arg(type)
    
    if (type == "scores") scale <- "prob"
    scale <- match.arg(scale)
    
    ## Extract the model matrices from the new data set with initial formula
    args <- extract_model(call = call("PLNfixedVsup", formula = private$formula, data = newdata), envir = envir)
    ## Remove intercept to prevent interference with binary coding of the grouping factor
    args$X <- args$X[ , colnames(args$X) != "(Intercept)", drop = FALSE]
    
    ## Problem dimensions
    n.new  <- nrow(args$Y)
    p      <- ncol(args$Y)
    groups <- levels(private$grouping)
    K <- length(groups)
 
    ## Initialize priors
    if (is.null(prior)) {
      prior <- table(private$grouping)
    } else {
      names(prior) <- groups
    }
    if (any(prior <= 0) || anyNA(prior)) stop("Prior group proportions should be positive.")
    prior <- prior / sum(prior)
    
    ## Compute conditional log-likelihoods of new data, using previously estimated parameters
    cond.log.lik <- matrix(0, n.new, K)
    ve_step_list <- NULL
    for (k in 1:K) { ## One VE-step to estimate the conditional (variational) likelihood of each group
      grouping <- factor(rep(groups[k], n.new), levels = groups)
      X <- cbind(args$X, grouping)
      ve_step_list[[k]] <- super$optimize_vestep(X, args$O, args$Y, args$w,
                                       B = self$model_par$B,
                                       L = self$model_par$L,
                                       control = control)
      cond.log.lik[, k] <- ve_step_list[[k]]$Ji
    }
    
    ## Compute (unnormalized) posterior probabilities
    log.prior <- rep(1, n.new) %o% log(prior)
    log.posterior <- cond.log.lik + log.prior
    
    res <- log.posterior
    ## trick to avoid rounding errors before exponentiation
    row_max <- apply(res, 1, max)
    res <- exp(sweep(res, 1, row_max, "-"))
    res <- sweep(res, 1, rowSums(res), "/")
    rownames(res) <- rownames(newdata)
    colnames(res) <- groups
    g <- apply(res, 1, which.max)-1

    if (type == "posterior") {
      return (res)
    }
    if (type == "group") {
      return (g)
    }
    
    if (type == "latent") {
      latent <- matrix(NA, nrow=n.new, ncol=p)
      for (i in 1:n.new){
        if (g[i]==0) {latent[i,] <- ve_step_list[[1]]$Z[i,]}
        else {latent[i,] <- ve_step_list[[2]]$Z[i,]}
      }
      return (latent)
    }

  },
  
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## Print methods -------------------------
  #' @description User friendly print method
  show = function() {
    super$show(paste0("Supervised Model \n"))
  }
),

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## PRIVATE MEMBERS
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
private = list(
  C        = NULL,
  P        = NULL,
  Mu       = NULL,
  grouping = NULL,
  svdLDA   = NULL
),

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  ACTIVE BINDINGS ----
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
active = list(
  #' @field rank the dimension of the current model
  rank = function() {nlevels(private$grouping) - 1},
  #' @field nb_param number of parameters in the current PLN model
  nb_param = function() {self$p * (self$d + self$rank)},
  #' @field model_par a list with the matrices associated with the estimated parameters of the PLN model: B (covariates), Sigma (latent covariance), C (latent loadings), P (latent position) and Mu (group means)
  model_par = function() {
    par <- super$model_par
    par$C  <- private$C
    par$P  <- private$P
    par$Mu <- private$Mu
    par
  },
  #' @field group_means a matrix of group mean vectors in the latent space.
  group_means = function() {
    self$model_par$Mu
  }
)

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  END OF THE CLASS ----
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)