#' An R6 Class to represent a PLNfit in a standard, general framework
#'
#' @description The function [PLN()] fit a model which is an instance of a object with class [`PLNfit`].
#'
#' Fields are accessed via active binding and cannot be changed by the user.
#'
## Parameters common to all PLN-xx-fit methods (shared with PLNfit but inheritance does not work)
#' @param responses the matrix of responses (called Y in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param covariates design matrix (called X in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param offsets offset matrix (called O in the model). Will usually be extracted from the corresponding field in PLNfamily-class
#' @param weights an optional vector of observation weights to be used in the fitting process.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which PLN is called.
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param control a list-like structure for controlling the fit, see [PLN_param()].
#' @param config part of the \code{control} argument which configures the optimizer
#' @param nullModel null model used for approximate R2 computations. Defaults to a GLM model with same design matrix but not latent variable.
#' @param B matrix of regression matrix
#' @param Sigma variance-covariance matrix of the latent variables
#' @param Omega precision matrix of the latent variables. Inverse of Sigma.
#'
#' @inherit PLN details
#'
#' @rdname PLNfit
#' @include PLNfit-class.R
#' @importFrom R6 R6Class
#'
#' @examples
#' \dontrun{
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLN <- PLN(Abundance ~ 1, data = trichoptera)
#' class(myPLN)
#' print(myPLN)
#' }
PLNfit <- R6Class(
  classname = "PLNfit",
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    ## PRIVATE INTERNAL FIELDS
    formula    = NA    , # the formula call for the model as specified by the user
    B          = NA    , # regression parameters of the latent layer
    Sigma      = NA    , # covariance matrix of the latent layer
    Omega      = NA    , # precision matrix of the latent layer. Inverse of Sigma
    S          = NA    , # variational parameters for the variances
    M          = NA    , # variational parameters for the means
    Z          = NA    , # matrix of latent variable
    A          = NA    , # matrix of expected counts (under variational approximation)
    Ji         = NA    , # element-wise approximated loglikelihood
    R2         = NA    , # approximated goodness of fit criterion
    optimizer  = list(), # list of links to the functions doing the optimization
    monitoring = list(), # list with optimization monitoring quantities
    B.new      = NA    , # regression parameters of the latent layer
    Sigma.new  = NA    , # covariance matrix of the latent layer
    Omega.new  = NA    , # precision matrix of the latent layer. Inverse of Sigma
    S.new      = NA    , # variational parameters for the variances
    M.new      = NA    , # variational parameters for the means
    Z.new      = NA    , # matrix of latent variable
    A.new      = NA    , # matrix of expected counts (under variational approximation)
    Ji.new     = NA    , # element-wise approximated loglikelihood

  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## CONSTRUCTOR
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    #' @description Initialize a [`PLNfit`] model
    #' @importFrom stats lm.wfit lm.fit poisson residuals coefficients runif
    initialize = function(responses, covariates, offsets, weights, formula, control) {
      ## problem dimensions
      n <- nrow(responses); p <- ncol(responses); d <- ncol(covariates)
      ## set up various quantities
      private$formula <- formula # user formula call
      ## initialize the variational parameters
      if (isPLNfit(control$inception)) {
        if (control$trace > 1) cat("\n User defined inceptive PLN model")
        stopifnot(isTRUE(all.equal(dim(control$inception$model_par$B), c(d,p))))
        stopifnot(isTRUE(all.equal(dim(control$inception$var_par$M)  , c(n,p))))
        private$Sigma <- control$inception$model_par$Sigma
        private$B     <- control$inception$model_par$B
        private$M     <- control$inception$var_par$M
        private$S     <- control$inception$var_par$S
      } else {
        if (control$trace > 1) cat("\n Use LM after log transformation to define the inceptive model")
        fits <- lm.fit(weights * covariates, weights * log((1 + responses)/(1 + exp(offsets))))
        private$B <- matrix(fits$coefficients, d, p)
        private$M <- matrix(fits$residuals, n, p)
        private$S <- matrix(1, n, p)
      }
      private$optimizer$main   <- ifelse(control$backend == "nlopt", nlopt_optimize, private$torch_optimize)
      private$optimizer$vestep <- nlopt_optimize_vestep
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## SETTER METHOD
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    #' @description
    #' Update a [`PLNfit`] object
    #' @param M     matrix of variational parameters for the mean
    #' @param S     matrix of variational parameters for the variance
    #' @param Ji    vector of variational lower bounds of the log-likelihoods (one value per sample)
    #' @param R2    approximate R^2 goodness-of-fit criterion
    #' @param Z     matrix of latent vectors (includes covariates and offset effects)
    #' @param A     matrix of fitted values
    #' @param monitoring a list with optimization monitoring quantities
    #' @return Update the current [`PLNfit`] object
    update = function(B=NA, Sigma=NA, Omega=NA, M=NA, S=NA, Ji=NA, R2=NA, Z=NA, A=NA, monitoring=NA) {
      if (!anyNA(B))      private$B  <- B
      if (!anyNA(Sigma))      private$Sigma  <- Sigma
      if (!anyNA(Omega))      private$Omega  <- Omega
      if (!anyNA(M))          private$M      <- M
      if (!anyNA(S))          private$S      <- S
      if (!anyNA(Z))          private$Z      <- Z
      if (!anyNA(A))          private$A      <- A
      if (!anyNA(Ji))         private$Ji     <- Ji
      if (!anyNA(R2))         private$R2     <- R2
      if (!anyNA(monitoring)) private$monitoring <- monitoring
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## GENERIC OPTIMIZER
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    #' @description Call to the NLopt or TORCH optimizer and update of the relevant fields
    optimize = function(responses, covariates, offsets, weights, config) {
      args <- list(data   = list(Y = responses, X = covariates, O = offsets, w = weights),
                   params = list(B = private$B, M = private$M, S = private$S),
                   config = config)
      optim_out <- do.call(private$optimizer$main, args)
      do.call(self$update, optim_out)
    },

    #' @description Result of one call to the VE step of the optimization procedure: optimal variational parameters (M, S) and corresponding log likelihood values for fixed model parameters (Sigma, B). Intended to position new data in the latent space.
    #' @param B Optional fixed value of the regression parameters
    #' @param Sigma variance-covariance matrix of the latent variables
    #' @return A list with three components:
    #'  * the matrix `M` of variational means,
    #'  * the matrix `S2` of variational variances
    #'  * the vector `log.lik` of (variational) log-likelihood of each new observation
    optimize_vestep = function(covariates, offsets, responses, weights,
                      B = self$model_par$B,
                      Omega = self$model_par$Omega,
                      control = PLN_param(backend = "nlopt")) {
      n <- nrow(responses); p <- ncol(responses)
      args <- list(data   = list(Y = responses, X = covariates, O = offsets, w = weights),
                   ## Initialize the variational parameters with the new dimension of the data
                   params = list(M = matrix(0, n, p), S = matrix(1, n, p)),
                   B = as.matrix(B),
                   Omega = as.matrix(Omega),
                   config = control$config_optim)
      optim_out <- do.call(private$optimizer$vestep, args)
      optim_out
    },

    #' @description Predict position, scores or observations of new data.
    #' @param newdata A data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
    #' @param type Scale used for the prediction. Either `link` (default, predicted positions in the latent space) or `response` (predicted counts).
    #' @param envir Environment in which the prediction is evaluated
    #' @return A matrix with predictions scores or counts.
    predict = function(newdata, type = c("link", "response"), envir = parent.frame()) {

      ## Extract the model matrices from the new data set with initial formula
      X <- model.matrix(formula(private$formula)[-2], newdata, xlev = attr(private$formula, "xlevels"))
      O <- model.offset(model.frame(formula(private$formula)[-2], newdata))

      ## mean latent positions in the parameter space
      EZ <- X %*% private$B
      if (!is.null(O)) EZ <- EZ + O
      EZ <- sweep(EZ, 2, .5 * diag(self$model_par$Sigma), "+")
      colnames(EZ) <- colnames(private$Sigma)

      type <- match.arg(type)
      results <- switch(type, link = EZ, response = exp(EZ))
      attr(results, "type") <- type
      results
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Print functions -----------------------
    #' @description User friendly print method
    #' @param model First line of the print output
    show = function(model = paste("A multivariate Poisson Lognormal fit with", self$vcov_model, "covariance model.\n")) {
      cat(model)
      cat("==================================================================\n")
      print(as.data.frame(round(self$criteria, digits = 3), row.names = ""))
      cat("==================================================================\n")
      cat("* Useful fields\n")
      cat("    $model_par, $latent, $latent_pos, $var_par, $optim_par\n")
      cat("    $loglik, $BIC, $ICL, $loglik_vec, $nb_param, $criteria\n")
      cat("* Useful S3 methods\n")
      cat("    print(), coef(), sigma(), vcov(), fitted()\n")
      cat("    predict(), predict_cond(), standard_error()\n")
    },

    #' @description User friendly print method
    print = function() { self$show() }

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Other functions ----------------
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field n number of samples
    n = function() {nrow(private$M)},
    #' @field q number of dimensions of the latent space
    q = function() {ncol(private$M)},
    #' @field p number of species
    p = function() {ncol(private$B)},
    #' @field d number of covariates
    d = function() {nrow(private$B)},
    #' @field nb_param number of parameters in the current PLN model
    nb_param   = function() {as.integer(self$p * self$d + self$p * (self$p + 1)/2)},
    #' @field model_par a list with the matrices of the model parameters: B (covariates), Sigma (covariance), Omega (precision matrix), plus some others depending on the variant)
    model_par  = function() {list(B = private$B, Sigma = private$Sigma, Omega = private$Omega, Theta = t(private$B), L = private$L, V = private$V,
                                  B.new = private$B.new, Sigma.new = private$Sigma.new, Omega.new = private$Omega.new, Theta.new = t(private$B.new))},
    #' @field var_par a list with the matrices of the variational parameters: M (means) and S2 (variances)
    var_par    = function() {list(M = private$M, S2 = private$S**2, S = private$S,
                                  M.new = private$M.new, S2.new = private$S.new**2, S.new = private$S.new)},
    #' @field optim_par a list with parameters useful for monitoring the optimization
    optim_par  = function() {c(private$monitoring, backend = private$backend)},
    #' @field latent a matrix: values of the latent vector (Z in the model)
    latent     = function() {private$Z},
    #' @field latent_pos a matrix: values of the latent position vector (Z) without covariates effects or offset
    latent_pos = function() {private$M},
    #' @field fitted a matrix: fitted values of the observations (A in the model)
    fitted     = function() {private$A},
    #' @field vcov_coef matrix of sandwich estimator of the variance-covariance of B (need fixed -ie known- covariance at the moment)
    vcov_coef = function() {attr(private$B, "vcov_variational")},
    #' @field vcov_model character: the model used for the residual covariance
    vcov_model = function() {"full"},
    #' @field weights observational weights
    weights     = function() {as.numeric(attr(private$Ji, "weights"))},
    #' @field loglik (weighted) variational lower bound of the loglikelihood
    loglik     = function() {sum(self$weights[self$weights > .Machine$double.eps] * private$Ji[self$weights > .Machine$double.eps]) },
    #' @field loglik_vec element-wise variational lower bound of the loglikelihood
    loglik_vec = function() {private$Ji},
    #' @field BIC variational lower bound of the BIC
    BIC        = function() {self$loglik - .5 * log(self$n) * self$nb_param},
    #' @field entropy Entropy of the variational distribution
    entropy    = function() {.5 * (self$n * self$q * log(2*pi*exp(1)) + sum(log(self$var_par$S2)))},
    #' @field ICL variational lower bound of the ICL
    ICL        = function() {self$BIC - self$entropy},
    #' @field R_squared approximated goodness-of-fit criterion
    R_squared  = function() {private$R2},
    #' @field criteria a vector with loglik, BIC, ICL and number of parameters
    criteria   = function() {data.frame(nb_param = self$nb_param, loglik = self$loglik, BIC = self$BIC, ICL = self$ICL)}
  )
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  END OF THE CLASS PLNfit
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS PLNfit_fixedv
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
#' @rdname PLNfit_fixedv
#' @importFrom R6 R6Class
#'
#' }
PLNfit_fixedv <- R6Class(
  classname = "PLNfit_fixedv",
  inherit = PLNfit,
  private = list(V = NULL, L = NULL, q = NULL),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public  = list(
    #' @description Initialize a [`PLNfit_fixedv`] model
    initialize = function(responses, covariates, offsets, weights, formula, control) {
      super$initialize(responses, covariates, offsets, weights, formula, control)
      private$optimizer$main <- ifelse(control$backend == "nlopt", nlopt_optimize_fixedv, private$torch_optimize)
      private$optimizer$vestep <- nlopt_optimize_vestep_fixedv
      private$optimizer$mstep <- nlopt_optimize_mstep_fixedv
      private$q <- ifelse(is.null(control$rank), p, control$rank)
      private$V <- control$V[,1:private$q]
      private$L <- control$config_optim$L0

      if (private$q < self$p) {
        private$M  <- matrix(0, self$n, private$q)
        private$S  <- matrix(1, self$n, private$q)
        private$B  <- matrix(0, ncol(covariates), private$q)
      }
      
    },
    #' @description Update a [`PLNfit_fixedv`] object
    update = function(B=NA, Sigma=NA, Omega=NA, L=NA, M=NA, S=NA, Z=NA, A=NA, Ji=NA, R2=NA, monitoring=NA) {
      super$update(B = B, Sigma = Sigma, Omega = Omega, M = M, S = S, Z = Z, A = A, Ji = Ji, R2 = R2, monitoring = monitoring)
      if (!anyNA(L)) private$L <- L
    },
    #' @description Call to the NLopt or TORCH optimizer and update of the relevant fields
    optimize = function(responses, covariates, offsets, weights, config) {
      args <- list(data   = list(Y = responses, X = covariates, O = offsets, w = weights),
                   params = list(B = private$B, M = private$M, S = private$S, L = private$L, V = private$V),
                   config = config)
      optim_out <- do.call(private$optimizer$main, args)
      do.call(self$update, optim_out)
    },
    
    #' @description Result of one call to the M step of teh optimization procedure: optimal model parameters (C,B) for fixed variational parameters (M,S)
    optimize_mstep = function(responses, covariates, offsets, weights, config) {
      # problem dimension
      n <- nrow(responses); p <- ncol(responses); q <- private$q
      
      # turn offset vector to offset matrix
      offsets <- as.matrix(offsets)
      if (ncol(offsets) == 1) offsets <- matrix(offsets, nrow = n, ncol = p)
      
      # get smoothed M and S
      #self$smooth(Y = responses, M = private$M, S = private$S)
      # Neural network for M
      model.m <- keras_model_sequential()
      model.m %>% 
        layer_dense(units = 30, activation = "relu", input_shape = c(p)) %>%
        layer_dense(units = 30, activation = "relu") %>%
        layer_dense(units = q, activation = "linear")
      model.m %>% compile(loss = "mean_squared_error", optimizer = optimizer_adam(lr = 0.001))
      history.m <- model.m %>% fit(
        log(responses+1), private$M,
        epochs = 150,
        batch_size = 32,
        validation_split = 0.3,
        callbacks = list(callback_early_stopping(monitor = "val_loss", 
                                                 patience = 10, 
                                                 restore_best_weights = TRUE))
      )
      private$M.new <- predict(model.m, log(responses+1))
      
      # Neural network for S
      model.s <- keras_model_sequential()
      model.s %>%
        layer_dense(units = 30, activation = "relu", input_shape = c(p)) %>%
        layer_dense(units = 30, activation = "relu") %>%
        layer_dense(units = q, activation = "softplus")
      model.s %>% compile(loss = "mean_squared_error", optimizer = optimizer_adam(lr = 0.001))
      history.s <- model.s %>% fit(
        log(responses+1), abs(private$S),
        epochs = 150,
        batch_size = 32,
        validation_split = 0.3,
        callbacks = list(callback_early_stopping(monitor = "val_loss",
                                                 patience = 10,
                                                 restore_best_weights = TRUE))
      )
      private$S.new <- predict(model.s, log(responses+1))
      
      args <- list(data   = list(Y = responses, X = covariates, O = offsets, w = weights),
                   ## Initialize the model parameters with the current values
                   params = list(B = private$B, L = private$L, V = private$V),
                   M = private$M.new,
                   S = private$S.new,
                   config = config)
      optim_out <- do.call(private$optimizer$mstep, args)
      private$B.new <- optim_out$B
      private$Sigma.new <- optim_out$Sigma
      
    },
    
    #' @description Result of one call to the VE step of the optimization procedure: optimal variational parameters (M, S) and corresponding log likelihood values for fixed model parameters (C, B). Intended to position new data in the latent space for further use with PCA.
    #' @return A list with three components:
    #'  * the matrix `M` of variational means,
    #'  * the matrix `S2` of variational variances
    #'  * the vector `log.lik` of (variational) log-likelihood of each new observation
    #' @description Result of one call to the VE step of the optimization procedure: optimal variational parameters (M, S) and corresponding log likelihood values for fixed model parameters (Sigma, B). Intended to position new data in the latent space.
    #' @param B Optional fixed value of the regression parameters
    #' @param L diagonal matrix with sqrt(Lambda)
    #' @return A list with three components:
    #'  * the matrix `M` of variational means,
    #'  * the matrix `S2` of variational variances
    #'  * the vector `log.lik` of (variational) log-likelihood of each new observation
    optimize_vestep = function(covariates, offsets, responses, weights,
                               B = self$model_par$B,
                               L = self$model_par$L,
                               control = PLN_param(backend = "nlopt")) {
      n <- nrow(responses); q <- private$q
      args <- list(data   = list(Y = responses, X = covariates, O = offsets, w = weights),
                   ## Initialize the variational parameters with the new dimension of the data
                   params = list(M = matrix(0, n, q), S = matrix(1, n, q), V = private$V),
                   B = as.matrix(B),
                   L = as.matrix(L),
                   config = control$config_optim)
      optim_out <- do.call(private$optimizer$vestep, args)
      optim_out
    },
    
    #' @description Update R2, fisher and std_err fields after optimization
    #' @param config a list for controlling the post-treatments (optional bootstrap, jackknife, R2, etc.). See details
    #' @details The list of parameters `config` controls the post-treatment processing, with the following entries:
    #' * trace integer for verbosity. should be > 1 to see output in post-treatments
    #' * jackknife boolean indicating whether jackknife should be performed to evaluate bias and variance of the model parameters. Default is FALSE.
    #' * bootstrap integer indicating the number of bootstrap resamples generated to evaluate the variance of the model parameters. Default is 0 (inactivated).
    #' * variational_var boolean indicating whether variational Fisher information matrix should be computed to estimate the variance of the model parameters (highly underestimated). Default is FALSE.
    #' * rsquared boolean indicating whether approximation of R2 based on deviance should be computed. Default is TRUE
    postTreatment = function(responses, covariates, offsets, weights = rep(1, nrow(responses)), config, nullModel = NULL) {
      super$postTreatment(responses, covariates, offsets, weights, config, nullModel)
      ## 6. compute and store matrix of standard variances for B with sandwich correction approximation
      if (config$sandwich_var) {
        if(config$trace > 1) cat("\n\tComputing sandwich estimator of the variance...")
        private$vcov_sandwich_B(responses, covariates)
      }
    }
  ),

  
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## END OF TORCH METHODS FOR OPTIMIZATION
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  active = list(
    #' @field nb_param number of parameters in the current PLN model
    nb_param   = function() {as.integer(self$p * self$d + self$p)},
    #' @field vcov_model character: the model used for the residual covariance
    vcov_model = function() {"fixedV"},
    #' @field vcov_coef matrix of sandwich estimator of the variance-covariance of B (needs known covariance at the moment)
    vcov_coef = function() {attr(private$B, "vcov_sandwich")},
    #' @field model_par a list with the matrices associated with the estimated parameters of the pPCA model: B (covariates), Sigma (covariance), Omega (precision) and C (loadings)
    model_par = function() {
      par <- super$model_par
      par$L <- private$L
      par
    }
  )
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  END OF THE CLASS PLNfit_fixedv
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
)

