## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS PLNfixedVunsupfit
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' An R6 Class to represent a PLNfit in a mixture framework
#'
#' @description The function \code{\link{PLNfixedVunsup}} produces a collection of models which are instances of object with class \code{PLNmixturefit}.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for ...
#'
#' @param responses the matrix of responses common to every models
#' @param covariates the matrix of covariates common to every models
#' @param offsets the matrix of offsets common to every models
#' @param weights an optional vector of observation weights to be used in the fitting process.
#' @param control a list for controlling the optimization.
#' @param clusters the dimensions of the successively fitted models
#' @param formula model formula used for fitting, extracted from the formula in the upper-level call
#' @param cluster the number of clusters of the current model
#' @param nullModel null model used for approximate R2 computations. Defaults to a GLM model with same design matrix but not latent variable.
#'
#' @include PLNfit-class.R
#'
#' @importFrom R6 R6Class
#' @importFrom purrr map2_dbl map2
#' @importFrom purrr map_int map_dbl
#' @seealso The function \code{\link{PLNfixedVunsup}}
PLNfixedVunsupfit <-
  R6Class(classname = "PLNfixedVunsupfit",
          inherit = PLNfit_fixedv,
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## PRIVATE MEMBERS
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    private = list(
      tau        = NA, # posterior probabilities of cluster belonging
      monitoring = NA, # a list with optimization monitoring quantities
      A0         = NA    , 
      A1         = NA    ,
      Z0         = NA    , 
      Z1         = NA    ,
      Ji0         = NA    , # element-wise approximated loglikelihood
      Ji1         = NA     # element-wise approximated loglikelihood
    ),
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## PUBLIC MEMBERS
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    public  = list(
      #' @description Initialize a [`PLNmixturefit`] model
      #'@param posteriorProb matrix of posterior probability for cluster belonging
      initialize = function(responses, covariates, offsets, weights, formula, control) {
        private$formula <- formula
        n <- nrow(responses)
        p <- ncol(responses)
        private$q <- control$rank
        private$tau   <- control$config_optim$tau0
        private$L <- control$config_optim$L0
        private$B <- matrix(0, ncol(covariates)+1, private$q)
        private$V <- control$V[,1:private$q]
        private$M  <- matrix(0, n, private$q)
        private$S  <- matrix(1, n, private$q)
        
        private$optimizer$main <- nlopt_optimize_fixedv_unsup
        private$optimizer$vestep <- nlopt_optimize_vestep_fixedv_unsup
        
      },
      
      #' @description Update a [`PLNmixturefit`] object
      update = function(B=NA, Sigma=NA, Omega=NA, L=NA, M=NA, S=NA, Z0=NA, Z1=NA, A0=NA, A1=NA, Ji0=NA, Ji1=NA, R2=NA, monitoring=NA) {
        private$B <- B
        private$Sigma <- Sigma
        private$Omega <- Omega
        private$L <- L
        private$M <- M
        private$S <- S
        private$Z0 <- Z0
        private$Z1 <- Z1
        private$A0 <- A0
        private$A1 <- A1
        private$Ji0 <- Ji0
        private$Ji1 <- Ji1
        private$monitoring <- monitoring
      },
      
      #' @description Optimize a [`PLNfixedVunsupfit`] model
      #' @param config a list for controlling the optimization
      optimize = function(responses, covariates, offsets, weights, config) {
        
        X0 <- cbind(covariates, 0)
        X1 <- cbind(covariates, 1)
        
        ## ===========================================
        ## INITIALISATION
        cond <- FALSE; iter <- 1
        objective   <- numeric(config$maxit_out); objective[iter]   <- Inf
        convergence <- numeric(config$maxit_out); convergence[iter] <- NA
      
        ## ===========================================
        ## OPTIMISATION
        while (!cond) {
          cat("\n # iter")
          cat(iter)
          iter <- iter + 1
          if (config$trace > 1) cat("", iter)

          ## ---------------------------------------------------
          ## M - STEP
          ## UPDATE parameters
          
          args <- list(data   = list(Y = responses, X0 = X0, X1 = X1, O = offsets, w = weights, w0 = private$tau[,1], w1 = private$tau[,2]),
                       params = list(B = private$B, M = private$M, S = private$S, L = private$L, V = private$V),
                       config = config)
          optim_out <- do.call(private$optimizer$main, args)
          do.call(self$update, optim_out)

          ## ---------------------------------------------------
          ## E - STEP
          ## UPDATE THE POSTERIOR PROBABILITIES
          if (self$k > 1) { # only needed when at least 2 components!
            private$tau <-
              #sapply(private$comp, function(comp) comp$loglik_vec) %>% # Jik
              cbind(private$Ji0, private$Ji1) %>%
              sweep(2, log(self$mixtureParam), "+") %>% # computation in log space
              apply(1, .softmax) %>%        # exponentiation + normalization with soft-max
              t() %>% .check_boundaries()   # bound away probabilities from 0/1
          }

          ## Assess convergence
          objective[iter]   <- -self$loglik
          convergence[iter] <- abs(objective[iter-1] - objective[iter]) /abs(objective[iter])
          if ((convergence[iter] < config$ftol_out) | (iter >= config$maxit_out)) cond <- TRUE
          # if ((iter >= config$maxit_out)) cond <- TRUE
        }

        ## ===========================================
        ## OUTPUT
        ## formatting parameters for output
        private$monitoring = list(objective        = objective[2:iter],
                                  convergence      = convergence[2:iter],
                                  outer_iterations = iter-1)
      },

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Prediction methods --------------------
      #' @description Predict group of new samples
      #' @param newdata A data frame in which to look for variables, offsets and counts with which to predict.
      #' @param type The type of prediction required. The default `posterior` are posterior probabilities for each group ,
      #'  `response` is the group with maximal posterior probability and `latent` is the averaged latent coordinate (without
      #'  offset and nor covariate effects),
      #'  with weights equal to the posterior probabilities.
      #' @param prior User-specified prior group probabilities in the new data. The default uses a uniform prior.
      #' @param control a list-like structure for controlling the fit. See [PLNmixture_param()] for details.
      #' @param envir Environment in which the prediction is evaluated
      predict = function(newdata,
                         type = c("posterior", "group", "latent"),
                         prior = matrix(rep(1/self$k, self$k), nrow(newdata), self$k, byrow = TRUE),
                         control = PLN_param(), envir = parent.frame()) {

        type  <- match.arg(type)

        ## Extract the model matrices from the new data set with initial formula
        args <- extract_model(call("PLNfixedVunsup", formula = private$formula, data = newdata), envir)
        ## Remove intercept to prevent interference with binary coding of the grouping factor
        args$X <- args$X[ , colnames(args$X) != "(Intercept)", drop = FALSE]
        
        n_new <- nrow(args$Y)

        ## Sanity checks
        stopifnot(ncol(args$X) == self$d-1)
        stopifnot(ncol(args$O) == self$p, ncol(args$Y) == self$p)
        stopifnot(nrow(prior)  == nrow(args$X), ncol(prior) == self$k)
        stopifnot(nrow(args$X) == n_new)

        ## Initialize the posteriorProbabilities
        tau <- prior

        ## ===========================================
        ## INITIALISATION
        cond <- FALSE; iter <- 1
        objective   <- numeric(control$config_optim$maxit_out); objective[iter]   <- Inf
        convergence <- numeric(control$config_optim$maxit_out); convergence[iter] <- NA

        ## ===========================================
        ## OPTIMISATION
        while(!cond) {
          iter <- iter + 1
          if (control$trace > 1) cat("", iter)
          
          ## ---------------------------------------------------
          ## VE step
          args_optim <- list(data   = list(Y = args$Y, X0 = cbind(args$X,0), X1 = cbind(args$X,1), O = args$O,
                                     w = args$w, w0 = tau[,1], w1 = tau[,2]),
                       params = list(M = matrix(0, n_new, private$q), S = matrix(1, n_new, private$q), V = private$V),
                       B = as.matrix(private$B),
                       L = as.matrix(private$L),
                       config = control$config_optim)
          optim_out <- do.call(private$optimizer$vestep, args_optim)

          ## E - STEP
          ## UPDATE THE POSTERIOR PROBABILITIES
          if (self$k > 1) { # only needed when at least 2 components!
            tau <-
              cbind(optim_out$Ji0, optim_out$Ji1) %>%
              sweep(2, log(colMeans(tau)), "+") %>% # computation in log space
              apply(1, .softmax) %>%        # exponentiation + normalization with soft-max
              t() %>% .check_boundaries()   # bound away probabilities from 0/1
          }

          ## Assess convergence
          J_ik <- cbind(optim_out$Ji0, optim_out$Ji1)
          J_ik[tau <= .Machine$double.eps] <- 0
          rowSums(tau * J_ik) - rowSums(.xlogx(tau)) + tau %*% log(colMeans(tau))
          objective[iter]   <- -sum(J_ik)
          convergence[iter] <- abs(objective[iter-1] - objective[iter]) /abs(objective[iter])
          if ((convergence[iter] < control$config_optim$ftol_out) | (iter >= control$config_optim$maxit_out)) cond <- TRUE
          # if ((iter >= control$config_optim$maxit_out)) cond <- TRUE
        }
        
        g <- apply(tau, 1, which.max)

        if (type == "latent") {
          mat <- matrix(NA, nrow=n_new, ncol=self$p)
          for (i in 1:n_new) {
            if (g[i]==1){mat[i,] <- optim_out$Z0[i,]}
            else {mat[i,] <- optim_out$Z1[i,]}
          }
        }
        
        switch(type,
          "posterior" = {rownames(tau) <- rownames(newdata); tau},
          "group"  = g,
          "latent" = mat
          )
      },

      ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ## Print methods ---------------------
      #' @description User friendly print method
      show = function() {
        cat("Unsupervised Model \n")
      },
      #' @description User friendly print method
      print = function() self$show()
    ),
    active = list(
      #' @field n number of samples
      n = function() {nrow(private$tau)},
      #' @field p number of dimensions of the Z
      p = function() {ncol(self$model_par$Sigma)},
      #' @field q number of dimensions of the low rank latent space
      q = function() {ncol(self$model_par$B)},
      #' @field k number of components
      k = function() {ncol(private$tau)},
      #' @field d number of covariates
      d = function() {nrow(self$model_par$B)},
      #' @field latent a matrix: values of the latent vector (Z in the model)
      latent = function(){
        mat = matrix(NA, nrow=self$n, ncol=self$q)
        for (i in 1:self$n) {
          if (self$memberships[i]==1){mat[i,] <- private$Z0[i,]}
          else {mat[i,] <- private$Z1[i,]}
        }
        mat
      },
      #' @field posteriorProb matrix ofposterior probability for cluster belonging
      posteriorProb = function(value) {if (missing(value)) return(private$tau) else private$tau <- value},
      #' @field memberships vector for cluster index
      memberships = function(value) {apply(private$tau, 1, which.max)},
      #' @field mixtureParam vector of cluster proportions
      mixtureParam  = function() {colMeans(private$tau)},
      #' @field optim_par a list with parameters useful for monitoring the optimization
      optim_par  = function() {private$monitoring},
      #' @field nb_param number of parameters in the current PLN model
      nb_param  = function() {(self$k-1) + self$q * self$d + sum(map_int(self$components, 'nb_param'))},
      #' @field loglik variational lower bound of the loglikelihood
      loglik = function() {sum(self$loglik_vec)},
      #' @field loglik_vec element-wise variational lower bound of the loglikelihood
      loglik_vec = function() {
        J_ik <- cbind(private$Ji0, private$Ji1)
        J_ik[private$tau <= .Machine$double.eps] <- 0
        rowSums(private$tau * J_ik) - rowSums(.xlogx(private$tau)) + private$tau %*% log(self$mixtureParam)
        },
      #' @field model_par a list with the matrices of the model parameters: B (covariates), Sigma (covariance), Omega (precision matrix), plus some others depending on the variant)
      model_par  = function() {list(B = private$B, Sigma = private$Sigma, Omega = private$Omega, Theta = t(private$B), L = private$L, V = private$V,Pi = self$mixtureParam)},
      #' @field var_par a list with the matrices of the variational parameters: M (means) and S2 (variances)
      var_par    = function() {list(M = private$M, S2 = private$S**2, S = private$S)},
      #' @field vcov_model character: the model used for the covariance (either "spherical", "diagonal" or "full")
      vcov_model = function() {private$covariance},
      #' @field fitted a matrix: fitted values of the observations (A in the model)
      fitted = function() {private$mix_up('fitted')},
      #' @field group_means a matrix of group mean vectors in the latent space.
      group_means = function() {
        self$components %>%
          map(function(C) C$model_par$Theta)  %>%
          setNames(paste0("group_", 1:self$k)) %>% as.data.frame()
      }
    )
)
