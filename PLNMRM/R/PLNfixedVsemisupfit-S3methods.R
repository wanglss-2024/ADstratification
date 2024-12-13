## =========================================================================================
##
## PUBLIC S3 METHODS FOR PLNfixedVsemisupfit
##
## =========================================================================================

## Auxiliary functions to check the given class of an objet
isPLNfixedVsemisupfit <- function(Robject) {inherits(Robject, "PLNfixedVsemisupfit"       )}



#' Prediction for a [`PLNfixedVsemisupfit`] object
#'
#' Predict either posterior probabilities for each group or latent positions based on new samples
#'
#' @param object an R6 object with class [`PLNfixedVsemisupfit`]
#' @param newdata A data frame in which to look for variables, offsets and counts with which to predict.
#' @param type The type of prediction required. The default `posterior` are posterior probabilities for each group ,
#'  `response` is the group with maximal posterior probability and `latent` is the averaged latent in the latent space,
#'  with weights equal to the posterior probabilities.
#' @param prior User-specified prior group probabilities in the new data. The default uses a uniform prior.
#' @param control a list-like structure for controlling the fit. See [PLN_param()] for details.
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A matrix of posterior probabilities for each group (if type = "posterior"), a matrix of (average) position in the
#' latent space (if type = "position") or a vector of predicted groups (if type = "response").
#' @export
predict.PLNfixedVsemisupfit <-
  function(object, newdata,
           type = c("posterior", "group", "latent"),
           prior = matrix(rep(1/object$k, object$k), nrow(newdata), object$k, byrow = TRUE),
           control = PLN_param(), ...) {
    
    stopifnot(isPLNfixedVsemisupfit(object))
    object$predict(newdata, type, prior, control, parent.frame())
    
  }

#' Extract model coefficients
#'
#' @description Extracts model coefficients from objects returned by [PLN()] and its variants
#'
#' @name coef.PLNfixedVsemisupfit
#'
#' @param object an R6 object with class [`PLNfixedVsemisupfit`]
#' @param type type of parameter that should be extracted. Either "main" (default) for \deqn{\Theta},
#' "means" for \deqn{\mu}, "mixture" for \deqn{\pi} or "covariance" for \deqn{\Sigma}
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A matrix of coefficients extracted from the PLNfit model.
#'
#' @export
coef.PLNfixedVsemisupfit <- function(object, type = c("main", "means", "covariance", "mixture"), ...) {
  stopifnot(isPLNfixedVsemisupfit(object))
  switch(match.arg(type),
         main       = object$model_par$Theta,
         means      = object$group_means,
         mixture    = object$mixtureParam,
         covariance = object$model_par$Sigma)
}

