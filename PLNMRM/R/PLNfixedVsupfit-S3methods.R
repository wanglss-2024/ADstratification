## =========================================================================================
##
## PUBLIC S3 METHODS FOR PLNLDAfit
##
## =========================================================================================

## Auxiliary functions to check the given class of an objet
isPLNfixedVsupfit <- function(Robject) {inherits(Robject, "PLNfixedVsupfit"       )}

#' Predict group of new samples
#'
#' @name predict.PLNfixedVsupfit
#'
#' @param object an R6 object with class [`PLNfixedVsupfit`]
#' @param newdata A data frame in which to look for variables, offsets and counts  with which to predict.
#' @param type The type of prediction required. The default are posterior probabilities for each group (in either unnormalized log-scale or natural probabilities, see "scale" for details), "response" is the group with maximal posterior probability and "scores" is the average score along each separation axis in the latent space, with weights equal to the posterior probabilities.
#' @param scale The scale used for the posterior probability. Either log-scale ("log", default) or natural probabilities summing up to 1 ("prob").
#' @param prior User-specified prior group probabilities in the new data. If NULL (default), prior probabilities are computed from the learning set.
#' @param control a list for controlling the optimization. See [PLN()] for details.
#' @param ... additional parameters for S3 compatibility. Not used
#' @return A matrix of posterior probabilities for each group (if type = "posterior"), a matrix of (average) scores in the latent space (if type = "scores") or a vector of predicted groups (if type = "response").
#' @export
predict.PLNfixedVsupfit <- function(object, newdata,
                              type = c("posterior", "group", "latent"),
                              scale = c("log", "prob"),
                              prior = NULL,
                              control = PLN_param(backend="nlopt"), ...) {
  stopifnot(isPLNfixedVsupfit(object))
  object$predict(newdata, type, scale, prior, control, parent.frame())
}

#' Extracts model coefficients from objects returned by [PLNfixedVsup()]
#'
#' @description The method for objects returned by [PLNfixedVsupfit()] only returns
#'              coefficients associated to the \deqn{\Theta} part of the model (see the PLNLDA vignette
#'              for mathematical details).
#'
#' @name coef.PLNfixedVsupfit
#'
#' @param object an R6 object with class PLNLDAfit
#' @param ... additional parameters for S3 compatibility. Not used
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPLNLDA <- PLNLDA(Abundance ~ Wind, grouping = Group, data = trichoptera)
#' coef(myPLNLDA)
#' @return Either NULL or a matrix of coefficients extracted from the PLNLDAfit model.
#'
#' @export
coef.PLNfixedVsupfit <- function(object, ...) {
  stopifnot(isPLNfixedVsupfit(object))
  n.covariates <- object$d - object$rank - 1
  if (n.covariates == 0) return(NULL)
  object$model_par$Theta[ , 1:n.covariates, drop = FALSE]
}

