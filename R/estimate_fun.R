#' Initial Estimation with sl3
#'
#' This function relies on the stacked ensemble learner in order to estimate relevant
#' parts of the likelihood as guided by the efficient influence curve.
#'
#' @param Y data.frame object containing the outcome variable.
#' @param X data.frame object containing the relevant covariates.
#' @param folds user-specified list of folds- it should correspond to an element of \code{origami}.
#' @param library list of \code{sl3} algorithms to be used for estimation.
#'
#' @return An object of class \code{tmle3opttx}.
#' \describe{
#' \item{cvRisk}{CV Risk for each algorithm and Super Learner.}
#' \item{valY}{Prediction based on the estimated SL fit.}
#' \item{library}{list of \code{sl3} algorithms to be used for estimation.}
#' \item{folds}{user-specified list of folds- it should correspond to an element of \code{origami}.}
#' \item{fullFit}{Overall SuperLearner fit.}
#' \item{cvFit}{Fit for each fold separately, for all folds.}
#' }
#'
#' @importFrom sl3 make_sl3_Task
#'
#' @export
#

estInit <- function(Y, X, folds=NULL, library) {

  covars<-names(X)
  outcome<-names(Y)
  all<-cbind.data.frame(Y,X)

  #Create sl3 task:
  task <- sl3::make_sl3_Task(all, covariates = covars, outcome = outcome, folds=folds)

  #Return sl3 trained object:
  sl3_fit<-sl3.fit(task=task, library=library)

  out <- list(cvRisk = sl3_fit$risk, valY = sl3_fit$pred, SL.library = library, folds = folds,
              fullFit = sl3_fit$sl.fit, cvFit=sl3_fit$cv.fit)

  class(out) <- c("initEst", "sl3")

  return(out)

}

############################################################################
# Create fold-specific predictions for QAW, Q1W, Q0W, B, Rule, K
############################################################################

#' Wrapper for split-specific predictions
#'
#' Function to make split-specific predictions more streamlined; it executes \code{cv_split}
#' over all the folds.
#'
#' @param folds user-specified list of folds. It should correspond to an element of \code{origami}.
#' @param Q data.frame containg the relevant outcome and covariates for estimating Q.
#' @param g data.frame containg the relevant outcome and covariates for estimating g.
#' @param estQ Q result of \code{initEst} format.
#' @param estg g result of \code{initEst} format.
#' @param nodes list contating names for the outcome, treatment, covariates and V-covariates as
#' observed in the provided data.
#'
#' @return An object of class \code{tmle3opttx}.
#' \describe{
#' \item{estSplit}{\cite{cv_split} results for each fold-specific fit. Contains fold-specific fit
#' predictions on all the samples.}
#' \item{valSplit}{Prediction for each sample based on the fold-specific fit when the said
#' sample was in the validation set. Contains prediction for QAW, Q1W, Q0W, B, Rule, K,
#' sample id and the fold for based on which the fit was generated.}
#' }
#'
#' @export
#

estSplit<-function(folds, Q, g, estQ, estg, nodes){

  #Observed data:
  Y<-Q[, nodes$Ynode, drop=FALSE]
  A<-Q[, nodes$Anode, drop=FALSE]
  W<-Q[, nodes$Wnodes, drop=FALSE]

  #Predict QAW, Q1W, Q0W, blip and rule based on CV-fold fits:
  estSplit <- origami::cross_validate(cv_split, folds, Y, A, W, estQ, estg, nodes, .combine = F)
  estSplit$errors<-NULL

  #Return predictions for the validation samples:
  estSplit_val<-extract_vals(folds, estSplit)

  estSplit$folds<-folds

  return(list(estSplit=estSplit,valSplit=estSplit_val))
}

#' Split-specific predictions
#'
#' Function to estimate split-specific predictions for QAW (E(Y|A,W)), Q1W (E(Y|A=1,W)),
#' Q0W (E(Y|A=0,W)), B, Rule based on maximizing B, and K.
#'
#' @param fold one of the folds from the folds object.
#' @param Y data.frame object for outcome.
#' @param A data.frame object for treatment.
#' @param W data.frame object for covaraites.
#' @param estQ object of class \code{initEst} with \code{sl3} results for Q.
#' @param estg object of class \code{initEst} with \code{sl3} results for g.
#' @param nodes list of node names corresponding to the outcome, treatment, covariate and V nodes.
#'
#' @return An object of class \code{tmle3opttx}.
#' \describe{
#' \item{Y}{Observed outcome (Y).}
#' \item{A}{Observed exposure (A).}
#' \item{QAW}{Split-specific prediction of Q (E(Y|A,W)) on all data.}
#' \item{Q1W}{Split-specific prediction of Q(1,W) (E(Y|A=1,W)) on all data.}
#' \item{Q0W}{Split-specific prediction of Q(0,W) (E(Y|A=0,W)) on all data.}
#' \item{pA1}{Split-specific prediction of g (P(A=1|W)) on all data.}
#' \item{B}{Split-specific prediction of (A/pA1 - (1 - A)/(1 - pA1)) * (Y - QAW) + Q1W - Q0W.}
#' \item{Rule}{Split-specific prediction of the maximizing rule.}
#' \item{K}{Split-specific absolute value of B.}
#' }
#'
#' @importFrom sl3 make_sl3_Task
#'
#' @export
#'

cv_split <- function(fold, Y,A,W, estQ, estg, nodes){

  Q=cbind.data.frame(Y,A,W)
  g=cbind.data.frame(A,W)

  #Get fold index
  v <- origami::fold_index()

  #Get the coefficients:
  coefQ<-estQ$fullFit$coefficients
  coefg<-estg$fullFit$coefficients

  #Get split-specific fits:
  splitQ <- estQ$cvFit[[v]]
  splitg <- estg$cvFit[[v]]

  #Predict on all the data:
  new_folds<-origami::make_folds(Q, fold_fun = origami::folds_resubstitution)[[1]]

  #QAW:
  covars<-c(names(A),names(W))
  outcome<-names(Y)

  task<-sl3::make_sl3_Task(Q, covariates = covars, outcome = outcome,folds=new_folds)
  QAW <- as.matrix(splitQ$predict(task)) %*% as.matrix(coefQ)

  #Q1W:
  Q1<-Q
  Q1[,"A"] = 1

  task<-sl3::make_sl3_Task(Q1, covariates = covars, outcome = outcome, folds=new_folds)
  Q1W <- as.matrix(splitQ$predict(task)) %*% as.matrix(coefQ)

  #Q0W:
  Q0<-Q
  Q0[,"A"] = 0

  task<-sl3::make_sl3_Task(Q0, covariates = covars, outcome = outcome, folds=new_folds)
  Q0W <- as.matrix(splitQ$predict(task)) %*% as.matrix(coefQ)

  #p1A:
  covars<-names(W)
  outcome<-names(A)

  task<-sl3::make_sl3_Task(g, covariates = covars, outcome = outcome, folds=new_folds)
  pA1 <- as.matrix(splitg$predict(task)) %*% as.matrix(coefg)

  B <- (A/pA1 - (1 - A)/(1 - pA1)) * (Y - QAW) + Q1W - Q0W

  #Maximize:
  Rule <- as.numeric(B > 0)

  K <- as.vector(abs(B))

  return(list(Y=Y$Y, A=A$A, QAW = QAW, Q1W = Q1W, Q0W = Q0W, pA1 = pA1,
              B = B$A, Rule = Rule, K = K$A))

}

############################################################################
#
############################################################################

#' Wrapper for split-specific Rule Estimation
#'
#' Function to make split-specific rule estimation more streamlined.
#'
#' @param folds user-specified list of folds; it should correspond to an element of \code{origami}.
#' @param Q data.frame containg the relevant outcome and covariates for estimating Q (E(Y|A,W)).
#' @param estSplt result object of running \code{estSplit}.
#' @param library list of \code{sl3} algorithms used in order to fit of the blip function.
#' @param nodes list contating names for the outcome, treatment, covariates and V-covariates as
#' observed in the provided data.
#'
#' @return An object of class \code{tmle3opttx}.
#' \describe{
#' \item{B}{}
#' \item{optA}{}
#' \item{coef}{}
#' \item{bSplit}{}
#' \item{x}{}
#' \item{y}{}
#' }
#'
#' @importFrom stats coef
#' @importFrom nnls nnls
#'
#' @export
#'

estBlip<-function(folds, estSplt, Q, library, nodes){

  #Observed data:
  Y<-Q[, nodes$Ynode, drop=FALSE]
  A<-Q[, nodes$Anode, drop=FALSE]
  W<-Q[, nodes$Wnodes, drop=FALSE]
  V<-Q[, nodes$Vnodes, drop=FALSE]

  bSplit<-origami::cross_validate(cv_split_b, folds, estSplt$estSplit, Y,A,V, library,
                                     .combine = F)

  #Construct SL prediction of the rule:
  #For now, just use non-negative linear least squares.

  #Combine all predictions of E(D1|V) for each validation fold:
  x<-do.call(rbind, bSplit$cvPred)
  y<-unlist(bSplit$B)

  fit_coef <- stats::coef(nnls::nnls(as.matrix(x), as.matrix(y)))
  fit_coef <- fit_coef/sum(fit_coef)

  pred<-as.matrix(x) %*% fit_coef

  #Actual rule:
  optA<-as.numeric(pred > 0)

  return(list(B=pred,optA=optA,coef=fit_coef,bSplit=bSplit,x=x,y=y))

}

#' Split-specific Rule Estimation
#'
#' Function to estimate split-specific predictions for the rule.
#'
#' @param fold one of the folds from the folds object.
#' @param split result of the cross-validated \code{cv_split} function.
#' @param Y data.frame object for outcome.
#' @param A data.frame object for treatment.
#' @param V data.frame object for covariates the rule is based on.
#' @param library list of \code{sl3} algorithms for the fit of the E(D1|V) function.
#'
#' @return An object of class \code{tmle3opttx}.
#' \describe{
#' \item{fullPred}{Split-specific fits for E(B|V), with full V-fold cross-validation and
#' Super Learner results.}
#' \item{cvPred}{Split-specific E(B|V), with prediction for only a specific fold that was
#' used as a validation fold.}
#' \item{valSet}{Validation samples used for each fold.}
#' \item{B}{Estimated B based on split-specific fits of Q and g and evaluated on all
#' validation samples. Returned B is of the dimension of the validation samples used for
#' the fold in question.}
#' }
#'
#' @importFrom sl3 make_sl3_Task
#' @importFrom origami fold_index
#'
#' @export
#'

cv_split_b<-function(fold, split, Y,A,V, library){

  v <- origami::fold_index()
  Q=cbind.data.frame(Y,V)

  #Get the corresponding B for the fold v
  #(B based on fit for fold arrangment v)
  B<-data.frame(B=split$B[[v]])
  folds<-split$folds

  #Estimate D1 based on covariates specified in V:
  res<-estInit(Y = B, X = V, folds = folds, library = library)

  #Grab just the v fit prediction:
  valset<-Q[fold$validation_set,]

  #Ok, because we set Y to be first, V all the rest.
  #It's ok this is not the proper outcome we want to estimate.
  covars <- names(valset[,-1,drop=FALSE])
  outcome <- names(valset[,1,drop=FALSE])

  task<-sl3::make_sl3_Task(valset, covariates = covars, outcome = outcome, folds=fold)

  #Extract split-specific bit for D1 estimated based on V:
  valset_res<-res$cvFit[[v]]$predict(task)
  row.names(valset_res)<-fold$validation_set

  return(list(fullPred=res,cvPred=valset_res,valSet=fold$validation_set,B=B[fold$validation_set,]))
}

