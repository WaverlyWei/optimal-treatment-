#' Wrapper function for the main TMLE procedure
#'
#' This function performs all the necessary steps in order to call the main TMLE function
#' in order to obtain a targeted estimate of the mean under the optimal individualized treatment
#' regime. It also performs all the calculations necessary to obtain inference for the parameter.
#'
#' @param obsA Observed exposure (A).
#' @param obsY Observed outcome (Y).
#' @param pA1 Estimate of g, (P(A=1|W)) obtained by calling \code{estSplit}.
#' @param Q0W Estimate of Q(0,W), (E(Y|A=0,W)), obtained by calling \code{estSplit}.
#' @param Q1W Estimate of Q(1,W), (E(Y|A=1,W)) obtained by calling \code{estSplit}.
#' @param ruleA user-specified rule for exposure under which we want to estimate the
#' mean of the outcome.
#' @param maxIter maximum number of iterations for the iterative TMLE.
#' @param Qbounds bounds for the Q estimates.
#'
#' @return An object of class \code{tmle3opttx}.
#' \describe{
#' \item{tmlePsi}{Context-specific mean of the outcome under user-specified \code{ruleA} estimated
#' using TMLE methodology.}
#' \item{tmleSD}{Standard deviation of the context-specific mean of the outcome under
#' user-specified \code{ruleA} estimated using TMLE methodology.}
#' \item{tmleCI}{Confidence interval of the context-specific mean of the outcome under
#' user-specified \code{ruleA} estimated using TMLE methodology.}
#' \item{IC}{Influence curve for the context-specific parameter under user-specified \code{ruleA}.}
#' \item{rule}{Used rule for the exposure.}
#' \item{steps}{Number of steps until convergence of the iterative TMLE.}
#' \item{initialData}{Initial estimates of g and Q, and observed A and Y.}
#' \item{tmleData}{Final updates estimates of g, Q and clever covariates.}
#' \item{all}{Full results of \code{tmleOPT}.}
#' }
#'
#' @importFrom stats qnorm
#'
#' @export
#

ruletmle<-function(obsA,obsY,pA1,Q0W,Q1W,ruleA,maxIter=1000,Qbounds=c(1e-4,1-1e-4)){

  #Get samples that follow the rule
  #(1 for samples that follow the rule in the observed data)
  A<-as.numeric(unlist(obsA)==ruleA)

  #Probability of success
  #(if rule is A=0, get p(A=0))
  pA<-unlist(pA1)*ruleA + (1-unlist(pA1))*(1-ruleA)

  #Need E(Y|A=0,W) if rule is 0
  Q<-unlist(Q1W)*ruleA + unlist(Q0W)*(1-ruleA)

  #Create appropriate data.frame and run tmle
  tmledata<-data.frame(A=A, Y=unlist(obsY), gk=pA, Qk=Q)
  res<-tmleOIT(tmledata, maxIter, Qbounds)

  #Get inference:
  psi<-res$psi
  sd<-sqrt(res$ED2)/sqrt(length(obsA))
  lower <- psi - stats::qnorm(0.975) * sd
  upper <- psi + stats::qnorm(0.975) * sd

  #Get used rule:
  if(length(ruleA)==1){
    ruleA<-rep(ruleA,nrow(tmledata))
  }

  return(list(tmlePsi=psi, tmleSD=sd, tmleCI=c(lower,upper), rule=ruleA,
              IC=res$IC, steps=res$steps, initialData=tmledata, tmleData=res$tmledata,all=res))
}

#' Main TMLE Calculations for the Mean under the Optimal Individualized Treatment
#'
#' This function performs all the main TMLE targeting steps for the mean under the
#' optimal individualized treatment.
#'
#' @param tmledata \code{data.frame} containing all observed values for the A and Y node,
#' as well as the estimate of g, Q.
#' @param maxIter maximum number of iterations for the iterative TMLE.
#' @param Qbounds bounds for the Q estimates.
#'
#' @return An object of class \code{tmle3opttx}.
#' \describe{
#' \item{tmledata}{Final updates estimates of g, Q and clever covariates.}
#' \item{psi}{Average treatment effect estimated using TMLE.}
#' \item{steps}{Number of steps until convergence of the iterative TMLE.}
#' \item{IC}{Influence function.}
#' \item{ED}{Mean of the final influence function.}
#' \item{ED2}{Mean of the squared final influence function.}
#' }
#'
#' @export
#

tmleOIT <- function(tmledata, maxIter=1000, Qbounds=c(1e-4,1-1e-4)){

  order <- 1/nrow(tmledata)

  #Initial estimate:
  eststep <- estOIT(tmledata)

  for (iter in seq_len(maxIter)){

    updatestep <- updateOIT(eststep$tmledata, Qbounds)
    eststep <- estOIT(updatestep$tmledata)

    ED <- sapply(eststep$Dstar, mean)

    if (all(abs(ED) < order)) {
      converge <- T
      break
    }
  }

  ED2 <- sapply(eststep$Dstar, function(x) mean(x^2))

  return(list(tmledata = eststep$tmledata, psi = eststep$ests, steps = iter,
              IC = eststep$Dstar, ED = ED, ED2 = ED2))
}

#################
# gentmle setup
#################

#' Estimate function of the gentmle setup
#'
#' This function estimates the mean under an user-specified rule (or optimal individualized rule)
#' and the coresponding relevant parts of the influence curve for the single intervention.
#'
#' @param tmledata \code{data.frame} containing all observed values for the A and Y node,
#' and either initial or updated estimates for g and Q.
#'

estOIT <- function(tmledata) {

  psi <- mean(tmledata$Qk)

  tmledata$H1 <- with(tmledata, (1/gk))
  tmledata$HA <- with(tmledata, (A*H1))

  # influence curves
  Dstar_psi <- with(tmledata, HA * (Y - Qk) + Qk - psi)

  list(tmledata = tmledata, ests = c(psi = psi), Dstar = list(Dstar_psi = Dstar_psi))

}

#' Update function of the gentmle setup
#'
#' This function updates the Q part of the likelihood using the specified fluctuation model for the
#' mean under an user-specified rule (or optimal individualized rule)
#'
#' @param data \code{data.frame} containing all observed values for the A and Y node,
#' and either initial or updated estimates for g and Q.
#' @param Qbounds bounds for the Q estimates.
#'
#' @importFrom stats plogis qlogis
#'

updateOIT <- function(data, Qbounds) {

  #Check Qk is still bounded, some subjects receive treatment
  subset <- with(data, which(0 < Qk & Qk < 1 & A == 1))
  eps <- 0

  if (length(subset) > 0) {
    #fluctuate Q
    data$Qktrunc <- bound(with(data, Qk), Qbounds)
    qfluc <- fluctuate(data, Y ~ -1 + offset(stats::qlogis(Qktrunc)) + HA)
    eps <- qfluc$eps
    data$Qk <- with(data, stats::plogis(stats::qlogis(Qktrunc) + HA * eps))
  }

  list(tmledata = data, coefs = c(eps))

}

#' Fluctuate function of the gentmle setup
#'
#' This function uses a logistic parametric submodel to fluctuate the initial fit of Q.
#'
#' @param tmledata \code{data.frame} containing all observed values for the A and Y node,
#' and either initial or updated estimates for g and Q.
#' @param flucmod fluctuation model used.
#' @param subset subset of the data used for fluctuation.
#'
#' @importFrom stats glm
#'

fluctuate <- function(tmledata, flucmod, subset = seq_len(nrow(tmledata))) {
  suppressWarnings({
    fluc <- stats::glm(flucmod, data = tmledata[subset, ], family = "binomial")
  })
  list(eps = stats::coef(fluc)[1])
}
