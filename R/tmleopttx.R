#' Optimal Individualized Treatment using gentmle
#'
#' This function estimates the optimal individualized treatment effect for a binary treatment A.
#' Estimation of the Optimal treatment rule is performed using the sl3 package (Super Learner algorithm)
#' and mean performance under the optimal treatment rule is estimated using CV-TMLE and gentmle.
#'
#' @param data data.frame object containing outcome labeled as "Y", treatment labeled as "A"
#' and covariates labeled as "W".
#' @param Wnodes vector of column names indicating covariates.
#' @param Anode column name of treatment.
#' @param Ynode column name of outcome.
#' @param Vnodes vector of column names to base the treatment on.
#' @param V number of cross-validation folds used.
#' @param stratifyAY logical: should we stratify the cross-validation based on (A,Y) pairs
#' @param Q_library list of \code{sl3} algorithms for the fit of E(Y|A,Cy) (Q)
#' @param g_library list of \code{sl3} algorithms for the fit of P(A|Ca) (g)
#' @param b_library list of \code{sl3} algorithms for the fit of the blip function.
#' @param gbounds bounds for the q estimates.
#' @param Qbounds bounds for the Q estimates.
#' @param verbose controls the verbosity of the output.
#'
#' @return An object of class \code{tmle3opttx}.
#' \describe{
#' \item{tmlePsi}{TMLE estimate of the mean under the user-specified \code{ruleA}.
#' In particular, it returns psi under the learned optimal regime, observed exposure,
#' A=1 and A=0.}
#' \item{tmleSD}{Standard deviation of the mean under user-specified \code{ruleA} estimated
#' using the TMLE. In particular, it returns the standard deviation of the mean under the
#' optimal regime, observed exposure, A=1 and A=0.}
#' \item{tmleCI}{Confidence interval of the mean under user-specified \code{ruleA}
#' estimated using TMLE. It returns CI under the optimal regime, observed exposure, A=1 and A=0.}
#' \item{IC}{Influence curve for the mean under user-specified \code{ruleA}.
#' It returns IC under the optimal regime, observed exposure, A=1 and A=0.}
#' \item{rule}{Used rule for the exposure.}
#' \item{steps}{Number of steps until convergence of the iterative TMLE for each rule.}
#' \item{initialData}{Initial estimates of g and Q, and observed A and Y.}
#' \item{tmleData}{Final updates estimates of g, Q and clever covariates.}
#' \item{all}{Full results of \code{tmleOPT}.}
#' }
#'
#' @export
#

tmleopttx <- function(data, V=10, stratifyAY=TRUE,
                       Wnodes=grep("^W", names(data), value = T),
                       Anode = "A", Ynode = "Y", Vnodes = Wnodes,
                       Q_library = list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet"),
                       g_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet"),
                       b_library=list("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_glmnet","Lrnr_randomForest","Lrnr_xgboost"),
                       gbounds=c(1e-4,1-1e-4), Qbounds=c(1e-4,1-1e-4), verbose = TRUE){

  #Make folds:
  if (stratifyAY) {
    #Learn all possible combinations of A and Y:
    AYstrata <- sprintf("%s %s", data[, Anode], data[, Ynode])
    #Stratified folds:
    folds <- origami::make_folds(strata_ids = AYstrata, V = V)
  } else {
    #Regular folds:
    folds <- origami::make_folds(V = V)
  }

  #Store all nodes and libraries in one file:
  nodes <- list(Wnodes = Wnodes, Anode = Anode, Ynode = Ynode, Vnodes = Vnodes)
  SL.library <- list(Q_library = Q_library, g_library = g_library, b_library = b_library)

  #Initial estimate of Q:
  if(verbose){message("Fitting Q:")}
  estQ<-estInit(Y=data[, Ynode, drop=FALSE],
                X=data[,c(nodes$Anode,nodes$Wnodes), drop=FALSE],
                folds=folds,library=SL.library$Q_library)
  estQ$valY<-bound(estQ$valY, Qbounds)

  #Initial estimate of g:
  if(verbose){message("Fitting g:")}
  estg<-estInit(Y=data[, Anode, drop=FALSE],
                X=data[,nodes$Wnodes, drop=FALSE],
                folds=folds,library=SL.library$g_library)
  estg$valY<-bound(estg$valY, gbounds)

  #Split-specific predictions:
  message("Generating split-specific predictions:")
  estSplit<-estSplit(folds,
                    Q=data,
                    g=data[,c(nodes$Anode,nodes$Wnodes), drop=FALSE],
                    estQ, estg, nodes)

  #Fit the rule:
  if(verbose){message("Estimate the rule:")}
  estBlp<-estBlip(folds, estSplit, Q=data, library=SL.library$b_library, nodes)

  #Run TMLE for all rules (old way of doing it):
  rule0 <- ruletmle(estSplit$valSplit$A, estSplit$valSplit$Y, estSplit$valSplit$pA1,
                    estSplit$valSplit$Q0W, estSplit$valSplit$Q1W, ruleA=0)
  rule1 <- ruletmle(estSplit$valSplit$A, estSplit$valSplit$Y, estSplit$valSplit$pA1,
                    estSplit$valSplit$Q0W, estSplit$valSplit$Q1W, ruleA=1)
  ruleA <- ruletmle(estSplit$valSplit$A, estSplit$valSplit$Y, estSplit$valSplit$pA1,
                    estSplit$valSplit$Q0W, estSplit$valSplit$Q1W, ruleA=unlist(estSplit$valSplit$A))
  ruleOpt <- ruletmle(estSplit$valSplit$A, estSplit$valSplit$Y, estSplit$valSplit$pA1,
                      estSplit$valSplit$Q0W, estSplit$valSplit$Q1W, ruleA=estBlp$optA)

  res <- list(rule0=rule0,rule1=rule1,ruleA=ruleA,ruleOpt=ruleOpt)

  #Extract psi for each rule:
  res_fin<-extract_res(res)

  return(list(tmlePsi=res_fin$tmlePsi,tmleSD=res_fin$tmleSD,tmleCI=res_fin$tmleCI,
              IC=res_fin$IC,rule=res_fin$rule,steps=res_fin$steps,initialData=res_fin$initialData,
              tmleData=res_fin$tmleData,all=res_fin$all))
}
