#' Pleiotropy test
#'
#' The R function only performs pleiotropy test when CausalEstimate is known.
#' @param BetaOutcome beta of outcome
#' @param BetaExposure beta of exposure
#' @param SdOutcome standard error of outcome
#' @param SdExposure standard error of exposure
#' @param data dataset to store summary statistics
#' @param rho correlation between betas of outcome and exposure. It can be calculated using the genome wide summary statistics. (Zhu et al. AJHG 2015)
#' @param CausalEstimate Causal effect estimated from MR analysis
#' @param SdCausalEstimate the standard error of the causal effect estimate from MR  
#' @keywords Pleiotropy test
#' @export
#' @examples
#' SearchPleio()

# Find pleiotropic variants
SearchPleio<-function(BetaOutcome, BetaExposure, SdOutcome, SdExposure,data,rho,CausalEstimate, SdCausalEstimate){
  
  X<-data[, c(BetaOutcome)]-CausalEstimate*data[, c(BetaExposure)]
  Y<-sqrt(data[, c(SdOutcome)]^2+CausalEstimate^2*data[, c(SdExposure)]^2+data[,c(BetaExposure)]^2*SdCausalEstimate^2
          -2*CausalEstimate*rho*data[,c(SdExposure)]*data[,c(SdOutcome)])
  
  pval<- apply(cbind(2*pnorm(-abs(X/Y)), 1), 1, min)

  return(list(pleio_p=pval))
}
