#' Inverse variance weighted MR 
#'
#'  R function IVW was obtained from the paper by Bowden, J., Davey Smith, G., & Burgess, S. (2015). Mendelian randomization with invalid instruments: effect estimation and bias
#'  detection through Egger regression. International Journal of Epidemiology, 44(2), 512-525
#' @param BetaYG Beta of outcome
#' @param BetaXG Beta of exposure
#' @param seBetaYG standard error of outcome
#' @param data dataset for exposure and outcome
#' @keywords Mendelian Randomization
#' @export
#' @examples
#' IVW()

IVW<-function(BetaYG, BetaXG,seBetaYG, data){
  
  data$Weights = 1/data[, seBetaYG]^2
  IVWfit   = summary(lm(as.formula(paste0(BetaYG, " ~ -1+ ",BetaXG)),weights=Weights,data=data))
  DF      = nrow(data)-1
  IVWBeta = IVWfit$coef[1,1]
  SE      = IVWfit$coef[1,2]/min(1,IVWfit$sigma)
  IVW_p   = 2*(1-pt(abs(IVWBeta/SE),DF))
  IVW_CI  = IVWBeta + c(-1,1)*qt(df=DF, 0.975)*SE
  
  # IVWResults = (point estimate, corrected standard error,
  # 95% Confidence interval, t-statistic, p-value)
  res=c(IVWBeta,SE,IVWBeta/SE,IVW_p)
  names(res)=c("Beta","SE","T-stat","P")
  return(res)
}