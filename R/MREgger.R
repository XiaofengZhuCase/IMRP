#' Mendelian randomization analysis MREgger
#'
#' R function MREgger was obtained from the paper by Bowden, J., Davey Smith, G., & Burgess, S. (2015). Mendelian randomization with invalid instruments: effect estimation and bias
#' detection through Egger regression. International Journal of Epidemiology, 44(2), 512-525

#' @param BetaYG beta of outcome
#' @param BetaXG beta of exposure
#' @param seBetaYG standard error of outcome
#' @param data dataset of exposure and outcome
#' @keywords Mendelian randomization 
#' @export
#' @examples
#' MREgger()

MREgger<-function(BetaYG, BetaXG,seBetaYG, data){
  data$Weights = 1/data[, seBetaYG]^2
  data[, c(BetaYG,BetaXG)] <-   data[, c(BetaYG,BetaXG)] * sign(data[,BetaXG])
  MREggerFit  = summary(lm(as.formula(paste0(BetaYG, " ~ ",BetaXG)),weights=Weights,data=data))
  
  # Inference with correct standard errors
  
  MREggerBeta0   = MREggerFit$coef[1,1]
  MREggerBeta1   = MREggerFit$coef[2,1]
  SE0            = MREggerFit$coef[1,2]/min(1,MREggerFit$sigma)
  SE1            = MREggerFit$coef[2,2]/min(1,MREggerFit$sigma)
  DF             = nrow(data)-2
  MRBeta0_p      = 2*(1-pt(abs(MREggerBeta0/SE0),DF))
  MRBeta1_p      = 2*(1-pt(abs(MREggerBeta1/SE1),DF))
  MRBeta0_CI     = MREggerBeta0 + c(-1,1)*qt(df=DF, 0.975)*SE0
  MRBeta1_CI     = MREggerBeta1 + c(-1,1)*qt(df=DF, 0.975)*SE1
  
  # MREggerResults = (point estimate, corrected standard error,
  # 95% Confidence interval, t-statistic, p-value) for
  # intercept (row 1) and slope (row 2).
  
  MREggerResults     = matrix(nrow = 2,ncol = 4)
  MREggerResults[1,] = c(MREggerBeta0,SE0,MREggerBeta0/SE0,MRBeta0_p)
  MREggerResults[2,] = c(MREggerBeta1,SE1,MREggerBeta1/SE1,MRBeta1_p)
  colnames(MREggerResults)=c("Beta","SE","T-stat","P")
  rownames(MREggerResults)=c("Intercept",BetaXG)
  return(MREggerResults)
}
