#' Pleiotropy test and Mendelian Randomization analysis
#'
#' Tpleio performs pleiotropy test and one-step Mendelian Randomization analysis 
#' @param BetaOutcome beta of outcome
#' @param BetaExposure beta of exposure
#' @param SdOutcome standard error of outcome
#' @param SdExposure standard error of exposure 
#' @param data dataset for storing summary statistics
#' @param rho correlation between betas of outcome and exposure. It can be calculated using the genome wide summary statistics. (Zhu et al. AJHG 2015)
#' @param methods either IVW or MREgger
#' @param Outlier the index of pleiotropic variant outliers
#' @param SignifThreshold The threshold to define pleiotropic variants
#' @keywords Mendelian randomization
#' @export
#' @examples
#' Tpleio()

Tpleio<-function(BetaOutcome, BetaExposure, SdOutcome, SdExposure,data,rho,methods,Outlier, SignifThreshold){
  data$Weights <- 1/data[, SdOutcome]^2
  
  if(!is.null(Outlier)){subdata=data[-Outlier,]
  } else if(is.null(Outlier)){subdata=data }
  
  if(methods=="MR_Egger"){
    MR_Egger=MREgger(BetaYG=BetaOutcome, BetaXG=BetaExposure, seBetaYG =SdOutcome, data=subdata)
    CausalEstimate<-MR_Egger[2,1]
    SdCausalEstimate<-MR_Egger[2,2]
  } else if(methods=="IVW"){
    ivw=IVW(BetaYG=BetaOutcome, BetaXG=BetaExposure, seBetaYG =SdOutcome, data=subdata)
    CausalEstimate<-ivw[1]
    SdCausalEstimate<-ivw[2]
  }
  
  X<-data[, c(BetaOutcome)]-CausalEstimate*data[, c(BetaExposure)]
  Y<-sqrt(data[, c(SdOutcome)]^2+CausalEstimate^2*data[, c(SdExposure)]^2+data[,c(BetaExposure)]^2*SdCausalEstimate^2
          -2*CausalEstimate*rho*data[,c(SdExposure)]*data[,c(SdOutcome)])
  
  GlobalP_pre<-pchisq(sum((X/Y)^2), df=nrow(data)-1, lower.tail=FALSE)
  pval<- apply(cbind(2*pnorm(-abs(X/Y)), 1), 1, min)
  Outlier=which(pval<=SignifThreshold)
  Xaft <- X[-Outlier]
  Yaft <- Y[-Outlier]
  GlobalP_aft<-pchisq(sum((Xaft/Yaft)^2), df=length(Xaft)-1, lower.tail=FALSE)

  return(list(CausalEstimate=CausalEstimate,SdCausalEstimate=SdCausalEstimate,Causal_p=2*pnorm(-abs(CausalEstimate/SdCausalEstimate)), SNPPvalue=pval,PleioOutlier=Outlier,GlobalPvalue_pre=GlobalP_pre,GlobalPvalue_aft=GlobalP_aft,chisquare=sum((X/Y)^2)))
}