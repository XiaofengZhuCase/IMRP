#' MR_pleio performs iterative Mendelian Randomization and Pleiotropy analysis (IMRP)
#'
#' MR_pleio is a R program to perform Iterative Memdelian Randomization and Pleiotropy analysis (IMRP, Zhu et al. An iterative approach to detect pleiotropy and perform Mendelian 
#' Randomization analysis using GWAS summary statistics", Bioinformatics, 2020). To call the program, the summary statistics of instrumental variables for the exposure and outcome 
#' need to be combined into a single data set, with matched reference and effect alleles. We suggested to standardize the summary statistics before performing IMRP, although this is 
#' not a necessary step. We provide an example of MR analysis of HDL on CAD analyzed in the manuscript. The HLD data was downloaded from http://csg.sph.umich.edu/abecasis/
#' public/lipids2013/; CAD data was downloaded from http://www.cardiogramplusc4d.org/data-downloads/. If you have any questions, please contact Xiaofeng Zhu (xxz10@case.edu).
#' @param BetaOutcome the beta of outcome
#' @param BetaExposure the beta of exposure
#' @param SdOutcome the standard error of beta for the outcome
#' @param SdExposure the standard error of beta for the exposure
#' @param data the dataset of exposure and outcome
#' @param SignifThreshold the threshold to define pleiotropic outlier
#' @param rho the correlation between exposure beta and outcome beta. It can be estimated by the genome wide summary statistics. (Zhu et al. AJHG 2015)
#' @param methods The IVW or MEegger approaches in iterative algorithm
#' @keywords Mendelian Randomization
#' @export
#' @examples
#' MR_pleio()

MR_pleio <- function(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data, SignifThreshold,rho,methods="IVW"){
  Outlier=NULL
  if(SignifThreshold > 1)
    stop("The significance threshold cannot be greater than 1")
  
  if(length(BetaExposure) != length(SdExposure))
    stop("BetaExposure and SdExposure must have the same number of elements")
  
  if(class(data)[1] != "data.frame")
    stop("data must be an object of class data.frame, try to converse data to a data.frame \'data = as.data.frame(data)\'")
  
  # Functions
  
  
  # 0- Transforming the data + checking number of observations
  data <- data[, c(BetaOutcome, BetaExposure, SdOutcome, SdExposure)]
  data <- data[rowSums(is.na(data)) == 0, ]
  data[, c(BetaOutcome, BetaExposure)] <- data[, c(BetaOutcome, BetaExposure)] * sign(data[, BetaExposure[1]])
  data$Weights <- 1/data[, SdOutcome]^2
  
  if(nrow(data) <= length(BetaExposure) + 2)
    stop("Not enough intrumental variables")
  
    # 1- Computing the observed residual sum of squares (RSS)
    betadiff=1.0

    Test<-Tpleio(BetaOutcome, BetaExposure,SdOutcome, SdExposure, data=data,rho,methods="MR_Egger", Outlier, SignifThreshold)
    beta0=Test$CausalEstimate
    Outlier<- Test$PleioOutlier
    i=0
    while (betadiff>0.000001 & length(Outlier)<nrow(data) & i<500) 
    {
      i=i+1
      if(length(Outlier)==0) Outlier=NULL
      Test<-Tpleio(BetaOutcome, BetaExposure,SdOutcome, SdExposure, data=data,rho,methods, Outlier, SignifThreshold)
      Outlier<- Test$PleioOutlier
      betadiff=abs(beta0-Test$CausalEstimate)
      beta0=Test$CausalEstimate
    }
    return(Test)
}