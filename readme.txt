###  MR_pleio is a R program to perform Iterative Memdelian Randomization and Pleiotropy analysis (IMRP, Zhu et al. An iterative approach to detect pleiotropy and perform Mendelian 
###  Randomization analysis using GWAS summary statistics", Bioinformatics, 2020). To call the program, the summary statistics of instrumental variables for the exposure and outcome 
###  need to be combined into a single data set, with matched reference and effect alleles. We suggested to standardize the summary statistics before performing IMRP, although this is 
###  not a necessary step. We provide an example of MR analysis of HDL on CAD analyzed in the manuscript. The HLD data was downloaded from http://csg.sph.umich.edu/abecasis/
###  public/lipids2013/; CAD data was downloaded from http://www.cardiogramplusc4d.org/data-downloads/. If you have any questions, please contact Xiaofeng Zhu (xxz10@case.edu). 

To run the program, type

library(devtools)
install_github("XiaofengZhuCase/IMRP")
library("IMRP")

### The combined HDL and CAD data is named HDLCAD1.csv with each columns representing "SNP","Chr","Position","A1","A2","HDL_beta","HDL_se","HDL_N", "CAD_freq","CAD_beta","CAD_se". 
### Here HDL is a continuous trait and CAD is a binary trait. This approach is the same as suggested in the software MRmix (Qi and Chatterjee, Nat Commun 2019)

### Standardizing summary statistics

HDLCAD1$x1=HDLCAD1$HDL_beta/HDLCAD1$HDL_se/sqrt(HDLCAD1$HDL_N)
HDLCAD1$x2=HDLCAD1$CAD_beta*sqrt(2*HDLCAD1$CAD_freq*(1-HDLCAD1$CAD_freq))
HDLCAD1$x1_se=1/sqrt(HDLCAD1$HDL_N)
HDLCAD1$x2_se=HDLCAD1$CAD_se*sqrt(2*HDLCAD1$CAD_freq*(1-HDLCAD1$CAD_freq))

### to call MR_pleio, 

HDLCAD1_out=MR_pleio("x2","x1","x2_se","x1_se",as.data.frame(HDLCAD1),SignifThreshold=0.05,rho=0.03, method="IVW")

### In above function, 0.05 is the threshold to define outliers in pleiotropic variants. This value can be changed depending on that more or less potential pleiiotropic variants will be excluded. 
### rho=0.03 is the correlation coefficient estimated using the genome wide summary statistics of HDL and CAD.

### output
 
> HDLCAD1_out
$CausalEstimate
       Beta 
-0.07497658  ### causal effect estimate

$SdCausalEstimate
        SE ### standard error of causal effect estimate
0.02682728 

$Causal_p   ### p value for testing causal effect=0
       Beta 
0.005193358 

$SNPPvalue          ### P-values of testing for pleiotropic effects of each instrumental variables
  [1] 7.540308e-01 1.213973e-01 1.765340e-02 7.072475e-01 1.371316e-01
  [6] 1.614926e-01 4.872331e-01 4.480403e-05 1.724450e-04 4.013943e-01
 [11] 8.098048e-01 5.771981e-01 1.376007e-02 6.629121e-01 3.385389e-01
 [16] 8.395957e-02 9.293667e-02 8.796184e-01 6.673686e-01 6.327465e-01
 [21] 6.923009e-01 5.909383e-22 6.007826e-01 1.087855e-05 6.387741e-01
 [26] 8.385293e-01 6.110703e-01 8.234167e-02 1.686317e-02 2.583850e-01
 [31] 4.460922e-01 3.500133e-01 4.801952e-04 7.470540e-04 8.292815e-03
 [36] 7.654568e-03 4.717011e-01 1.604432e-01 1.487996e-01 1.520034e-02
 [41] 4.610298e-01 7.022568e-01 6.825129e-03 2.505968e-08 5.600326e-01
 [46] 1.813336e-01 4.364929e-01 3.351774e-01 5.669551e-01 1.036621e-01
 [51] 3.854098e-03 3.962321e-01 8.243925e-01 2.062745e-02 1.742970e-02
 [56] 3.078318e-03 6.822560e-06 9.713822e-01 3.447723e-01 4.875966e-01
 [61] 6.659414e-02 2.844937e-01 4.063577e-01 2.046420e-01 3.032017e-02
 [66] 3.145641e-02 6.247301e-01 9.638166e-01 8.065951e-01 8.263397e-02
 [71] 3.664172e-01 7.131909e-01 9.402118e-01 6.119039e-01 6.127570e-01
 [76] 5.716073e-01 5.651548e-02 8.898345e-01 7.170033e-02 1.803775e-01
 [81] 2.607218e-02 7.432054e-03 6.443119e-01 4.776510e-05 5.723056e-01
 [86] 5.056964e-01 5.830633e-01 9.254104e-04 8.278921e-02 7.677411e-03
 [91] 8.723933e-01 7.311049e-01 8.262254e-01 7.928084e-01 3.879858e-01
 [96] 9.880206e-01 1.752683e-03 1.509789e-02 9.984555e-03 3.057486e-03
[101] 1.711400e-01 2.661010e-01 3.290598e-01 4.517062e-03 1.185573e-02
[106] 6.981939e-01 1.150645e-01 7.323114e-01 2.170176e-01 9.441846e-01
[111] 5.357276e-01 3.590136e-01 3.993563e-04 8.469080e-03 6.381885e-03
[116] 1.837256e-09 1.260571e-07 3.778442e-01 1.708251e-01 2.067299e-01
[121] 1.443810e-03 2.713742e-02 2.776027e-01 3.117851e-02 6.351858e-01
[126] 2.614562e-01 7.124180e-03 2.271693e-01 8.221943e-02 3.328679e-02
[131] 4.093434e-01 8.268081e-01 2.037193e-02 5.976885e-03 8.657520e-02
[136] 1.010016e-01 2.227204e-01 7.530752e-01 6.208741e-04 5.392791e-01
[141] 6.976097e-01 4.913947e-05 8.697273e-01

$PleioOutlier            #### pleiotropic variants outliers with P<0.05
 [1]   3   8   9  13  22  24  29  33  34  35  36  40  43  44  51  54  55  56  57
[20]  65  66  81  82  84  88  90  97  98  99 100 104 105 113 114 115 116 117 121
[39] 122 124 127 130 133 134 139 142

$GlobalPvalue_pre       #### global P value for testing if tere are any pleiotropic variants among the instrumental variables. 
[1] 1.393238e-65

$GlobalPvalue_aft       ### global P value for testing if there are still pleiotropic variants after excluding pleiotropic variants.
[1] 0.6019805

$chisquare              ### chisqure statistics for testing if tere are any pleiotropic variants among the instrumental variables, following a chisqure distribution with df=number of the instrumental variables-1. 
[1] 647.4692
