\name{HiCblock-package}
\alias{HiCblock-package}
\alias{HiCblock}
\docType{package}
\title{
\packageTitle{HiCblock}
}
\description{
\packageDescription{HiCblock}
}
\author{
Raphael Mourad

Maintainer: Raphael Mourad, raphael.mourad@ibcg.biotoul.fr
}
\references{
TO ADD when published.
}

\examples{
# Load data
# The Hi-C matrix is at 20kb resolution (low resolution only for example)
data(dataExample)
genomicFeatureList.GR=dataExample$GenomicFeatureList.GR
annotNames=dataExample$AnnotNames
HTCList=dataExample$HTC
distInter=c(100e3,140e3)
IBP=c("BEAF32","dCTCF","dTFIIIC","GAF","SuHw")

# Annotate Hi-C data with genomic features
HRPD=HiCblockProcData(genomicFeatureList.GR, annotNames, HTCList, distInter,verbose=TRUE)

# Model 1
modelBlock1=as.formula(paste0("Count~logDist+len+GC+map+BEAF32"))
HRM_Block1=HiCblockModel(HRPD,modelBlock1,"BEAF32",regressionMode="NB")
print(HRM_Block1)

# Model 2
modelBlock2=as.formula(paste0("Count~logDist+len+GC+map+",paste(IBP,collapse='+')))
HRM_Block2=HiCblockModel(HRPD,modelBlock2,IBP,regressionMode="NB")
print(HRM_Block2)


}
