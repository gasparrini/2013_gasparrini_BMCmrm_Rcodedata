###############################################################################
# Updated version of the code for the analysis in:
#
#   "Reducing and meta-analyzing estimates from distributed lag non-linear models"
#   Gasparrini and Armstrong 
#   BMC Medical Research Methodology - 2013
#   http://www.ag-myresearch.com/bmcmrm2013.html
#
# Update: 14 March 2016
# For any problem with this code, please contact antonio.gasparrini@lshtm.ac.uk
# Please refer to the original code for any copyright issue
#
#  See www.ag-myresearch.com for future updates
###############################################################################

###############################################################################
# CREATE 3 OBJECTS:
#
# 1) A VECTOR WITH NAMES OF REGIONS OF ENGLAND AND WALES
#
# 2) A LIST WITH THE DATA FOR EACH REGION, INCLUDING:
#   - DATE, YEAR, MONTH, DAY, TIME, DAY OF THE YEAR, DAY OF THE WEEK
#   - REGION NUMBERS AND NAMES
#   - MEAN, MINIMUM AND MAXIMUM TEMPERATURE
#   - DEW-POINT TEMPERATURE AND RELATIVE HUMIDITY
#   - MORTALITY (ALL-CAUSE)
#
# 3) A FUNCTION TO COMPUTE THE Q-AIC
#
###############################################################################

# LOAD PACKAGES (ASSUMED ALREADY INSTALLED)
library(dlnm) ; library(mvmeta) ; library(splines)

# CHECK VERSION OF THE PACKAGE
  if(packageVersion("dlnm")<"2.2.0")
    stop("update dlnm package to version >= 2.2.0")

# LOAD THE DATASET
regEngWales <- read.csv("regEngWales.csv",row.names=1)
dim(regEngWales)
head(regEngWales)

# REGIONS
regions <- as.character(unique(regEngWales$regnames))
  
# CREATE A LIST WITH THE REGIONAL SERIES
data <- lapply(regions,function(x) regEngWales[regEngWales$regnames==x,])
names(data) <- regions
m <- length(regions)

# TEMPERATURE RANGES
ranges <- t(sapply(data, function(x) range(x$tmean,na.rm=T)))

####################################################################

# FUNCTION TO COMPUTE THE Q-AIC IN QUASI-POISSON MODELS
fqaic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + 2*summary(model)$df[3]*phi
  return(qaic)
}

#
