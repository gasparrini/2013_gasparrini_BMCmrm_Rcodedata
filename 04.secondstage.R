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

####################################################################
# SECOND STAGE
# - RUN THE MULTIVARIATE META-ANALYTICAL MODELS WITH mvmeta
# - CREATE BASIS VARIABLES USING onebasis, TO BE USED FOR PREDICTION
# - OBTAIN PREDICTIONS THROUGH crosspred (dlnm)
####################################################################

####################################################################
# PERFORM MULTIVARIATE META-ANALYSIS

# LOAD THE PACKAGES (mvmeta PACKAGE IS ASSUMED TO BE INSTALLED)
library(mvmeta)

# SELECT THE ESTIMATION METHOD
method <- "reml"
# IN THE CURRENT VERSION, SET control=list(showiter=T) TO 
#   INSPECT THE OPTIMIZATION SEARCH

# OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
mvall <- mvmeta(yall~1,Sall,method=method)
summary(mvall)

# OVERALL CUMULATIVE SUMMARY FOR THE ALTERNATIVE MODELS
mvall2 <- mvmeta(yall2~1,Sall2,method=method)
summary(mvall2)
mvall3 <- mvmeta(yall3~1,Sall3,method=method)
summary(mvall3)

# PREDICTOR-SPECIFIC SUMMARY FOR 22C (MAIN MODEL)
mvhot <- mvmeta(yhot~1,Shot,method=method)
summary(mvhot)
# NOTE THE PROBLEM FOR ESTIMATED (CO)VARIANCE MATRIX

# PREDICTOR-SPECIFIC SUMMARY FOR 0C (MAIN MODEL)
mvcold <- mvmeta(ycold~1,Scold,method=method)
summary(mvcold)

####################################################################
# CREATE BASES FOR PREDICTION

# BASES OF TEMPERATURE AND LAG USED TO PREDICT, EQUAL TO THAT USED FOR ESTIMATION
# COMPUTED USING THE ATTRIBUTES OF THE CROSS-BASIS USED IN ESTIMATION
xvar <- seq(bound[1],bound[2],by=0.1)
bvar <- do.call("onebasis",c(list(x=xvar),attr(cb,"argvar")))
xlag <- 0:210/10
blag <- do.call("onebasis",c(list(x=xlag),attr(cb,"arglag")))

####################################################################
# REGION-SPECIFIC FIRST-STAGE SUMMARIES

regall <- lapply(seq(nrow(yall)),function(i) crosspred(bvar,coef=yall[i,],
  vcov=Sall[[i]],model.link="log",cen=17))
reghot <- lapply(seq(nrow(yhot)),function(i) crosspred(blag,coef=yhot[i,],
  vcov=Shot[[i]],model.link="log",cen=17))
regcold <- lapply(seq(nrow(ycold)),function(i) crosspred(blag,coef=ycold[i,],
  vcov=Scold[[i]],model.link="log",cen=17))

####################################################################
# PREDICTION FOR A GRID OF TEMPERATURE AND LAG VALUES

# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL
cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall),
  model.link="log",by=0.1,from=bound[1],to=bound[2],cen=17)

# OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR ALTERNATIVE MODELS
cpall2 <- crosspred(bvar,coef=coef(mvall2),vcov=vcov(mvall2),
  model.link="log",by=0.1,from=bound[1],to=bound[2],cen=17)
cpall3 <- crosspred(bvar,coef=coef(mvall3),vcov=vcov(mvall3),
  model.link="log",by=0.1,from=bound[1],to=bound[2],cen=17)

# PREDICTOR-SPECIFIC SUMMARIES FOR 22C (MAIN MODEL)
cphot <- crosspred(blag,coef=coef(mvhot),vcov=vcov(mvhot),
  model.link="log",at=0:210/10,cen=17)

# PREDICTOR-SPECIFIC SUMMARIES FOR 22C (MAIN MODEL)
cpcold <- crosspred(blag,coef=coef(mvcold),vcov=vcov(mvcold),
  model.link="log",at=0:210/10,cen=17)

#

