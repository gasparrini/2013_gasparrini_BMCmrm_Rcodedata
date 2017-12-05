###############################################################################
# Updated version of the code for the analysis in:
#
#   "Reducing and meta-analyzing estimates from distributed lag non-linear models"
#   Gasparrini and Armstrong 
#   BMC Medical Research Methodology - 2013
#   http://www.ag-myresearch.com/2013_gasparrini_bmcmrm.html
#
# Update: 05 December 2017
# * an updated version of this code, compatible with future versions of the
#   software, is available at:
#   http://www.ag-myresearch.com/2013_gasparrini_bmcmrm.html
###############################################################################

####################################################################
# FIRST STAGE
# - DEFINE THE CROSS-BASIS MATRICES FOR THE 3 MODELS
# - BUILD OBJECTS TO STORE THE RESULTS
# - RUN THE POISSON TIME SERIES MODELS
# - REDUCE THE FITTED MAIN MODEL TO SUMMARIES
# - STORE THE RESULTS
# COMPUTING TIME IS ~40SEC (IN A 2.66GHz-4GBRAM PC WITH WINDOWS)
####################################################################

####################################################################
# DEFINE THE CROSS-BASIS MATRICES
# NB: THE USER CAN MODIFY THE CHOICES BELOW TO RUN ALTERNATIVE MODELS

# MAIN MODEL
# - PREDICTOR SPACE: QUADRATIC SPLINE WITH SPECIFIC KNOT SELECTION
# - LAG SPACE: NATURAL CUBIC SPLINE WITH DF AT EQUALLY-SPACED LOG-VALUES
lag <- c(0,21)
bound <- colMeans(ranges)
varknots <- equalknots(bound,fun="bs",degree=2,df=4)
lagknots <- logknots(21,df=5,int=T)
argvar <- list(fun="bs",degree=2,knots=varknots,bound=bound)
arglag <- list(fun="ns",knots=lagknots)

# ALTERNATIVE MODELS
# - IDENTICAL BASIS FOR PREDICTOR SPACE
# LAG SPACE: CONSTANT FOR LAG 0-3 AND LAG 0-21
lag2 <- c(0,3)
lag3 <- c(0,21)
arglag2 <- arglag3 <- list(fun="strata",df=1)

####################################################################
# BUILT OBJECTS WHERE RESULTS WILL BE STORED
#   y- IS THE MATRIX FOR THE OUTCOME PARAMETERS
#   S- IS THE LISTS OF (CO)VARIANCE MATRICES

# OVERALL CUMULATIVE SUMMARIES
yall <- matrix(NA,length(data),4,dimnames=list(regions,paste("b",seq(4),sep="")))
yall2 <- yall3 <- yall

# PREDICTOR-SPECIFIC SUMMARIES FOR MAIN MODEL
yhot <- matrix(NA,length(data),5,dimnames=list(regions,paste("b",seq(5),sep="")))
ycold <- matrix(NA,length(data),5,dimnames=list(regions,paste("b",seq(5),sep="")))

# (CO)VARIANCE MATRICES
Sall <- vector("list",length(data))
names(Sall) <- regions
Shot <- Scold <- Sall2 <- Sall3 <- Sall

# Q-AIC
qaic <- qaic2 <- qaic3 <- 0

####################################################################
# RUN THE MODEL FOR EACH CITY

# WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
options(warn=-1)

# LOOP FOR CITIES
for(i in seq(data)) {
  
  # PRINT
  cat(i,"")

  # LOAD
  sub <- data[[i]]
  
  # DEFINE THE CROSS-BASES
  cb <- crossbasis(sub$tmean,lag=lag,argvar=argvar,arglag=arglag)
  cb2 <- crossbasis(sub$tmean,lag=lag2,argvar=argvar,arglag=arglag2)
  cb3 <- crossbasis(sub$tmean,lag=lag3,argvar=argvar,arglag=arglag3)
  
  # SET THE FIRST 21 RECORDS FOR cb2 AS MISSING
  # THIS MAKES THE 3 MODELS COMPARABLE THROUGH AIC (SAME OBS)
  cb2[0:21,] <- NA
  
  # RUN THE FIRST-STAGE MODELS
  mfirst <- glm(death ~ cb+dow+ns(time,df=10*14),family=quasipoisson(),sub)
  mfirst2 <- glm(death ~ cb2+dow+ns(time,df=10*14),family=quasipoisson(),sub)
  mfirst3 <- glm(death ~ cb3+dow+ns(time,df=10*14),family=quasipoisson(),sub)
  
####################################################################
  # REDUCTION TO SUMMARY ASSOCIATIONS

  # TO OVERALL CUMULATIVE SUMMARY
  # NB: CENTERING NOT REALLY NEEDED HERE, AS COEF-VCOV (EXTRACTED BELOW) IN THE 
  #   VAR SPACE DO NOT DEPEND ON CENTERING VALUE
  crall <- crossreduce(cb,mfirst,cen=17)
  crall2 <- crossreduce(cb2,mfirst2,cen=17)
  crall3 <- crossreduce(cb3,mfirst3,cen=17)

  # TO PREDICTOR-SPECIFIC SUMMARY FOR 22C AND 0C
  # NB: CENTERING NEEDED HERE, AS COEF-VCOV (EXTRACTED BELOW) IN THE LAG SPACE
  #   DO DEPEND ON CENTERING VALUE
  crhot <- crossreduce(cb,mfirst,type="var",value=22,cen=17)
  crcold <- crossreduce(cb,mfirst,type="var",value=0,cen=17)

####################################################################
  # STORE THE RESULTS
  
  # OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
  yall[i,] <- coef(crall)
  Sall[[i]] <- vcov(crall)
  
  # OVERALL CUMULATIVE SUMMARY FOR THE ALTERNATIVE MODELS
  yall2[i,] <- coef(crall2)
  yall3[i,] <- coef(crall3)
  Sall2[[i]] <- vcov(crall2)
  Sall3[[i]] <- vcov(crall3)
  
  # PREDICTOR-SPECIFIC SUMMARY FOR 22C (MAIN MODEL)
  yhot[i,] <- coef(crhot)
  Shot[[i]] <- vcov(crhot)
  # PREDICTOR-SPECIFIC SUMMARY FOR 0C (MAIN MODEL)
  ycold[i,] <- coef(crcold)
  Scold[[i]] <- vcov(crcold)
  
  # Q-AIC
  qaic[i] <- fqaic(mfirst)
  qaic2[i] <- fqaic(mfirst2)
  qaic3[i] <- fqaic(mfirst3)
  
}

####################################################################

# TEST: REDUCTION OF ALTERNATIVE MODELS TO THE SPACE OF THE PREDICTOR RETURNS
# THE SAME PARAMETERS APART FROM SCALING (SUMMED UPON 22 LAGS)
coef(crosspred(cb3,mfirst3,cen=17))
coef(crossreduce(cb3,mfirst3,cen=17))/22

# GRAND Q-AIC
sum(qaic) ; sum(qaic2) ; sum(qaic3)

# RESET WARNING
options(warn=0)

#
