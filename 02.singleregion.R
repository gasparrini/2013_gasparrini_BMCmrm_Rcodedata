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
# RUN AN EXAMPLE OF FIRST-STAGE ANALYSIS FOR A SINGLE REGION
####################################################################

# SELECT REGION
reg <- "N-East"

# ARGUMENTS AND LISTS FOR CROSS-BASIS DEFINITION
bound <- colMeans(ranges)
varknots <- equalknots(bound,fun="bs",degree=2,df=4)
lagknots <- logknots(21,df=5,int=T)
argvar <- list(fun="bs",degree=2,knots=varknots,bound=bound)
arglag <- list(fun="ns",knots=lagknots)

# BASIS FOR TEMPERATURE:
# - QUADRATIC SPLINE FOR PREDICTOR, WITH SPECIFIC KNOT SELECTION
# - NATURAL CUBIC SPLINE FOR LAG, WITH DF AT EQUALLY-SPACED LOG-VALUES
# WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
suppressWarnings(
cb <- crossbasis(data[[reg]]$tmean,lag=21,argvar=argvar,arglag=arglag)
)
summary(cb)

# RUN THE MODEL
model <- glm(death ~ cb + dow + ns(time,df=10*14),
  family=quasipoisson(),data[[reg]])

# PREDICTION USING:
#   crosspred FOR BI-DIMENSIONAL RELATIONSHIP
#   crossreduce FOR UNI-DIMENSIONAL SUMMARIES
# (NB: CENTERING AT SPECIFIC TEMPERATURE VALUE)
# (NB: ALL THE ESTIMATES ARE ALREADY REPORTED BY crosspred ALONE)

cp <- crosspred(cb,model,from=bound[1],to=bound[2],by=1,cen=17)
crall <- crossreduce(cb,model,from=bound[1],to=bound[2],by=0.2,cen=17)
crlag <- crossreduce(cb,model,type="lag",value=4,from=bound[1],to=bound[2],
  bylag=0.2,cen=17)
crvar <- crossreduce(cb,model,type="var",value=22,from=bound[1],to=bound[2],
  bylag=0.2,cen=17)

# PLOTS

pdf("figure1.pdf",height=6,width=8.5)
par(mar=c(1.5,1,0,0)+0.1,cex.axis=0.9,cex.lab=1)
layout(matrix(rep(1:4,each=2),2,4,byrow=TRUE))

# 3D PLOT WITH DIFFERENT NON-DEFAULT PERSPECTIVE AND GREY SURFACE LINES
d3 <- plot(cp,xlab="Temperature (C)",zlab="RR",phi=35,theta=205,ltheta=170,
  shade=0.4)

# LINES IN THE SURFACE CORRESPONDING TO THE EFFECTS IN THE PLOTS BELOW
lines(trans3d(x=17,y=0:21,z=cp$matRRfit[as.character(17),],
  pmat=d3),lwd=2)
lines(trans3d(x=22,y=0:21,z=cp$matRRfit[as.character(22),],
  pmat=d3),lwd=2,col=2)
lines(trans3d(x=cp$predvar,y=4,z=cp$matRRfit[,"lag4"],
  pmat=d3),lwd=2,col=2)

par(mar=c(5,4,1,1)+0.1,mgp=c(2.5,1,0))

# PLOTS FOR PREDICTOR-SPECIFIC, LAG-SPECIFIC AND OVERALL CUMULATIVE SUMMARIES
plot(crvar,xlab="Lag",ylab="RR",col=2,lwd=2)
mtext(text=paste("Predictor-specific association at temperature ",22,
  "C",sep=""),cex=0.7)
plot(crlag,xlab="Temperature (C)",ylab="RR",col=2,ylim=c(.96,1.06),lwd=2)
mtext(text="Lag-specific association at lag 4",cex=0.7)
plot(crall,xlab="Temperature (C)",ylab="RR",ylim=c(.8,2),col=2,lwd=2)
mtext(text="Overall cumulative association",cex=0.7)

dev.off()

#