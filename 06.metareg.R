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

# INPUT THE META-VARIABLE: LATITUDE
lat <- c(54.84815,53.58832,53.72352,52.85539,52.53304,52.03734,51.50583,
  51.24213,51.05361,52.02615)

# MULTIVARIATE META-REGRESSION
(mvalllat <- update(mvall,.~lat))
summary(mvalllat)
(mvhotlat <- update(mvhot,.~lat))
summary(mvhotlat)
# NOTE THE PROBLEM FOR ESTIMATED (CO)VARIANCE MATRIX
(mvcoldlat <- update(mvcold,.~lat))
summary(mvcoldlat)

####################################################################
# PREDICTION FROM META-REGRESSION

val <- round(quantile(lat,c(25,75)/100),1)

predall <- predict(mvalllat,data.frame(lat=val),vcov=T)
cpalllat25 <- crosspred(bvar,coef=predall[[1]]$fit,vcov=predall[[1]]$vcov,
  model.link="log",by=0.2,cen=17)
cpalllat75 <- crosspred(bvar,coef=predall[[2]]$fit,vcov=predall[[2]]$vcov,
  model.link="log",by=0.2,cen=17)

predhot <- predict(mvhotlat,data.frame(lat=val),vcov=T)
cphotlat25 <- crosspred(blag,coef=predhot[[1]]$fit,vcov=predhot[[1]]$vcov,
  model.link="log",by=0.2)
cphotlat75 <- crosspred(blag,coef=predhot[[2]]$fit,vcov=predhot[[2]]$vcov,
  model.link="log",by=0.2)

predcold <- predict(mvcoldlat,data.frame(lat=val),vcov=T)
cpcoldlat25 <- crosspred(blag,coef=predcold[[1]]$fit,vcov=predcold[[1]]$vcov,
  model.link="log",by=0.2)
cpcoldlat75 <- crosspred(blag,coef=predcold[[2]]$fit,vcov=predcold[[2]]$vcov,
  model.link="log",by=0.2)

####################################################################
# RESULTS FROM META-REGRESSION

# Q TEST AND I-SQUARE
(qalllat <- qtest(mvalllat))
round(((qalllat$Q-qalllat$df)/qalllat$Q)[1]*100,1)
(qhotlat <- qtest(mvhotlat))
round(((qhotlat$Q-qhotlat$df)/qhotlat$Q)[1]*100,1)
(qcoldlat <- qtest(mvcoldlat))
round(((qcoldlat$Q-qcoldlat$df)/qcoldlat$Q)[1]*100,1)

# WALD TEST
fwald <- function(model,var) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}
round(fwald(mvalllat,"lat"),3)
round(fwald(mvhotlat,"lat"),3)
round(fwald(mvcoldlat,"lat"),3)

# OVERALL EFFECTS AT THESE TWO PREDICTOR LEVELS
round(with(cpalllat25,cbind(allRRfit,allRRlow,allRRhigh)["22",]),3)
round(with(cpalllat75,cbind(allRRfit,allRRlow,allRRhigh)["22",]),3)
round(with(cpalllat25,cbind(allRRfit,allRRlow,allRRhigh)["0",]),3)
round(with(cpalllat75,cbind(allRRfit,allRRlow,allRRhigh)["0",]),3)

# PLOT 
pdf("figure4.pdf",height=6,width=8.5)
par(mar=c(5,4,1,1)+0.1,cex.axis=0.9,mgp=c(2.5,1,0))
layout(matrix(c(0,1,1,0,2,2,3,3),2,4,byrow=TRUE))

plot(cpalllat75,type="n",ci="n",ylab="RR",ylim=c(.8,2),xlab="Temperature (C)")
lines(cpalllat25,col=2,lty=4,lwd=2,ci="area",
  ci.arg=list(density=20,col=2))
lines(cpalllat75,col=4,lty=5,lwd=2,ci="area",
  ci.arg=list(density=20,angle=-45,col=4))
legend("top",paste(val),lty=4:5,col=c(2,4),lwd=1.5,bty="n",
  title="Latitude (degrees North)",inset=0.1,cex=0.9)
mtext("Overall cumulative summary",cex=0.7)

plot(cphotlat75,type="n",ci="n",ylab="RR",ylim=c(.95,1.12),xlab="Lag")
lines(cphotlat25,col=2,lty=4,lwd=2,ci="area",
  ci.arg=list(density=20,col=2))
lines(cphotlat75,col=4,lty=5,lwd=2,ci="area",
  ci.arg=list(density=20,angle=-45,col=4))
legend("top",paste(val),lty=4:5,col=c(2,4),lwd=1.5,bty="n",
  title="Latitude (degrees North)",inset=0.1,cex=0.9)
mtext(text=paste("Predictor-specific summary for temperature = ",22,
  "C",sep=""),cex=0.7)

plot(cpcoldlat75,type="n",ci="n",ylab="RR",ylim=c(.95,1.12),xlab="Lag")
lines(cpcoldlat25,col=2,lty=4,lwd=2,ci="area",
  ci.arg=list(density=20,col=2))
lines(cpcoldlat75,col=4,lty=5,lwd=2,ci="area",
  ci.arg=list(density=20,angle=-45,col=4))
legend("top",paste(val),lty=4:5,col=c(2,4),lwd=1.5,bty="n",
  title="Latitude (degrees North)",inset=0.1,cex=0.9)
mtext(text=paste("Predictor-specific summary for temperature = ",0,
  "C",sep=""),cex=0.7)

dev.off()


