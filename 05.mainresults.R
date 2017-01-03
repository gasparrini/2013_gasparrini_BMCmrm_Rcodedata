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
# OVERALL CUMULATIVE SUMMARY ASSOCIATION

# PLOT
pdf("figure2.pdf",height=5,width=13)

par(mar=c(5,4,1,1)+0.1,cex.axis=0.9,mgp=c(2.5,1,0))
layout(matrix(1:2,ncol=2))

plot(cpall,type="n",ylab="RR",ylim=c(.8,2),xlab="Temperature (C)")
for(i in seq(regall)) lines(regall[[i]],ptype="overall",col=grey(0.5),lty=2)
abline(h=1)
lines(cpall,col=2,lwd=2)
mtext("Main model: first-stage and pooled estimates",cex=1)
legend ("top",c("Pooled (with 95%CI)","First-stage region-specific"),
  lty=c(1,2),lwd=1.5,col=c(2,grey(0.7)),bty="n",inset=0.1,cex=0.8)

plot(cpall,ylab="RR",col=2,lwd=2,ylim=c(.8,2),xlab="Temperature (C)")
lines(cpall2,col=3,lty=2,lwd=2)
lines(cpall3,col=4,lty=4,lwd=2)
mtext("Comparison of alternative models",cex=1)
legend ("top",c("B-spline of lag 0-21 (with 95%CI)","Constant of lag 0-3",
  "Constant of lag 0-21"),lty=c(1,2,4),lwd=1.5,col=2:4,bty="n",inset=0.05,
  cex=0.8,title="Function for the lag space:")

dev.off()

# POINT OF MINIMUM MORTALITY
cpall$predvar[which.min(cpall$allRRfit)]
round(sum(regEngWales$tmean<17.1)/nrow(regEngWales)*100,1)

# Q TEST AND I-SQUARE
(qall <- qtest(mvall))
round(((qall$Q-qall$df)/qall$Q)[1]*100,1)
(qall2 <- qtest(mvall2))
round(((qall2$Q-qall2$df)/qall2$Q)[1]*100,1)
(qall3 <- qtest(mvall3))
round(((qall3$Q-qall3$df)/qall3$Q)[1]*100,1)

####################################################################
# PREDICTOR-SPECIFIC SUMMARIES

# PLOT
pdf("figure3.pdf",height=5,width=13)

par(mar=c(5,4,1,1)+0.1,cex.axis=0.9,mgp=c(2.5,1,0))
layout(matrix(1:2,ncol=2))

plot(cphot,type="n",ylab="RR",ylim=c(.95,1.12),xlab="Lag")
for(i in seq(reghot)) lines(reghot[[i]],ptype="overall",col=grey(0.5),lty=2)
abline(h=1)
lines(cphot,col=2,lwd=2)
legend ("top",c("Pooled (with 95%CI)","First-stage region-specific"),
  lty=c(1,2),lwd=1.5,col=c(2,grey(0.7)),bty="n",inset=0.1,cex=0.8)
mtext(text=paste("Predictor-specific summary for temperature = ",22,
  "C",sep=""),cex=1)

plot(cpcold,type="n",ylab="RR",ylim=c(.95,1.12),xlab="Lag")
for(i in seq(regcold)) lines(regcold[[i]],ptype="overall",col=grey(0.5),lty=2)
abline(h=1)
lines(cpcold,col=2,lwd=2)
legend ("top",c("Pooled (with 95%CI)","First-stage region-specific"),
  lty=c(1,2),lwd=1.5,col=c(2,grey(0.7)),bty="n",inset=0.1,cex=0.8)
mtext(text=paste("Predictor-specific summary for temperature = ",0,
  "C",sep=""),cex=1)

dev.off()

# OVERALL EFFECTS AT THESE TWO PREDICTOR LEVELS
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["22",]),3)
round(with(cpall,cbind(allRRfit,allRRlow,allRRhigh)["0",]),3)

# TESTS AND STATISTICS
(qhot <- qtest(mvhot))
(qcold <- qtest(mvcold))
round(((qhot$Q-qhot$df)/qhot$Q)[1]*100,1)
round(((qcold$Q-qcold$df)/qcold$Q)[1]*100,1)

####################################################################
# COMPARISON OF RANDOM VS. FIXED EFFECT MODEL FOR SUMMARY AT 22C

mvhot2 <- update(mvhot,method="fixed")
coef(mvhot); coef(mvhot2)

cphot2 <- crosspred(blag,coef=coef(mvhot2),vcov=vcov(mvhot2),
  model.link="log",at=0:210/10)

plot(cphot,ci="lines",ylab="RR",ylim=c(.95,1.12),xlab="Lag",col=2)
lines(cphot2,ci="lines",col=4)
mtext(text=paste("Predictor-specific summary for temperature = ",22,
  "C",sep=""),cex=1)
legend ("top",c("Random-effects model","Fixed-effects model"),
  lty=1,col=c(2,4),bty="n",inset=0.1,cex=0.8)

  
#