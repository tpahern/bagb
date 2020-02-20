# Simple bias analysis for unmeasured confounder - Mørch et. al.

rm(list=ls())
cat("\014") # clears console

p1=.0774 #prev of family history in US population
rr.eu <- seq(0.5, 2, length=40) #vary association between hormone contraception use and famhx
p0=p1/rr.eu #estimate family history prevalence in nonusers of hormone contraception
rr.ud <- 2 #association between family history and breast cancer
rr.obs=1.2 #crude hormone contraception/breast cancer association
bf.s=(rr.eu+(rr.ud-1)*rr.eu*p1)/(rr.eu+(rr.ud-1)*p1) # bias factor (Lin formula)
rr.ba <- rr.obs/bf.s

data <- data.frame(cbind(rr.eu, rr.ba))


plot(data$rr.eu, data$rr.ba, type="l", 
     xlab="Contraception/family history prevalence ratio", 
     ylab="Bias-adjusted HR",
     cex.lab=1.3,
     cex.axis=1.2)
abline(h=1.2, col="gray", lwd=3)
lines(data$rr.eu, data$rr.ba, col="black", lwd=5)