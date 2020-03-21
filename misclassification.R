# Chien example - bias analysis for misclassification 
## - Based on Greenland's approach (PMID: 3043623)
## - Tom Ahern (02tahern@med.uvm.edu)
## - 2020-02-18

# These computations can also be completed with an Excel spreadsheet at:
## https://sites.google.com/site/biasanalysis/exposure-misclassification-variance-correction

rm(list=ls())
cat("\014") #clears console

# assign values to crude 2x2 table
a <- 118 #cases, exposed
b <- 832 #cases, unexposed
c <- 103 #controls, exposed
d <- 884 #controls, unexposed

# crude association
or_crude <- ((a/c)/(b/d))
se_crude <- sqrt(1/a+1/b+1/c+1/d)
lcl95_crude <- exp(log(or_crude) - (1.96*se_crude))
ucl95_crude <- exp(log(or_crude) + (1.96*se_crude))

# internal validation data from Boudreau et al (PMID: 14742292)
## among cases
av1 <- 24 #classified+/validated+
bv1 <- 2 #classified+/validated-
cv1 <- 19 #classified-/validated+
dv1 <- 144 #classified-/validated-

  # terms for Greenland formulas
  P1 <- av1/(av1+cv1)
  VarP1 <- P1*(1-P1)/(av1+cv1)
  R1 <- dv1/(bv1+dv1)
  VarR1 <- R1*(1-R1)/(bv1+dv1)
  D1 <- P1+R1-1
  M1 <- a+b
  A1 <- (M1*R1-b)/D1
  B1 <- M1-A1
  Estar1 <- a/M1
  E1 <- A1/M1
  Fstar1 <- b/M1
  F1 <- B1/M1
  C1 <- Estar1*Fstar1/(M1*E1^2*F1^2)

## among controls
av0 <- 18 #classified+/validated+
bv0 <- 4 #classified+/validated-
cv0 <- 13 #classified-/validated+
dv0 <- 130 #classified-/validated-

  # terms for Greenland formulas
  P0 <- av0/(av0+cv0)
  VarP0 <- P0*(1-P0)/(av0+cv0)
  R0 <- dv0/(bv0+dv0)
  VarR0 <- R0*(1-R0)/(bv0+dv0)
  D0 <- P0+R0-1
  M0 <- c + d
  A0 <- (M0*R0-d)/D0
  B0 <- M0-A0
  Estar0 <- c/M0
  E0 <- A0/M0
  Fstar0 <- d/M0
  F0 <- B0/M0
  C0 <- Estar0*Fstar0/(M0*E0^2*F0^2)

# Greenland's variance (formula 1 in PMID 3043623)
Var_OR_g <- ((VarP1/F1^2 + VarR1/E1^2 + C1) / D1^2) + ((VarP0/F0^2 + VarR0/E0^2 + C0) / D0^2)
SE_OR_g <- sqrt(Var_OR_g)

# corrected cell frequencies
ac <- A1
bc <- B1
cc <- A0
dc <- B0

# bias-adjusted association
OR_g <- round(((ac/cc)/(bc/dc)),2)
lcl95_g <- round(exp(log(OR_g)-1.96*SE_OR_g),2)
ucl95_g <- round(exp(log(OR_g)+1.96*SE_OR_g),2)

print(paste0(OR_g, " (", lcl95_g, ", ", ucl95_g,")"))