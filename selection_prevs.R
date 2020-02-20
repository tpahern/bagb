# DiForti - probabilistic bias analysis for selection bias
# - probabilistic correction based on estimated exposure pravelences in nonparticipating cases and controls
# - Tom Ahern (02tahern@med.uvm.edu)
# - 2020-02-19

library(tidyverse)
library(reshape2)
library(RColorBrewer)

rm(list=ls())
cat("\014") # clears console

# DiForti et. al. observed data
a <- 266 #exposed cases
b <- 635 #unexposed cases
c <- 84 #exposed controls
d <- 1153 #unexposed controls

n_cases <- a+b
n_ctrls <- c+d
n_exp <- a+c
n_unexp <- b+d

# DiForti et. al. site exclusions
ax <- 67.6
bx <- 161.4
cx <- 17.8
dx <- 244.2

# DiForti et. al. refusals
ar <- 126.3
br <- 262.7
cr <- 361.9
dr <- 2303.0

# DiForti et. al. all non-participants
anp <- ax+ar
bnp <- bx+br
cnp <- cx+cr
dnp <- dx+dr

# total nonparticipants by case/control status
np_cases <- round(anp+bnp, 0)
np_ctrls <- round(cnp+dnp, 0)

# crude association (observed data)
or_crude <- (a/c)/(b/d)
se_crude <- sqrt(1/a+1/b+1/c+1/d)
lcl95_crude <- exp(log(or_crude) - (1.96*se_crude))
ucl95_crude <- exp(log(or_crude) + (1.96*se_crude))

# calculate exposure prevalences in nonparticipating cases (p1) and controls (p0)
p1 <- anp/np_cases
p0 <- cnp/np_ctrls

# set some important parameters
set.seed(33)
niter <- 1e5

# sample prevalences from beta distributions
p1_draw <- rbeta(niter, anp, bnp)
p0_draw <- rbeta(niter, cnp, dnp)

prevdraws <- data.frame(p1=p1_draw, p0=p0_draw) %>% 
  melt()

# plot overlaid densities of drawn p1, p0 values
ggplot(prevdraws, 
       aes(x=value, fill=variable)) + 
  geom_density(alpha=0.25) + 
  theme_minimal() + 
  scale_fill_brewer(palette="Set1")

# vector of random standard normal deviates for simulation intervals
rannorm <- rnorm(niter)

#compute crude OR using simulation
or_crude_c <- exp(log(or_crude) - (rannorm*se_crude))

# calculated bias-adjusted cell frequencies (first source of uncertainty)
n_exp_cases <- a+p1_draw*np_cases
n_unexp_cases <- b+(1-p1_draw)*np_cases
n_exp_ctrls <- c+p0_draw*np_ctrls
n_unexp_ctrls <- d+(1-p0_draw)*np_ctrls

# calculate OR and its standard error (first source of uncertainty)
or_ba <- ((n_exp_cases/n_exp_ctrls)/(n_unexp_cases/n_unexp_ctrls))
se_or_ba <- sqrt(1/n_exp_cases+1/n_exp_ctrls+1/n_unexp_cases+1/n_unexp_ctrls)
## add random error uncertainty from conventional SE and bias adjusted SE
or_bac <- exp(log(or_ba)-(rannorm*se_crude))
or_baca <- exp(log(or_ba)-(rannorm*se_or_ba))

# calculate bias-adjusted OR (second soure of uncertainty)
## calculate bias-adjusted cell frequencies
n_exp_cases_b <- a+rbinom(niter, np_cases, p1_draw)
n_unexp_cases_b <- n_cases+np_cases-n_exp_cases_b
n_exp_ctrls_b <- c+rbinom(niter, np_ctrls, p0_draw)
n_unexp_ctrls_b <- n_ctrls+np_ctrls-n_exp_ctrls_b
## calculate OR and its standard error (second source of uncertainty)
or_ba_b <- ((n_exp_cases_b/n_exp_ctrls_b)/(n_unexp_cases_b/n_unexp_ctrls_b))
se_or_ba_b <- sqrt(1/n_exp_cases_b+1/n_exp_ctrls_b+1/n_unexp_cases_b+1/n_unexp_ctrls_b)
## add random error uncertainty from conventional SE and bias adjusted SE
or_bac_b <- exp(log(or_ba_b)-(rannorm*se_crude))
or_baca_b <- exp(log(or_ba_b)-(rannorm*se_or_ba_b))

# summarize bias analysis results
summary(n_exp_cases)
summary(n_exp_cases_b)
summary(n_exp_ctrls)
summary(n_exp_ctrls_b)

# summarize the crude odds ratio (simulation)
print("Crude association: simulated")
summary(or_crude_c)
quantile(or_crude_c, c(0.025, 0.5, 0.975), na.rm = TRUE)

# summarize the bias-adjusted odds ratios
print("Bias-adjusted (first source of uncertainty)")
summary(or_ba)
quantile(or_ba, c(0.025, 0.5, 0.975), na.rm = TRUE)

print("Bias-adjusted (first source of uncertainty): random error based on se_crude")
summary(or_bac)
quantile(or_bac, c(0.025, 0.5, 0.975), na.rm = TRUE)

print("Bias-adjusted (first source of uncertainty): random error based on se_or_ba")
summary(or_baca)
quantile(or_baca, c(0.025, 0.5, 0.975), na.rm = TRUE)

print("Bias-adjusted (second source of uncertainty)")
summary(or_ba_b)
quantile(or_ba_b, c(0.025, 0.5, 0.975), na.rm = TRUE)

print("Bias-adjusted (second source of uncertainty): randon error based on se_crude")
summary(or_bac_b)
quantile(or_bac_b, c(0.025, 0.5, 0.975), na.rm = TRUE)

print("Bias-adjusted (second source of uncertainty): randon error based on se_or_ba_b")
summary(or_baca_b)
quantile(or_baca_b, c(0.025, 0.5, 0.975), na.rm = TRUE)