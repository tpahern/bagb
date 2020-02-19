# Chien example - probabilistic bias analysis for misclassification 
# - Tom Ahern (02tahern@med.uvm.edu)
# - 2020-02-18

library(tidyverse)
library(reshape2)
library(MASS)

setwd("/Users/tomahern/Dropbox/Projects/NLM/reproducibility/Manuscripts/Bias analysis gone bad/bagb_shiny")
rm(list=ls())
cat("\014") #clears console

# assign values to crude 2x2 table
a <- 118 #cases, exposed
b <- 832 #cases, unexposed
c <- 103 #controls, exposed
d <- 884 #controls, unexposed

n_case <- a + b
n_ctrl <- c + d

# crude association
or_crude <- ((a/c)/(b/d))
se_crude <- sqrt(1/a+1/b+1/c+1/d)
lcl95_crude <- exp(log(or_crude) - (1.96*se_crude))
ucl95_crude <- exp(log(or_crude) + (1.96*se_crude))

# bias parameters derived from Boudreau (PMCID: 14742292)

se_case_alpha <- 24
se_case_beta <- 19

sp_case_alpha <- 144
sp_case_beta <- 2

se_ctrl_alpha <- 18
se_ctrl_beta <- 13

sp_ctrl_alpha <- 130
sp_ctrl_beta <- 4

# set some important parameters
niter <- 100000 #number of iterations
corr <- 0.8 #desired correlation between case/control Se, Sp draws

# set up 3 pairs of correlated random variables for beta distibutions
# (NORTA approach)
set.seed(33) #because Larry Bird AND Kareem Abdul-Jabbar
mu <- c(0,0)
sig <- matrix(c(1,corr,corr,1), nrow=2)

# pair for Se
z1 <- mvrnorm(niter, mu, sig)
u1 <- pnorm(z1)

# pair for Sp
z2 <- mvrnorm(niter, mu, sig)
u2 <- pnorm(z2)

# uncorrelated vectors for use with exposure prevalences
h1 <- runif(niter)
h2 <- runif(niter)

cor(z1) #original variable - should be close to 'corr'
cor(z2)

cor(u1) #should be a bit less than 'corr'
cor(u2)

cor(h1,h2) #no correlation

# draw sensitivities for cases and controls
se_case <- qbeta(u1[,1],se_case_alpha,se_case_beta)
se_ctrl <- qbeta(u1[,2],se_ctrl_alpha,se_ctrl_beta)

# draw specificities for cases and controls
sp_case <- qbeta(u2[,1],sp_case_alpha,sp_case_beta)
sp_ctrl <- qbeta(u2[,2],sp_ctrl_alpha,sp_ctrl_beta)

# check out distributions and correlations
summary(se_case)
summary(se_ctrl)
cor(se_case, se_ctrl)

summary(sp_case)
summary(sp_ctrl)
cor(sp_case, sp_ctrl)

# calculate bias-adjusted cell frequencies (wholesale application of Se, Sp)
ac <- (a - n_case*(1-sp_case))/(se_case - (1-sp_case)) #bias-adjusted cases, exposed 
bc <- n_case - ac #bias-adjusted cases, unexposed
cc <- (c - n_ctrl*(1-sp_ctrl))/(se_ctrl - (1-sp_ctrl)) #bias-adjusted controls, exposed
dc <- n_ctrl - cc #bias-adjusted controls, unexposed

# vector of standard normal deviates for simulation of random error
set.seed(47) 
rannorm <- rnorm(niter, mean=0, sd=1)

# calculate bias-adjusted odds ratios - only Se, Sp uncertainty
or_b <- (ac/cc)/(bc/dc)
se_or_b <- sqrt(1/ac+1/bc+1/cc+1/dc)

# add uncertainty due to bias-adjusted standard error
or_br <- exp(log(or_b) - (rannorm*se_or_b))

# calculate bias-adjusted odds ratio - individual application of Se, Sp
# 1. calculate prevalence of exposure in cases and controls, accounting for sampling error
PrevE_cases <- qbeta(h1, ac, n_case-ac)
PrevE_controls <- qbeta(h2, cc, n_ctrl-cc)
PrevE_cases_c <- ac/n_case
PrevE_controls_c <- cc/n_ctrl

# 2. calculate PPV and NPV of exposure classification in cases and controls
PPV_case <- (se_case*PrevE_cases)/((se_case*PrevE_cases)+(1-sp_case)*(1-PrevE_cases))
PPV_control <- (se_ctrl*PrevE_controls)/((se_ctrl*PrevE_controls)+(1-sp_ctrl)*(1-PrevE_controls))
NPV_case <- (sp_case*(1-PrevE_cases))/((1-se_case)*PrevE_cases+sp_case*(1-PrevE_cases))
NPV_control <- (sp_ctrl*(1-PrevE_controls))/((1-se_ctrl)*PrevE_controls+sp_ctrl*(1-PrevE_controls))

# 3. calculate the expected number of exposed cases and controls
ad <- rbinom(niter, a, PPV_case) + rbinom(niter, b, 1-NPV_case)
bd <- n_case - ad
cd <- rbinom(niter, c, PPV_control) + rbinom(niter, d, 1-NPV_control)
dd <- n_ctrl - cd

# 4. calculate bias-adjusted OR with second source of uncertainty
or_bi <- (ad/cd)/(bd/dd)
se_or_bi <- sqrt(1/ad+1/bd+1/cd+1/dd)
or_bri <- exp(log(or_bi) - rannorm*se_or_bi)

# create dataframes with variables from each bias adjustment approach
tb_whole <- as_tibble(data.frame(se_case, se_ctrl, sp_case, sp_ctrl, ac, bc, cc, dc, or_b, or_br))
tb_indiv <- as_tibble(data.frame(se_case, se_ctrl, sp_case, sp_ctrl, ad, bd, cd, dd, or_bi, or_bri,
                                 PrevE_cases, PrevE_controls, PrevE_cases_c, PrevE_controls_c,
                                 PPV_case, PPV_control, NPV_case, NPV_control))

# retain rows with positive bias-adjusted cell frequencies
chien_whole <- tb_whole %>%
  filter(ac > 0) %>%
  filter(bc > 0) %>%
  filter(cc > 0) %>%
  filter(dc > 0) 

chien_indiv <- tb_indiv %>% 
  filter(ad > 0) %>% 
  filter(bd > 0) %>% 
  filter(cd > 0) %>% 
  filter(dd > 0)

# report the number proportion of iterations discarded
n_discard_whole <- niter - length(chien_whole$or_br)
p_discard_whole <- (n_discard_whole/niter)*100
print(paste0(n_discard_whole, " iterations discarded in wholesale correction", " (", p_discard_whole, "%)"))

n_discard_indiv <- niter - length(chien_indiv$or_bri)
p_discard_indiv <- (n_discard_indiv/niter)*100
print(paste0(n_discard_indiv, " iterations discarded in individual correction", " (", p_discard_indiv, "%)"))

# plot distribution of drawn Se, Sp in retained records
draws_se <- chien_indiv[,1:2] %>% 
melt()
draws_se_plot <- ggplot(draws_se, 
                        aes(x=value, fill=variable)) + geom_density(alpha=0.25)

draws_sp <- chien_indiv[,3:4] %>%
melt()
draws_sp_plot <- ggplot(draws_sp,
                        aes(x=value, fill=variable)) + geom_density(alpha=0.25)

# summarize the wholesale bias-adjusted odds ratios
summary(chien_whole$or_br)
hist(chien_whole$or_br, breaks=10000, xlim=c(0.5,5), main="Distribution of wholesale misclassification-adjusted ORs (with random error)")
# get 2.5th, 50th, and 97.5th percentiles
cat("Wholesale correction - systematic error:", "\n")
quantile(chien_whole$or_b, probs=c(0.025, 0.5, 0.975))
cat("Wholesale correction - systematic + random error:", "\n")
quantile(chien_whole$or_br, probs=c(0.025, 0.5, 0.975))

# summarize the individualized bias-adjusted odds ratios
summary(chien_indiv$or_bri)
hist(chien_indiv$or_bri, breaks=10000, xlim=c(0.5,5), main="Distribution of individualized misclassification-adjusted ORs (with random error)")
# get 2.5th, 50th, and 97.5th percentiles
cat("Individualized correction - systematic error:", "\n")
quantile(chien_indiv$or_bi, probs=c(0.025, 0.5, 0.975))
cat("Individualized correction - systematic + random error:", "\n")
quantile(chien_indiv$or_bri, probs=c(0.025, 0.5, 0.975))