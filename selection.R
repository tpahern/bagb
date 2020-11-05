# DiForti - probabilistic bias analysis for selection bias
# - Tom Ahern (02tahern@med.uvm.edu) & Tim Lash (tlash@emory.edu)
# - 2020-11-05

rm(list=ls())
cat("\014") # clears console
library(triangle)


# DiForti et. al. observed data
a <- 266
b <- 635
c <- 84
d <- 1153

or_crude <- (a/c)/(b/d)
se_crude <- sqrt(1/a+1/b+1/c+1/d)

set.seed(2718)
niter <- 100000

#draw participation proportion ratio controls to cases (fixed 2.5 originally)
part_ratio <- rltriangle(niter,0.9,3.9,2.5)
case_p_rate=389/(389+229+901)

#total controls
control_t <- 1499*part_ratio*case_p_rate/(1-part_ratio*case_p_rate)

#draw exposure prevalence ratio for cases (fixed 1.1 originally)
case_e_ratio <- rltriangle(niter,0.23,2,1.1)

#compute number of exposed cases among non-participants
cases_e_nonp <- a/(a+b)*389*case_e_ratio

#draw exposure prevalence ratio for controls (fixed 2.0 originally)
control_e_ratio <- rltriangle(niter,0.9,8.7,2)
#compute number of exposed controls among non-participants
controls_e_nonp <- c/(c+d)*control_t*control_e_ratio

#compute selection factor
x1y1 <- a/(a+229*a/(a+b)+cases_e_nonp)
x1y0 <- c/(c+262*c/(c+d)+controls_e_nonp)
x0y1 <- b/(b+229*b/(a+b)+(389-cases_e_nonp))
x0y0 <- d/(d+262*d/(c+d)+(control_t-1499-controls_e_nonp))

# calculate the bias factor
or_sel_draw <- (x1y1/x1y0)/(x0y1/x0y0)

# calculate vector of bias-adjusted odds ratios
or_bias <- or_crude/or_sel_draw
hist(or_bias, breaks=100, xlim=c(0,7), main="Distribution of selection bias-adjusted ORs")
summary(or_bias)

# re-incorporate random error
rannorm <- rnorm(niter, mean=0, sd=1)
or_bias_rand <- exp(log(or_bias) + se_crude*rannorm)
summary(or_bias_rand)

# get 2.5th, 50th, and 97.5th percentiles
cat("Systematic error:", "\n")
quantile(or_bias, probs=c(0.025, 0.5, 0.975))

cat("Systematic + random error:", "\n")
quantile(or_bias_rand, probs=c(0.025, 0.5, 0.975))
ucl95_or_bias_rand <- quantile(or_bias_rand, 0.975)
lcl95_or_bias_rand <- quantile(or_bias_rand, 0.025)

# Now correct for bias due to confounding:
#  apply the Greenland & Mickey variance adjustment
rrc <- 1.94

# variance of the confounding-adjusted OR
v_log_oradj <- ((log(4.1)-log(2.2))/(2*1.96))^2

# variance of the crude OR (=6.2; see paper for explanation)
v_log_orcrude <- ((log(8.2)-log(4.8))/(2*1.96))^2

# calculate variance of the RRc
v_rrc <- v_log_oradj - v_log_orcrude

# variance of or_bias_rand
v_or_bias_rand <- ((log(ucl95_or_bias_rand)-log(lcl95_or_bias_rand))/(2*1.96))^2

# variance when adjusting for RRc
v_final <- v_or_bias_rand + v_rrc

# adjust point estimate for confounding
or_bias_rand_conf <- quantile(or_bias_rand, 0.5)/rrc
lcl95_or_bias_rand_conf <- exp(log(or_bias_rand_conf) - (1.96*sqrt(v_final)))
ucl95_or_bias_rand_conf <- exp(log(or_bias_rand_conf) + (1.96*sqrt(v_final)))
cat("Systematic + random error; adjusted for confounding:", "\n")
cat("OR =",round(or_bias_rand_conf,2),", 95% SI: ",round(lcl95_or_bias_rand_conf,2),",",round(ucl95_or_bias_rand_conf,2))
