# DiForti - probabilistic bias analysis for selection bias
# - Tom Ahern (02tahern@med.uvm.edu)
# - 2020-02-19

rm(list=ls())
cat("\014") # clears console

# - probabilistic bias analysis 
# - use log-normal distribution for selection odds ratio

# DiForti et. al. observed data
a <- 266
b <- 635
c <- 84
d <- 1153

or_crude <- (a/c)/(b/d)
se_crude <- sqrt(1/a+1/b+1/c+1/d)

# user-defined selection proportions in joint exposure (x) and outcome (y) categories
x1y1 <- 0.58
x1y0 <- 0.18
x0y1 <- 0.60
x0y0 <- 0.31

# calculate the bias factor
or_sel <- (x1y1/x1y0)/(x0y1/x0y0)
sd <- 0.21 # adopting value used by DiForti, but can be modified

# calculate location and shape to use with rlnorm
# https://msalganik.wordpress.com/2017/01/21/making-sense-of-the-rlnorm-function-in-r/
location <- log(or_sel^2 / sqrt(sd^2 + or_sel^2))
shape <- sqrt(log(1 + (sd^2 / or_sel^2)))

# draw or.sel values from the distribution, display, and summarize
set.seed(2718)
niter <- 100000
or_sel_draw <- rlnorm(niter, location, shape)
hist(or_sel_draw, breaks=100, xlim=c(0.8, 3), main="Distribution of drawn selection ORs")
summary(or_sel_draw)

# calculate vector of bias-adjusted odds ratios
or_bias <- or_crude/or_sel_draw
hist(or_bias, breaks=100, xlim=c(1.5,7), main="Distribution of selection bias-adjusted ORs")
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