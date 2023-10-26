# Q1
n <- 70
s <- 22
f <- n-s
prob <- s/n
sd_bern <- prob*(1-prob)


posteriorBetaPlot <- function(a0, b0){
  n <- 70
  s <- 22
  f <- n-s
  mean_post <- c()
  sd_post <- c()
  for (i in 1:10000) {
    posterior_vals <- rbeta(n = i, shape1 = a0+s, shape2 = b0+f)
    mean_post <- c(mean_post, mean(posterior_vals))
    sd_post <- c(sd_post, sd(posterior_vals))
  }
  return(list(mean_post, sd_post))
}
a0 <- 8
b0 <- 8

posterior <- posteriorBetaPlot(a0 = 8, b0 = 8)
plot(posterior[[1]], type = 'l')
abline(h=(a0+s)/(a0+s+b0+f), lty=2, col='red')
plot(posterior[[2]], type = 'l')
abline(h=sqrt(((a0+s)*(b0+f))/(((a0+s+b0+f)^2)*(a0+s+b0+f+1))),lty=2, col='red')

draw_posterior <- rbeta(10000,a0+s,b0+f)
draw_posterior_prob <- draw_posterior[draw_posterior>0.3]
posterior_prob <- length(draw_posterior_prob)/length(draw_posterior)
exact_prob <- pbeta(q = 0.3, shape1 = a0+s, shape2 = b0+f, lower.tail = FALSE)

posterior_odds <- c()
for (i in 1:length(draw_posterior)) {
  odds_calc <- draw_posterior[i]/(1-draw_posterior[i])
  posterior_odds <- c(posterior_odds, odds_calc)
}
hist(posterior_odds)
mean(posterior_odds)
# Q2

y <- c(33, 24, 48, 32, 55, 74, 23, 17)
hist(y, breaks = 4)
mean(y)
draw_x <- rchisq(n = 10000, df = length(y))
hist(draw_x)


tau_sq <- sum((log(y)-3.6)^2)/length(y)  

tau_sq

si=((length(y))*tau_sq)/draw_x
mean(si)
hist(si, xlim = c(0,1), breaks = 100)




# theta <- c()
# sig_sq <- c()
# for (i in 1:length(draw_x)) {
#   sig_sq <- c(sig_sq, ((length(y)-1)*tau_sq)/draw_x[i])
#   theta <- c(theta,rnorm(1, mean = mean(y),sd = sqrt(sig_sq[i]/length(y))))
# }
# hist(sig_sq)
# hist(theta)
# mean(theta)

z <- sqrt(si)/sqrt(2)

gini <- 2*pnorm(q = z, mean = 0, sd = 1) - 1
mean(gini)
hist(gini, freq = FALSE, breaks = 100)
abline(v = quantile(x = gini, probs = c(0.025,0.975))[1], lty = 2, col = 'red')
abline(v = quantile(x = gini, probs = c(0.025,0.975))[2], lty = 2, col = 'red')

dens <- density(x = gini)
max(dens$y)

hdpi <- hdi(dens, credMass = 0.95)
hdpi
abline(v=hdpi[1], col='blue')
abline(v=hdpi[2], col = 'blue')
print(hdpi)

plot(dens)

f <- approxfun(dens$x, dens$y)
tot_area <- integrate(f, min(dens$x), max(dens$x))

mid_ind <- order(dens$y, decreasing = TRUE)[1]
i <- mid_ind
j <- mid_ind

curv_area <- c()
breadth <- abs(dens$x[i] - dens$x[i-1])
len <- dens$y[i]
curv_area <- c(curv_area, len*breadth/tot_area$value)

while (sum(curv_area) <= 0.95) {
  if (dens$y[i] > 0.4) {
    i <- i-1
    breadth <- abs(dens$x[i] - dens$x[i-1])
    len <- dens$y[i]
    curv_area <- c(curv_area, len*breadth/tot_area$value)
  }
  if (j < length(dens$y)) {
    j <- j+1
    breadth <- abs(dens$x[j] - dens$x[j-1])
    len <- dens$y[j]
    curv_area <- c(curv_area, len*breadth/tot_area$value)
  }
}
i <- i-1
j <- j+1
plot(dens)
print(paste("Lower:", dens$x[i]))
print(paste("Upper:", dens$x[j]))
abline(v = c(dens$x[i], dens$x[j]), col= 'green')

# Question3

k = seq(0.001, 10, 0.1)

posterior <- function(k){
  y <- c(-2.79, 2.33, 1.83, -2.44, 2.23, 2.33, 2.07, 2.02, 2.14, 2.54)
  mu <- 2.4
  return(exp((k*sum(cos(y-mu)))-(k/2))/(besselI(k, nu = 0))^length(y))
}
normaliz_const <- integrate(posterior, min(k), max(k))
norm_pdf <- posterior(k)/normaliz_const$value

f <- approxfun(k, norm_pdf)
tot_area <- integrate(f, min(k), max(k))

plot(norm_pdf)
posterior_mode <- k[which.max(norm_pdf)]
