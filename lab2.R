#question1
library(readxl)
library(mvtnorm)

data = read_xlsx(path = 'C:/Users/Varun/Desktop/LiU_STIMA/Semester_2/Period_4/Bayesian Learning/Lab/Lab2/Linkoping2022.xlsx')


data_df = as.data.frame(data)

date = as.Date(data_df[,2])

time = as.numeric(format(date, "%j"))/365
temp = data_df[,3]

resp_cov_df = data.frame(covariate = time, response = temp)

#Step1: draws from chi squared distribution
n_draws <- 20
nu_0 <- 3
# n <- nrow(resp_cov_df)
X <- rchisq(n_draws, df = nu_0)

#Step2: Compute sigma^2 - draw from inverse chi-squared
sigma_0 <- 12
nu_0 <- 3
inv_chi <- (nu_0 * sigma_0)/X

omega_0 <- diag(x = 4.5, nrow = 3)
mu_0 <- t(matrix(data = c(-7.532457, 83.578936, -78.298468), nrow = 1))
# rmvnorm(n = 1, mean = mu_0, sigma = inv_chi[1]*(omega_0))
beta_prior <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(beta_prior) <- c('beta_0', 'beta_1', 'beta_2')

set.seed(1234)
for (i in 1:n_draws) {
  beta_prior <- rbind(beta_prior, rmvnorm(n = 1, mean = mu_0, sigma = inv_chi[i]*solve(omega_0)))
}
colnames(beta_prior) <- c('beta_0', 'beta_1', 'beta_2')


pred_df <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(pred_df) <- c('beta_0', 'beta_1', 'beta_2')

sim_temp <- list()
for (j in 1:nrow(beta_prior)) {
  sim_temp[[j]] <- beta_prior[j,1] + beta_prior[j,2]*resp_cov_df$covariate + beta_prior[j,3]*(resp_cov_df$covariate^2)
}

pred_df <- rbind(pred_df, beta_prior)
pred_df <- cbind(pred_df, I(sim_temp))

plot(resp_cov_df$covariate, resp_cov_df$response, type = 'p', col = 'black', ylim = c(-10, 30))

colors <- rainbow(nrow(pred_df))
for (j in 1:nrow(pred_df)) {
  lines(resp_cov_df$covariate, pred_df$sim_temp[[j]], col = colors[j])
}
# how to select the optimal prior params

# get_params <- function(x){
#   nu_0 = x[1]
#   sigma_0 = x[2]
#   omega_0 = x[3]
#   mu_0 = x[4:6]
#   n_draws <- 20
#   X <- rchisq(20, df = 1)
#   inv_chi <- (nu_0 * sigma_0)/X
#   precision_mat <- diag(x = omega_0, nrow = 3)
#   mu_mat <- t(matrix(data = c(mu_0), nrow = 1))
#   beta_prior <- data.frame(matrix(nrow = 0, ncol = 3))
#   colnames(beta_prior) <- c('beta_0', 'beta_1', 'beta_2')
# 
#   for (i in 1:n_draws) {
#     beta_prior <- rbind(beta_prior, rmvnorm(n = 1, mean = mu_mat, sigma = inv_chi[i]*solve(precision_mat)))
#   }
#   colnames(beta_prior) <- c('beta_0', 'beta_1', 'beta_2')
# 
# 
#   pred_df <- data.frame(matrix(nrow = 0, ncol = 3))
#   colnames(pred_df) <- c('beta_0', 'beta_1', 'beta_2')
# 
#   sim_temp <- list()
#   min_err <- c()
#   for (j in 1:nrow(beta_prior)) {
#     sim_temp[[j]] <- beta_prior[j,1] + beta_prior[j,2]*resp_cov_df$covariate + beta_prior[j,3]*(resp_cov_df$covariate^2)
#     min_err <- mean((resp_cov_df[,2]-sim_temp[[j]])^2)
#   }
#   
#   # pred_df <- rbind(pred_df, beta_prior)
#   # pred_df <- cbind(pred_df, I(sim_temp))
#   return(min_err)
# 
# }
# # pred_vals <- get_params(nu_0 = 3, omega_0 = 0.9, sigma_0 = 3, mu_0 = c(0, 75, -75))
# # # 
# # # 
# # # optimal_params <- optim(par = c(1, 1, 1, 1), fn = get_params, gr = 'BFGS')
# # # optimal_params <- optim(par = c(nu_0 = 3, sigma_0 = 3, omega_0 = 0.9, mu_0 = c(0, 75, -75)),fn = get_params, gr = 'BFGS')
# # # 
# optimal_params <- optim(par = c(1, 1, 1, 1, 100, 100), fn = get_params, method = 'BFGS')



# 1b i.
library(LaplacesDemon)
library(ggplot2)
covariates <- cbind(1, resp_cov_df$covariate, resp_cov_df$covariate**2)
response <- resp_cov_df$response
dim(ols_estimate)
ols_estimate <- solve(t(covariates) %*% covariates) %*% t(covariates) %*% response

mu_n <- solve(t(covariates)%*%covariates+omega_0) %*% (t(covariates) %*% covariates %*% ols_estimate + omega_0 %*% mu_0)
omega_n <- t(covariates) %*% covariates + omega_0
nu_n <- nu_0 + nrow(covariates)

sigma_n <- (nu_0*sigma_0 + (t(response) %*% response + t(mu_0) %*% omega_0 %*% mu_0 - t(mu_n) %*% omega_n %*% mu_n))/(nu_n)

beta_posterior <- data.frame(matrix(nrow = 0, ncol = 3))

sigma_posterior <- rinvchisq(n = 1000, df = nu_n, scale = sigma_n)

# for (i in 1:1000) {
#   beta_posterior <- rbind(beta_posterior, rmvnorm(n = 1, mean = mu_n, sigma = sigma_posterior[i]*solve(omega_n)))
# }
# colnames(beta_posterior) <- c('beta_0', 'beta_1', 'beta_2')
n <- nrow(covariates)
k <- ncol(covariates)
# s2 <- (1/(n-k))*t(response - covariates %*% ols_estimate) %*% (response - covariates %*% ols_estimate)
for (i in 1:length(sigma_posterior)) {
  beta_posterior <- data.frame(mvtnorm::rmvnorm(n = 1000, mean = mu_n, sigma = sigma_posterior[i]*solve(omega_n)))
}

colnames(beta_posterior) <- c('beta_0', 'beta_1', 'beta_2')


ggplot(beta_posterior, aes(beta_0)) +
  geom_histogram(aes(y=..density..), bins = 50, color="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

ggplot(beta_posterior, aes(beta_1)) +
  geom_histogram(aes(y=..density..), bins = 50, color="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

ggplot(beta_posterior, aes(beta_2)) +
  geom_histogram(aes(y=..density..), bins = 50, color="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

sigma_posterior <- data.frame(sigma_posterior)
colnames(sigma_posterior) <- 'sigma_posterior'
ggplot(data.frame(sigma_posterior), aes(sigma_posterior)) +
  geom_histogram(aes(y=..density..), bins = 50, color="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# 1b ii.

covariates <- as.matrix(covariates)
beta_posterior <- as.matrix(beta_posterior)
# posterior_median <- c(median_beta0, median_beta1, median_beta2)

posterior_temp_median <- beta_posterior %*% t(covariates)
posterior_temp_median_new <- apply(posterior_temp_median, 2, median)

plot(resp_cov_df$covariate, resp_cov_df$response)
lines(resp_cov_df$covariate,posterior_temp_median_new, col = 'red')

# posterior_temp <- as.matrix(beta_posterior) %*% t(covariates)

post_prob_interval <- apply(posterior_temp_median, 2, quantile, probs = c(0.05, 0.95))
post_prob_interval <- t(post_prob_interval)

lines(resp_cov_df$covariate, post_prob_interval[,1], lty = 2, col = 'blue')
lines(resp_cov_df$covariate, post_prob_interval[,2], lty = 2, col = 'blue')

calc_x <- function(beta) {
  return(-beta[2] / (2 * beta[3]))
}

xTilde <- apply(beta_posterior, 1, calc_x)
hist(xTilde, breaks = 25)

# 2 a).

data_2 <- read.table("WomenAtWork.dat", header = TRUE)

covs <- c(2:ncol(data_2))
X <- as.matrix(data_2[,covs])
npar <- dim(X)[2]

y <- as.matrix(data_2[,1])

mu <- as.matrix(rep(0,npar))
tau <- 2
Sigma <- (tau^2)*diag(npar)

LogPostLogistic <- function(betas,y,X,mu,Sigma){
  linPred <- X%*%betas;
  logLik <- sum( linPred*y - log(1 + exp(linPred)) );
  #if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  logPrior <- dmvnorm(betas, mu, Sigma, log=TRUE);
  
  return(logLik + logPrior)
}

initVal <- matrix(0,npar,1)

OptimRes <- optim(initVal,LogPostLogistic,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)
beta_posteriorq2 <- OptimRes$par
hessian_post <- OptimRes$hessian
# approxPostStd <- sqrt(diag(solve(-OptimRes$hessian)))

posterior_dist_param <- data.frame(mvtnorm::rmvnorm(n = 1000, mean = beta_posteriorq2, sigma = solve(-hessian_post)))

NSmallChild_posterior <- matrix(posterior_dist_param[,6])

post_prob_interval_NSmallChild <- apply(NSmallChild_posterior, 2, quantile, probs = c(0.05, 0.95))
hist(NSmallChild_posterior)
abline(v=post_prob_interval_NSmallChild[1,1], lty = 2, col = 'red')
abline(v=post_prob_interval_NSmallChild[2,1], lty = 2, col = 'red')

glmModel<- glm(Work ~ 0 + ., data = data_2, family = binomial)


# b).

new_xVals <- matrix(c(1, 18, 11, 7, 40, 1, 1), nrow = 1)
as.matrix(posterior_dist_param)
post_linpred <- as.matrix(posterior_dist_param) %*% t(new_xVals)
logit_pred <- matrix(exp(post_linpred)/(1+exp(post_linpred)))
hist(logit_pred)
finalPred <- ifelse(logit_pred<0.5, 0, 1)

posterior_pred_draw <- function(X, beta){
  post_linpred <- as.matrix(beta) %*% t(X)
  logit_pred <- matrix(exp(post_linpred)/(1+exp(post_linpred)))
  # finalPred <- ifelse(logit_pred<0.5, 0, 1)
  return(1-logit_pred)
}

posterior_predictions <- posterior_pred_draw(new_xVals, posterior_dist_param)
hist(posterior_predictions)

# hist(sapply(posterior_predictions,function(x) rbinom(1, size = 13, prob = x)))

# c). 

posterior_pred_draw_binom <- function(prob,n){
  # post_linpred <- as.matrix(beta) %*% t(X)
  # logit_pred <- matrix(exp(post_linpred)/(1+exp(post_linpred)))
  # prob_not_work <- 1-logit_pred
  num_not_work <- rbinom(n = length(prob), size = n, prob = prob)
  return(num_not_work)
}
hist(posterior_pred_draw_binom(posterior_predictions ,n = 13))
