# Purpose: see how reparameterization (the "Matt trick") is affected by the
# magnitude of group-level variance

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(tidyr)
library(dplyr)

# Compile all stan models
for (i in 1:4) {
  for (sfx in c("cp", "ncp")) {
    nam <- paste0("model", i, "_", sfx)
    code_file <- paste0("/Users/susanna/projects/centering_testing/src/varying_intercept_only/",
                        nam, ".stan")
    assign(nam, stan_model(file = code_file))
  }
}

J.list <- c(5, 15, 30) # number of groups
tau0.list <- 10^c(-3:1) # group-level variance
n.obs.per.grp <- 30
sigma_y <- 0.05
w <- runif(J, min = -1, max = 1)
x <- runif(n, min = -1, max = 1)
mu0 <- 0.2

for (j in 1:length(J.list)) {
  for (t in 1:length(tau0.list)) {
    
  }
}
# Draw a single dataset; create data for other 4 scenarios by setting coeffs to
# be 0 as appropriate
j <- 2
J <- J.list[j]
n <- n.obs.per.grp*J
t <- 3
tau0 <- tau0.list[t]

ess.per.sec <- data.frame(model = rep(c(1:4), each = 2),
                          param.type = rep(c("cp", "ncp"), times = 4),
                          tau0.value = rep(tau0, times = 8),
                          J = rep(J, times = 8),
                          tau0.neff = rep(0, times = 8))
stan.inits <- vector(mode = "list", length = 8)
names(stan.inits) <- sort(c(paste0(paste0("model", c(1:4)), "_cp"),
                       paste0(paste0("model", c(1:4)), "_ncp")))
stan.draws <- vector(mode = "list", length = 8)
names(stan.draws) <- sort(c(paste0(paste0("model", c(1:4)), "_cp"),
                            paste0(paste0("model", c(1:4)), "_ncp")))
data.plot.list <- vector(mode = "list", length = 8)
# Now loop through each scenario and recreate the data as necessary by setting
# the true values of the appropriate coefficients to 0
# Only fit the model corresponding to the true data generating model
for (i in 1:4) {
  set.seed(1)
  if (i == 1 | i == 3) {
    gamma0 <- 0
  } else {
    gamma0 <- 0.8
  }
  beta0 <- rnorm(J, mean = mu0 + gamma0 * w, sd = tau0)
  beta0_rep <- rep(beta0, each = n.obs.per.grp)
  
  if (i == 1 | i == 2) {
    beta1 <- 0
  } else {
    beta1 <- 0.5
  }
  
  true.pars <- data.frame(mu0, gamma0, tau0, beta1, sigma_y)
  
  y <- rnorm(n, mean = beta0_rep + beta1 * x, sd = sigma_y)
  grp.df <- data.frame(id = c(1:J), w, beta0)
  unit.df <- data.frame(x, y, w = rep(w, each = n.obs.per.grp),
                        id = rep(c(1:J), each = n.obs.per.grp))
  
  # Run lmer
  lmer.res <- lmer(y ~ x + w + (1|id), data = unit.df)
  
  # Pull out BLUPS for random effects
  lmer.ranef <- ranef(lmer.res)$id
  names(lmer.ranef) <- "eta_j"
  lmer.ranef$id <- c(1:J)
  unit.df <- left_join(unit.df, lmer.ranef, by = "id")
  
  # Pull out coef ests for mu0, gamma0, beta
  lmer.coefs <- fixef(lmer.res)
  names(lmer.coefs) <- c("mu0", "beta1", "gamma0")
  lmer.coefs <- data.frame(t(lmer.coefs))
  cat("\n Coefs for", i, "\n")
  print(lmer.coefs)
  cat("\n Truth:\n")
  print(c(mu0, beta1, gamma0))
  lmer.coefs.gp <- lmer.coefs[rep(row.names(lmer.coefs), each = J), ]
  grp.df <- cbind(grp.df, lmer.coefs.gp)
  lmer.coefs.units <- lmer.coefs[rep(row.names(lmer.coefs), each = n.obs.per.grp*J), ]
  unit.df <- cbind(unit.df, lmer.coefs.units)
  
  # Plot data
  # data.plot.list[[2*i - 1]] <- ggplot(unit.df, aes(x = x, y = y)) +
  #   geom_point(aes(colour = factor(id))) +
  #   geom_line(aes(y = mu0 + gamma0 * w + beta1 * x + eta_j, colour = factor(id))) +
  #   scale_colour_discrete(guide = FALSE) +
  #   ggtitle(paste0("Model ", i)) +
  #   theme_bw()
  # data.plot.list[[2*i]] <- ggplot(grp.df, aes(x = w, y = beta0)) +
  #   geom_point() +
  #   geom_line(aes(y = mu0 + gamma0 * w)) +
  #   ggtitle(paste0("Model ", i)) +
  #   theme_bw()
  
  # Make stan data
  if (i == 1) {
    stan_data <- list(n = nrow(unit.df), J = J, grp_id = unit.df$id,
                      y = unit.df$y)
  } else if (i == 2) {
    stan_data <- list(n = nrow(unit.df), J = J, grp_id = unit.df$id,
                      y = unit.df$y, w = grp.df$w)
  } else if (i == 3) {
    stan_data <- list(n = nrow(unit.df), J = J, grp_id = unit.df$id,
                      y = unit.df$y, x = unit.df$x)
  } else {
    stan_data <- list(n = nrow(unit.df), J = J, grp_id = unit.df$id,
                      y = unit.df$y, x = unit.df$x, w = grp.df$w)
  }
  for (sfx in c("cp", "ncp")) {
    currmod <- paste0("model", i, "_", sfx)
    stan.res <- sampling(get(currmod),
                         data = stan_data, iter = 2000, chains = 4,
                         control = list(adapt_delta = 0.99))
    parnames <- stan.res@model_pars
    parnames <- parnames[!grepl("*eta0", parnames)]
    parnames <- parnames[!grepl("lp__", parnames)]
    parnames <- parnames[!grepl("ymean", parnames)]
    stan.draws[names(stan.draws) == currmod] <- extract(stan.res, pars = parnames)
    print(stan.res, pars = parnames)
    stan.summary <- summary(stan.res)$summary
    tau0.neff <- stan.summary["tau0", "n_eff"]
    stan.time <- mean(rowSums(data.frame(get_elapsed_time(stan.res))))
    ess.per.sec$tau0.neff[ess.per.sec$model == i &
                            ess.per.sec$param.type == sfx] <- tau0.neff / stan.time
    stan.inits[names(stan.inits) == currmod] <- get_inits(stan.res)
  }
}
#do.call(grid.arrange, c(data.plot.list, nrow = 4))




