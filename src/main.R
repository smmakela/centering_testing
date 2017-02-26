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
    tmp <- stanc(file = code_file, model_name = nam)
    assign(nam, stan_model(stanc_ret = tmp, model_name = nam))
  }
}

# Initialize some parameters
J.list <- c(5, 15, 30) # number of groups
tau0.list <- 10^c(-3:1) # group-level variance
n.obs.per.grp <- 30
sigma_y <- 0.05
mu0 <- 0.2

# Initialize lists/data frames where we'll be storing simulation results
# num.combs is number of models * cp/ncp * length(J.list) * length(tau0.list)
num.combs <- 4 * 2 * length(J.list) * length(tau0.list)
ess.per.sec <- expand.grid(model = c(1:4), param.type = c("cp", "ncp"),
                           tau0.value = tau0.list, J = J.list)
ess.per.sec$tau0.neff.per.sec <- rep(NA, length = nrow(ess.per.sec))
true.pars <- vector(mode = "list", length = num.combs)
# stan.inits <- vector(mode = "list", length = num.combs)
# names(stan.inits) <- sort(c(paste0(paste0("model", c(1:4)), "_cp"),
#                             paste0(paste0("model", c(1:4)), "_ncp")))
# stan.draws <- vector(mode = "list", length = num.combs)
# names(stan.draws) <- sort(c(paste0(paste0("model", c(1:4)), "_cp"),
#                             paste0(paste0("model", c(1:4)), "_ncp")))
# data.plot.list <- vector(mode = "list", length = 8)

counter <- 1

for (j in 1:length(J.list)) {
  for (t in 1:length(tau0.list)) {
    # Initialize parameters
    J <- J.list[j]
    tau0 <- tau0.list[t]
    n <- n.obs.per.grp*J
    w <- runif(J, min = -1, max = 1)
    x <- runif(n, min = -1, max = 1)
    
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
      y <- rnorm(n, mean = beta0_rep + beta1 * x, sd = sigma_y)
      
      true.pars[[counter]] <- list(J = J, mu0 = mu0, gamma0 = gamma0,
                                   tau0 = tau0, beta1 = beta1, sigma_y = sigma_y,
                                   beta0 = beta0, x = x, y = y, w = w)
      
      # Store group and unit data separately
      grp.df <- data.frame(id = c(1:J), w, beta0)
      unit.df <- data.frame(x, y, w = rep(w, each = n.obs.per.grp),
                            id = rep(c(1:J), each = n.obs.per.grp))
      
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
      
      # Loop through centered and noncentered parameterizations
      for (sfx in c("cp", "ncp")) {
        curr.mod <- paste0("model", i, "_", sfx)
        cat("\n===================================================================\n")
        cat("\n J=", J, "tau0 =", tau0, "CURRENT MODEL:", curr.mod, "\n")
        num.div.trans <- 1
        delta.counter <- 1
        adapt.delta.list <- c(0.8, 0.9, 0.99, 0.999, 0.9999)
        while (num.div.trans > 0 & delta.counter <= length(adapt.delta.list)) {
          adapt.delta <- adapt.delta.list[delta.counter]
          stan.res <- sampling(get(curr.mod),
                               data = stan_data, iter = 2000, chains = 4,
                               control = list(adapt_delta = adapt.delta))
          # Count number of divergent transitions after warmup
          samp.params <- get_sampler_params(stan.res)
          num.div.trans <- 0
          for (s in 1:length(samp.params)) {
            num.div.trans <- num.div.trans + sum(samp.params[[s]][1001:2000, 5])
          }
          delta.counter <- delta.counter + 1
        }
        #parnames <- stan.res@model_pars
        #parnames <- parnames[!grepl("*eta0", parnames)]
        #parnames <- parnames[!grepl("lp__", parnames)]
        #parnames <- parnames[!grepl("ymean", parnames)]
        #stan.draws[names(stan.draws) == currmod] <- extract(stan.res, pars = parnames)
        #print(stan.res, pars = parnames)
        stan.summary <- summary(stan.res)$summary
        tau0.neff <- stan.summary["tau0", "n_eff"]
        stan.time <- mean(rowSums(data.frame(get_elapsed_time(stan.res))))
        ess.ind <- ess.per.sec$model == i & ess.per.sec$param.type == sfx &
          ess.per.sec$J == J & ess.per.sec$tau0.value == tau0
        ess.per.sec$tau0.neff[ess.ind] <- tau0.neff
        ess.per.sec$time[ess.ind] <- stan.time
        ess.per.sec$tau0.neff.per.sec[ess.ind] <- tau0.neff / stan.time
        ess.per.sec$J[ess.ind] <- J
        ess.per.sec$tau0.value[ess.ind] <- tau0
        ess.per.sec$num.div.trans[ess.ind] <- num.div.trans
        ess.per.sec$adapt.delta[ess.ind] <- adapt.delta
        ess.per.sec$delta.counter.m1[ess.ind] <- delta.counter - 1
        #stan.inits[names(stan.inits) == currmod] <- get_inits(stan.res)
      } # end sfx loop (cp/ncp)
      counter <- counter + 1
    } # end i loop (models 1 -4)
  } # end tau0 loop
} # end J loop
# # Now loop through each scenario and recreate the data as necessary by setting
# # the true values of the appropriate coefficients to 0
# # Only fit the model corresponding to the true data generating model
# for (i in 1:4) {
#   set.seed(1)
#   if (i == 1 | i == 3) {
#     gamma0 <- 0
#   } else {
#     gamma0 <- 0.8
#   }
#   beta0 <- rnorm(J, mean = mu0 + gamma0 * w, sd = tau0)
#   beta0_rep <- rep(beta0, each = n.obs.per.grp)
#   
#   if (i == 1 | i == 2) {
#     beta1 <- 0
#   } else {
#     beta1 <- 0.5
#   }
#   
#   true.pars <- data.frame(mu0, gamma0, tau0, beta1, sigma_y)
#   
#   y <- rnorm(n, mean = beta0_rep + beta1 * x, sd = sigma_y)
#   grp.df <- data.frame(id = c(1:J), w, beta0)
#   unit.df <- data.frame(x, y, w = rep(w, each = n.obs.per.grp),
#                         id = rep(c(1:J), each = n.obs.per.grp))
#   
#   # Run lmer
#   lmer.res <- lmer(y ~ x + w + (1|id), data = unit.df)
#   
#   # Pull out BLUPS for random effects
#   lmer.ranef <- ranef(lmer.res)$id
#   names(lmer.ranef) <- "eta_j"
#   lmer.ranef$id <- c(1:J)
#   unit.df <- left_join(unit.df, lmer.ranef, by = "id")
#   
#   # Pull out coef ests for mu0, gamma0, beta
#   lmer.coefs <- fixef(lmer.res)
#   names(lmer.coefs) <- c("mu0", "beta1", "gamma0")
#   lmer.coefs <- data.frame(t(lmer.coefs))
#   cat("\n Coefs for", i, "\n")
#   print(lmer.coefs)
#   cat("\n Truth:\n")
#   print(c(mu0, beta1, gamma0))
#   lmer.coefs.gp <- lmer.coefs[rep(row.names(lmer.coefs), each = J), ]
#   grp.df <- cbind(grp.df, lmer.coefs.gp)
#   lmer.coefs.units <- lmer.coefs[rep(row.names(lmer.coefs), each = n.obs.per.grp*J), ]
#   unit.df <- cbind(unit.df, lmer.coefs.units)
#   
#   # Plot data
#   # data.plot.list[[2*i - 1]] <- ggplot(unit.df, aes(x = x, y = y)) +
#   #   geom_point(aes(colour = factor(id))) +
#   #   geom_line(aes(y = mu0 + gamma0 * w + beta1 * x + eta_j, colour = factor(id))) +
#   #   scale_colour_discrete(guide = FALSE) +
#   #   ggtitle(paste0("Model ", i)) +
#   #   theme_bw()
#   # data.plot.list[[2*i]] <- ggplot(grp.df, aes(x = w, y = beta0)) +
#   #   geom_point() +
#   #   geom_line(aes(y = mu0 + gamma0 * w)) +
#   #   ggtitle(paste0("Model ", i)) +
#   #   theme_bw()
#   
#   # Make stan data
#   if (i == 1) {
#     stan_data <- list(n = nrow(unit.df), J = J, grp_id = unit.df$id,
#                       y = unit.df$y)
#   } else if (i == 2) {
#     stan_data <- list(n = nrow(unit.df), J = J, grp_id = unit.df$id,
#                       y = unit.df$y, w = grp.df$w)
#   } else if (i == 3) {
#     stan_data <- list(n = nrow(unit.df), J = J, grp_id = unit.df$id,
#                       y = unit.df$y, x = unit.df$x)
#   } else {
#     stan_data <- list(n = nrow(unit.df), J = J, grp_id = unit.df$id,
#                       y = unit.df$y, x = unit.df$x, w = grp.df$w)
#   }
#   for (sfx in c("cp", "ncp")) {
#     currmod <- paste0("model", i, "_", sfx)
#     stan.res <- sampling(get(currmod),
#                          data = stan_data, iter = 2000, chains = 4,
#                          control = list(adapt_delta = 0.99))
#     parnames <- stan.res@model_pars
#     parnames <- parnames[!grepl("*eta0", parnames)]
#     parnames <- parnames[!grepl("lp__", parnames)]
#     parnames <- parnames[!grepl("ymean", parnames)]
#     stan.draws[names(stan.draws) == currmod] <- extract(stan.res, pars = parnames)
#     print(stan.res, pars = parnames)
#     stan.summary <- summary(stan.res)$summary
#     tau0.neff <- stan.summary["tau0", "n_eff"]
#     stan.time <- mean(rowSums(data.frame(get_elapsed_time(stan.res))))
#     ess.per.sec$tau0.neff[ess.per.sec$model == i &
#                             ess.per.sec$param.type == sfx] <- tau0.neff / stan.time
#     stan.inits[names(stan.inits) == currmod] <- get_inits(stan.res)
#   }
# }
#do.call(grid.arrange, c(data.plot.list, nrow = 4))

ess.per.sec$model.name <- paste0("Model ", ess.per.sec$model)
ess.per.sec$J.name <- paste0("J = ", ess.per.sec$J)
ess.per.sec$J.name <- factor(ess.per.sec$J.name,
                             levels = c("J = 5", "J = 15", "J = 30"))
ggplot(ess.per.sec, aes(x = factor(tau0.value), y = log10(tau0.neff.per.sec))) +
  geom_point(aes(colour = param.type)) +
  geom_line(aes(colour = param.type, group = param.type)) +
  scale_colour_discrete("Parameterization") +
  # geom_point(aes(colour = param.type, shape = factor(num.div.trans == 0)),
  #            size = 2) +
  # scale_shape_manual("Divergent transitions",
  #                    values = c(1, 16),
  #                    breaks = c(FALSE, TRUE),
  #                    labels = c("Yes", "No")) +
  xlab("tau0") +
  ylab("ESS per sec (log10)") +
  facet_grid(J.name ~ model.name) +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("/Users/susanna/projects/centering_testing/output/figures/tau0_ess_per_sec.pdf",
       width = 10, height = 8)
ggsave("/Users/susanna/projects/centering_testing/output/figures/tau0_ess_per_sec.png",
       width = 9, height = 6)

