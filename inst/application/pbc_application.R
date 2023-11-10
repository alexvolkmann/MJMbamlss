# Set up R session --------------------------------------------------------


library(survival)
library(JMbayes2)
library(bamlss)
library(MFPCA)
library(tidyverse)
library(parallel)
library(Rcpp)
library(Matrix)
library(sparseFLMM)
library(MJMbamlss)



# Prepare Data ------------------------------------------------------------


# Cox model for the composite event death or transplantation
pbc2$event <- as.numeric(pbc2$status != 'alive')

# Longitudinal format for pbc data
# Logarithmic scale for all three markers
p_long <- pbc2 %>%
  pivot_longer(c(serBilir, serChol, SGOT, albumin), names_to = "marker",
               values_to = "y") %>%
  mutate(survtime = years + sqrt(.Machine$double.eps),
         obstime = year + sqrt(.Machine$double.eps), marker = factor(marker),
         logy = log(y)) %>%
  select(id, survtime, event, sex, drug, age, marker, obstime, y, logy) %>%
  arrange(marker, id, obstime) %>%
  na.omit() %>%
  as.data.frame()

# Cap the survival times at the maximum observation times
max_obstime <- max(p_long$obstime)
p_long <- p_long %>%
  mutate(survtime = ifelse(survtime > max_obstime, max_obstime, survtime))

# Remove patients that do not have at least one serChol
which_mis <- p_long %>%
  filter(marker == "serChol") %>%
  group_by(id) %>%
  summarise(dplyr::n()) %>%
  select(id) %>%
  unlist()
p_long <- p_long %>%
  filter(id %in% which_mis) %>%
  droplevels()

# Data for JMbayes2
# Use logy for logarithmic transform
p_long_jmb <- p_long %>%
  pivot_wider(id_cols = c(id, obstime, survtime, event, sex, age, drug),
              values_from = logy, names_from = marker) 

p_long_id <- p_long_jmb %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup()



# Fit the JMbamlss Model --------------------------------------------------


# Choose the Number of Basis Functions ------------------------------------

# Estimate the model using estimated FPCs
few_obs <- apply(table(p_long$id, p_long$marker), 1,
                 function (x) any(x < 2))
long_obs <- p_long %>%
  group_by(id, marker) %>%
  summarize(maxobs = max(obstime), .groups = "drop_last") %>%
  ungroup(marker) %>%
  summarize(minmaxobs = min(maxobs), .groups = "drop_last") %>%
  ungroup() %>%
  filter(minmaxobs > 0.1 * max(p_long$survtime)) %>%
  select(id) %>%
  unlist() %>%
  paste()
take <- intersect(long_obs, paste(which(!few_obs)))

# Final choice
mfpca_w <- preproc_MFPCA(
  p_long %>% filter(id %in% take) %>% droplevels(),
  uni_mean = paste0("logy ~ 1 + sex + drug + s(obstime, k = 10, bs = 'ps') + ",
                    "s(age, k = 10, bs = 'ps')"),
  pve_uni = 0.99, nbasis = 7, weights = TRUE, save_uniFPCA = TRUE)

# Sensitivity
mfpca_w_4 <- preproc_MFPCA(
  p_long %>% filter(id %in% take) %>% droplevels(),
  uni_mean = paste0("logy ~ 1 + sex + drug + s(obstime, k = 10, bs = 'ps') + ",
                    "s(age, k = 10, bs = 'ps')"),
  pve_uni = 0.99, nbasis = 4, weights = TRUE, save_uniFPCA = TRUE)

mfpca_w_10 <- preproc_MFPCA(
  p_long %>% filter(id %in% take) %>% droplevels(),
  uni_mean = paste0("logy ~ 1 + sex + drug + s(obstime, k = 10, bs = 'ps') + ",
                    "s(age, k = 10, bs = 'ps')"),
  pve_uni = 0.99, nbasis = 10, weights = TRUE, save_uniFPCA = TRUE)

mfpca_u <- preproc_MFPCA(
  p_long %>% filter(id %in% take) %>% droplevels(),
  uni_mean = paste0("logy ~ 1 + sex + drug + s(obstime, k = 10, bs = 'ps') + ",
                    "s(age, k = 10, bs = 'ps')"),
  pve_uni = 0.99, nbasis = 7, weights = FALSE, save_uniFPCA = TRUE)


# Extract the information of first two Eigenfunctions and convert to data frame
refObj <- extractObs(mfpca_w$functions, 1:7)
fpcs <- lapply(list(mfpca_w_4, mfpca_w, mfpca_w_10, mfpca_u), function(x) {
  dat <- extractObs(x$functions, 1:7)
  dat <- flipFuns(refObj, dat)
  multifammPaper::funData2DataFrame(dat)
})
fpc_dat <- do.call(
  rbind, Map(cbind,
             nbasis = factor(1:4, labels = c("4B", "7B", "10B", "7Bu")),
             fpcs)) %>%
  mutate(dim = factor(dim,
                      labels = c("albumin", "serBilir", "serChol", "SGOT")),
         obs = factor(obs, labels = paste("MFPC", 1:7)))

ggplot(fpc_dat, aes(x = t, y = y, color = nbasis)) +
  geom_line(aes(linetype = nbasis)) +
  facet_grid(obs ~ dim, scales = "free_y") +
  theme_bw() +
  labs(y = NULL, x = "Time", color = NULL, linetype = NULL)
# 6x12

# Numerical results
sapply(attr(mfpca_w, "uniFPCA"), "[[", "npc")
sapply(attr(mfpca_w, "uniFPCA"), function (x) sum(x$evalues))
mfpca_w$values / sum(mfpca_w$values)
mfpca_u$values / sum(mfpca_u$values)


# Model using 7 Weighted Basis Functions ----------------------------------

vals <- which(mfpca_w$values > 0)
nfpc <- min(which(
  cumsum(mfpca_w$values[vals])/sum(mfpca_w$values[vals]) > 0.99))
mfpca_w_list <- lapply(vals, function (i, mfpca = mfpca_w) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})

# Prepare objects for model fit
p_long[, grepl("fpc", colnames(p_long))] <- NULL
p_long <- attach_wfpc(mfpca_w, p_long, n = length(vals))

# Model formula
f_est <- list(
  Surv2(survtime, event, obs = logy) ~ -1 + 
    s(survtime, k = 10, bs = "ps", xt = list("scale" = FALSE)),
  gamma ~ 1 + s(age, k = 10, bs = "ps") + sex + drug,
  as.formula(paste0(
    "mu ~ -1 + marker + sex:marker + drug:marker + ",
    "s(obstime, by = marker, xt = list('scale' = FALSE), k = 10, bs = 'ps') + ",
    "s(age, by = marker, xt = list('scale' = FALSE), k = 10, bs = 'ps') + ",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_w_list[[", x, "]], 'scale' = FALSE))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

set.seed(1209)
t_w <- system.time(
  bamlss_wei <- bamlss(f_est, family = mjm_bamlss, data = p_long,
                       timevar = "obstime", maxit = 1500,
                       n.iter = 12000L, burnin = 2000L, thin = 5)
)
# > t_w
# user   system  elapsed 
# 9250.846   14.366 9262.599 




# Model using unweighted MFPCA --------------------------------------------

vals <- which(mfpca_u$values > 0)
nfpc <- min(which(
  cumsum(mfpca_u$values[vals])/sum(mfpca_u$values[vals]) > 0.99))
mfpca_u_list <- lapply(vals, function (i, mfpca = mfpca_u) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})

# Prepare objects for model fit
p_long[, grepl("fpc", colnames(p_long))] <- NULL
p_long <- attach_wfpc(mfpca_u, p_long, n = length(vals))

# Model formula
f_est <- list(
  Surv2(survtime, event, obs = logy) ~ -1 + 
    s(survtime, k = 10, bs = "ps", xt = list("scale" = FALSE)),
  gamma ~ 1 + s(age, k = 10, bs = "ps") + sex + drug,
  as.formula(paste0(
    "mu ~ -1 + marker + sex:marker + drug:marker + ",
    "s(obstime, by = marker, xt = list('scale' = FALSE), k = 10, bs = 'ps') + ",
    "s(age, by = marker, xt = list('scale' = FALSE), k = 10, bs = 'ps') + ",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_u_list[[", x, "]], 'scale' = FALSE))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

set.seed(1209)
t_w <- system.time(
  bamlss_unw <- bamlss(f_est, family = mjm_bamlss, data = p_long,
                       timevar = "obstime", maxit = 1500,
                       n.iter = 12000L, burnin = 2000L, thin = 5)
)




# JMbayes2 Model ----------------------------------------------------------

# Get the quantile-based knots for comparability
kn <- mgcv::smoothCon(mgcv::s(survtime, k = 10, bs = "ps"), 
                      data = p_long)[[1]]$knots

set.seed(1205)
# Cox Model
t_jmb <- system.time({
  
  CoxFit <- coxph(Surv(survtime, event) ~ ns(age, df = 3) + sex + drug,
                  data = p_long_id)
  
  # Univariate longitudinal models
  fm1ns <- lme(albumin ~ ns(obstime, df = 3) + sex + drug + ns(age, df = 3),
               data = p_long_jmb, na.action = na.omit,
               random = ~ ns(obstime, df = 2) | id)
  fm2ns <- lme(serBilir ~ ns(obstime, df = 3) + sex + drug + ns(age, df = 3), 
               data = p_long_jmb, na.action = na.omit, 
               random = ~ ns(obstime, df = 2) | id)
  fm3ns <- lme(serChol ~ ns(obstime, df = 3) + sex + drug + ns(age, df = 3),
               data = p_long_jmb, na.action = na.omit, 
               random = ~ ns(obstime, df = 2) | id)
  fm4ns <- lme(SGOT ~ ns(obstime, df = 3) + sex + drug + ns(age, df = 3),
               data = p_long_jmb, na.action = na.omit, 
               random = ~ ns(obstime, df = 2) | id)
  
  
  # Multivariate Joint Model
  # the joint model that links all sub-models
  jointFit <- jm(CoxFit, list(fm1ns, fm2ns, fm3ns, fm4ns), time_var = "obstime",
                 n_iter = 12000L, n_burnin = 2000L, n_thin = 5L,
                 cores = 1, n_chains = 1, 
                 GK_k = 7, Bsplines_degree = 3, diff = 2, knots = list(kn),
                 save_random_effects = TRUE)
  
})



# Compare Fitted Observations ---------------------------------------------


# Create data set for prediction
id_grid <- sapply(p_long_id$survtime, function (x) {
  c(seq(0 + sqrt(.Machine$double.eps), x, by = 0.2), x)
})
id_grid_length <- sapply(id_grid, length)
newdat <- p_long_id[rep(seq_len(nrow(p_long_id)), id_grid_length), ] %>%
  mutate(obstime = unlist(id_grid)) %>%
  pivot_longer(cols = c(albumin, serBilir, serChol, SGOT), 
               names_to = "marker") %>%
  mutate(marker = as.factor(marker),
         value = NULL) %>%
  arrange(marker, id, obstime) %>%
  droplevels() %>%
  as.data.frame()

# Prediction of bamlss model fit
preds <- cbind(
  newdat,
  predict(bamlss_wei, newdata = newdat, model = "mu", FUN = c95))

# Prediction of JMbayes2 model fit
# ns() objects have to be predicted from original data set otherwise the values
# do not align
X <- lapply(split(newdat, newdat$marker), function(x) {
  dat <- cbind(predict(ns(p_long_jmb$obstime, df = 3), x$obstime), 
               model.matrix(~ sex + drug, x),
               predict(ns(p_long_jmb$age, df = 3), x$age))
  colnames(dat)[1:3] <- paste0("ns(obstime, df = 3)", 1:3)
  colnames(dat)[7:9] <- paste0("ns(age, df = 3)", 1:3)
  dat[, c(4, 1:3, 5:9)]
})
Z <- lapply(split(newdat, newdat$marker), function(x) {
  dat <- cbind(model.matrix(~ 1, x),
               predict(ns(p_long_jmb$obstime, df = 2), x$obstime))
  colnames(dat)[2:3] <- paste0("ns(obstime, df = 2)", 1:2)
  dat
})
B <- jointFit$mcmc$b[[1]]
n_re <- ncol(Z[[1]])
mcmc_mu <- do.call(rbind, lapply(seq_len(4), function (dim) {
  idL <- unclass(droplevels(split(newdat$id, newdat$marker)[[dim]]))
  tcrossprod(X[[dim]], jointFit$mcmc[[paste0("betas", dim)]][[1]]) +
    t(sapply(seq_len(nrow(Z[[dim]])), function (i) {
      Z[[dim]][i, ] %*% B[idL[i], (dim - 1)*n_re + seq_len(n_re), ]
    }))
}))
jmb_preds <- data.frame("jmb_2.5" = apply(mcmc_mu, 1, quantile, probs = 0.025),
                        "jmb_Mean" = rowMeans(mcmc_mu),
                        "jmb_97.5" = apply(mcmc_mu, 1, quantile, probs = 0.975))
preds <- cbind(preds, jmb_preds)

# Reform the data set for nicer plotting
preds_plot <- pivot_longer(preds, cols = c("Mean", "jmb_Mean"), 
                           names_to = "Model", values_to = "fit") %>%
  mutate(Model = factor(Model, levels = c("Mean", "jmb_Mean"), 
                        labels = c("bamlss", "JMB")))


# ids <- sample(levels(p_long$id), size = 5)
ids <- c(40, 54, 114, 204)

ggplot(preds_plot %>% 
         filter(id %in% ids) %>%
         droplevels() %>%
         mutate(id = factor(id, labels = paste("Subject", levels(id)))),
       aes(x = obstime, y = fit)) + 
  geom_line(aes(linetype = Model, color = Model)) +
  theme_bw() +
  facet_grid(marker ~ id, scales = "free_y") +
  geom_point(data = p_long %>% 
               filter(id %in% ids) %>%
               droplevels() %>%
               mutate(id = factor(id, labels = paste("Subject", levels(id)))),
             aes(y = logy), alpha = 0.5) +
  labs(y = NULL, x = "Time")
# save as 4x8




# Comparison Weighted - Unweighted Fit ------------------------------------

# Use sum of residual squares
p_long$wei <- predict(bamlss_wei, model = "mu")
p_long$unw <- predict(bamlss_unw, model = "mu")

tapply((p_long$wei - p_long$logy)^2, p_long$marker, sum)
# albumin  serBilir   serChol      SGOT 
# 16.60007 129.82868  24.19605 107.39504 

tapply((p_long$unw - p_long$logy)^2, p_long$marker, sum)
# albumin  serBilir   serChol      SGOT 
# 20.35439 120.98950  35.37122 115.17453 

# DIC
bamlss_wei$model.stats$sampler$DIC
# [1] -896.1151
bamlss_unw$model.stats$sampler$DIC
# [1] -579.9308



# Table for Estimated Parameters ------------------------------------------

# Association parameters
aj_mean <- jointFit$statistics$Mean$alphas
aj_low <- jointFit$statistics$CI_low$alphas
aj_up <- jointFit$statistics$CI_upp$alphas

alpha_bw <- summary(samples(bamlss_wei)[, grep("alpha\\.",
                                               colnames(samples(bamlss_wei)))])
abw_mean <- alpha_bw$statistics[1:4, 1]
abw_low <- alpha_bw$quantiles[1:4, 1]
abw_up <- alpha_bw$quantiles[1:4, 5]

cat(mapply(function(marker, b_mean, b_low, b_up, c_mean,
                    c_low, c_up){
  paste("\\textit{", marker,"}",  
        "&$", round(b_mean, digits = 2), "$ & $[", round(b_low, digits = 2),
        "$; $", round(b_up, digits = 2), "]$ & ",
        "$", round(c_mean, digits = 2), "$ & $[", round(c_low, digits = 2),
        "$;  $", round(c_up, digits = 2), "]$ \\\\ ")
}, c("albumin", "serBilir", "serChol", "SGOT"), abw_mean, abw_low, abw_up,
aj_mean, aj_low, aj_up))

# Other parameters

# JMB mean
j_gamma <- jointFit$statistics$Mean$gammas[c(4, 5)]
j_mu1 <- jointFit$statistics$Mean$betas1[c(1, 5, 6)]
j_mu2 <- jointFit$statistics$Mean$betas2[c(1, 5, 6)]
j_mu3 <- jointFit$statistics$Mean$betas3[c(1, 5, 6)]
j_mu4 <- jointFit$statistics$Mean$betas4[c(1, 5, 6)]
j_sigma <- log(jointFit$statistics$Mean$sigmas)
j_mean <- c(j_gamma, j_mu1, j_mu2, j_mu3, j_mu4, j_sigma)

# JMB CI
j_gamma_up <- jointFit$statistics$CI_up$gammas[c(4, 5)]
j_mu1_up <- jointFit$statistics$CI_up$betas1[c(1, 5, 6)]
j_mu2_up <- jointFit$statistics$CI_up$betas2[c(1, 5, 6)]
j_mu3_up <- jointFit$statistics$CI_up$betas3[c(1, 5, 6)]
j_mu4_up <- jointFit$statistics$CI_up$betas4[c(1, 5, 6)]
j_sigma_up <- log(jointFit$statistics$CI_up$sigmas)
j_up <- c(j_gamma_up, j_mu1_up, j_mu2_up, j_mu3_up, j_mu4_up, j_sigma_up)

j_gamma_low <- jointFit$statistics$CI_low$gammas[c(4, 5)]
j_mu1_low <- jointFit$statistics$CI_low$betas1[c(1, 5, 6)]
j_mu2_low <- jointFit$statistics$CI_low$betas2[c(1, 5, 6)]
j_mu3_low <- jointFit$statistics$CI_low$betas3[c(1, 5, 6)]
j_mu4_low <- jointFit$statistics$CI_low$betas4[c(1, 5, 6)]
j_sigma_low <- log(jointFit$statistics$CI_low$sigmas)
j_low <- c(j_gamma_low, j_mu1_low, j_mu2_low, j_mu3_low, j_mu4_low, j_sigma_low)

# Bamlss relevants
rels <- paste0("(\\.s\\.|edf$|accepted$|alpha|\\(Intercept\\)|logLik|",
               "logPost|DIC|pd)")
b <- summary(samples(bamlss_wei)[, -grep(rels, colnames(samples(bamlss_wei)))])
ords <- c(1:3, 7, 11, 4, 8, 12, 5, 9, 13, 6, 10, 14, 15:18)
b_mean <- b$statistics[ords, 1]
b_low <- b$quantiles[ords, 1]
b_up <- b$quantiles[ords, 5]


cat(mapply(function(predictor, b_mean, b_low, b_up, c_mean,
                    c_low, c_up){
  paste0("$", predictor, "$",
         " & $", round(b_mean, digits = 2), "$ & $[", round(b_low, digits = 2),
         "$; $", round(b_up, digits = 2), "]$ & ",
         "$", round(c_mean, digits = 2), "$ & $[", round(c_low, digits = 2),
         "$;  $", round(c_up, digits = 2), "]$ \\\\ ")
}, predictor = c("\\beta_{\\gamma1}", "\\beta_{\\gamma2}", "\\beta_{albumin0}",
                 "\\beta_{albumin1}", "\\beta_{albumin2}", 
                 "\\beta_{serBilir0}", "\\beta_{serBilir1}",
                 "\\beta_{serBilir2}", "\\beta_{serChol0}", "\\beta_{serChol1}",
                 "\\beta_{serChol2}", "\\beta_{SGOT0}", "\\beta_{SGOT1}",
                 "\\beta_{SGOT2}", "\\beta_{\\sigma albumin}", 
                 "\\beta_{\\sigma serBilir}", "\\beta_{\\sigma serChol}",
                 "\\beta_{\\sigma SGOT}"),
b_mean, b_low, b_up, j_mean, j_low, j_up))

