


# Specify location
# setwd()
# results_wd <- 


# Always
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


# Evaluate bamlss TRUE FPCs -----------------------------------------------

mnames_btru <- list.files(path = file.path(results_wd, "scen_II", 
                                           "bamlss_tru"))
preds_btru <- MJMbamlss:::sim_bamlss_predict(mnames_btru, results_wd, 
                                            "/scen_II/bamlss_tru/",
                                            "/scen_II/data/", rds = TRUE,
                                            old = TRUE)
saveRDS(preds_btru,
        file = file.path(results_wd, "scen_II", "preds_btru.rds"))

it_list <- MJMbamlss:::sim_results(lapply(preds_btru, "[[", "predictions"),
                                  lapply(preds_btru, "[[", "simulations"),
                                  name = "TRU")
eval_btru <- do.call(rbind, Map(cbind, it = sub("\\.rds", "", names(it_list)),
                               it_list))
saveRDS(eval_btru,
        file = file.path(results_wd, "scen_II", "eval_btru.rds"))

rm(preds_btru, eval_btru, it_list)



# Evaluate bamlss EST 1 FPCs ----------------------------------------------

mnames_best1 <- list.files(path = file.path(results_wd, "scen_II", 
                                            "bamlss_est1"))
preds_best1 <- MJMbamlss:::sim_bamlss_predict(mnames_best1, results_wd, 
                                             "/scen_II/bamlss_est1/",
                                             "/scen_II/data/", 
                                             rds = TRUE, old = TRUE)
saveRDS(preds_best1,
        file = file.path(results_wd, "scen_II", "preds_best1.rds"))


it_list <- MJMbamlss:::sim_results(lapply(preds_best1, "[[", "predictions"),
                                  lapply(preds_best1, "[[", "simulations"),
                                  name = "EST1")
eval_best1 <- do.call(rbind, Map(cbind, it = sub("\\.rds", "", names(it_list)),
                                 it_list))
saveRDS(eval_best1,
        file = file.path(results_wd, "scen_II", "eval_best1.rds"))

rm(preds_best1, eval_best1, it_list)


# Evaluate bamlss EST 99 FPCs ---------------------------------------------

mnames_best99 <- list.files(path = file.path(results_wd, "scen_II", 
                                             "bamlss_est99"))
preds_best99 <- MJMbamlss:::sim_bamlss_predict(mnames_best99, results_wd, 
                                              "/scen_II/bamlss_est99/",
                                              "/scen_II/data/",
                                              rds = TRUE)
saveRDS(preds_best99,
        file = file.path(results_wd, "scen_II", "preds_best99.rds"))

it_list <- MJMbamlss:::sim_results(lapply(preds_best99, "[[", "predictions"),
                                  lapply(preds_best99, "[[", "simulations"),
                                  name = "EST99")
eval_best99 <- do.call(rbind, Map(cbind, it = sub("\\.rds", "", names(it_list)),
                                  it_list))
saveRDS(eval_best99,
        file = file.path(results_wd, "scen_II", "eval_best99.rds"))
rm(preds_best99, eval_best99, it_list)





# Evaluate jmbayes  -------------------------------------------------------

mnames_jmb <- list.files(path = file.path(results_wd, "scen_II", "jmb"))
preds_jmb <- MJMbamlss:::sim_jmb_predict(mnames_jmb, results_wd, 
                                        "/scen_II/jmb/",
                                        "/scen_II/data/", rds = TRUE)
saveRDS(preds_jmb,
        file = file.path(results_wd, "scen_II", "preds_jmb.rds"))


it_list <- MJMbamlss:::sim_results(lapply(preds_jmb, "[[", "predictions"),
                                  lapply(preds_jmb, "[[", "simulations"),
                                  name = "JMB")
eval_jmb <- do.call(rbind, Map(cbind, it = sub("\\.rds", "", names(it_list)),
                               it_list))
saveRDS(eval_jmb,
        file = file.path(results_wd, "scen_II", "eval_jmb.rds"))

rm(preds_jmb, eval_jmb, it_list)




# Table for Paper ---------------------------------------------------------

e_btru <- readRDS(file.path(results_wd, "scen_II", "eval_btru.rds"))
e_best1 <- readRDS(file.path(results_wd, "scen_II", "eval_best1.rds"))
e_best99 <- readRDS(file.path(results_wd, "scen_II", "eval_best99.rds"))
e_jmb <- readRDS(file.path(results_wd, "scen_II", "eval_jmb.rds"))


r <- rbind(e_btru, e_best1, e_best99, e_jmb) %>%
  mutate(Model = factor(model, levels = c("TRU", "EST1", "EST99", "JMB"),
                        labels = c("TRUE", "EST", "TRUNC", "JMB")))

bias <- r %>%
  filter(type == "Bias", 
         !(predictor %in% c("mu_long", "lambga", "lambga_long"))) %>%
  group_by(Model, predictor, marker) %>%
  summarise(mean = mean(value), .groups = "drop") %>%
  mutate(mean = ifelse(predictor == "mu", mean*100, mean)) %>%
  pivot_wider(id_cols = c(predictor, marker), names_from = Model, 
              values_from = mean) %>%
  # rowid_to_column() %>%
  # mutate(rowid = ifelse(predictor == "sigma", 16, rowid),
  #        marker = ifelse(predictor == "sigma", "all", marker)) %>%
  # group_by(rowid, predictor, marker) %>%
  # summarise(across(1:4, mean), .groups = "drop") %>%
  slice(3, 1:2, 4:7)

mse <- r %>%
  filter(type == "rMSE", 
         !(predictor %in% c("mu_long", "lambga", "lambga_long"))) %>%
  group_by(Model, predictor, marker) %>%
  summarise(mean = mean(value), .groups = "drop") %>%
  mutate(mean = ifelse(predictor == "mu", mean*100, mean)) %>%
  pivot_wider(id_cols = c(predictor, marker), names_from = Model, 
              values_from = mean) %>%
  # rowid_to_column() %>%
  # mutate(rowid = ifelse(predictor == "sigma", 16, rowid),
  #        marker = ifelse(predictor == "sigma", "all", marker)) %>%
  # group_by(rowid, predictor, marker) %>%
  # summarise(across(1:4, mean), .groups = "drop") %>%
  slice(3, 1:2, 4:7)

coverage <- r %>%
  filter(type == "Coverage", 
         !(predictor %in% c("mu_long", "lambga", "lambga_long"))) %>%
  group_by(Model, predictor, marker) %>%
  summarise(mean = mean(value), .groups = "drop") %>%
  pivot_wider(id_cols = c(predictor, marker), names_from = Model, 
              values_from = mean) %>%
  # rowid_to_column() %>%
  # mutate(rowid = ifelse(predictor == "sigma", 16, rowid),
  #        marker = ifelse(predictor == "sigma", "all", marker)) %>%
  # group_by(rowid, predictor, marker) %>%
  # summarise(across(1:4, mean), .groups = "drop") %>%
  slice(3, 1:2, 4:7)

eval <- left_join(bias, mse, by = c("predictor", "marker"), 
                  suffix = c("bias", "mse")) %>%
  left_join(coverage, by = c("predictor", "marker")) %>%
  mutate(rowid = NULL,
         marker = NULL,
         predictor = c("$\\lambda + \\gamma$", paste0("$\\alpha_{", 1:2, "}$"),
                       paste0("$\\mu_{", 1:2, "}$"), 
                       paste0("$\\sigma_{", 1:2, "}$")))

xtable::print.xtable(xtable::xtable(eval, digits = 3), include.rownames = FALSE,
                     sanitize.text.function = identity)

# Time dependent evaluation -----------------------------------------------

ggplot(e_best1 %>% 
         filter(type == "MSE", predictor == "lambga_long"),
       aes(x = t, y = value))+
  geom_boxplot() +
  theme_bw()
ggplot(e_best1 %>% 
         filter(type == "Coverage", predictor == "mu_long"),
       aes(x = t, y = value))+
  geom_boxplot() +
  theme_bw() + 
  facet_wrap(~marker, nrow = 2) + 
  # ylim(c(-0.25, 0.1)) + ggtitle("Bias") # Bias
  # ylim(c(0, 0.05)) + ggtitle("MSE") # Change to type == MSE
  ggtitle("Coverage") # Change to type == Coverage

# JMbayes2
ggplot(e_jmb %>% 
         filter(type == "Coverage", predictor == "lambga_long"),
       aes(x = t, y = value))+
  geom_boxplot() +
  theme_bw() +
  ggtitle("Coverage")
ggplot(e_jmb %>% 
         filter(type == "Bias", predictor == "mu_long"),
       aes(x = t, y = value))+
  geom_boxplot() +
  theme_bw() + 
  facet_wrap(~marker, nrow = 2) + 
  ggtitle("Bias") # Bias
  # ylim(c(0, 1)) + ggtitle("MSE") # Change to type == MSE
  # ggtitle("Coverage") # Change to type == Coverage



# Plot longitudinal fits --------------------------------------------------

d_sim <- readRDS(file.path(results_wd, "scen_II", "data", "d100.rds"))
# set.seed(1444)
# ids <- c(sample(which(d_sim$data_short$survtime[1:150] > 0.2 & 
#                         d_sim$data_short$survtime[1:150] < 0.35), size = 1),
#          sample(which(d_sim$data_short$survtime[1:150] > 0.35 & 
#                         d_sim$data_short$survtime[1:150] < 0.7), size = 1),
#          sample(which(d_sim$data_short$survtime[1:150] > 0.7), size = 1))
ids <- c(16, 56, 136)

p_btru <- readRDS(file.path(results_wd, "scen_II", 
                            "preds_btru.rds"))$b100.rds$predictions
p_best1 <- readRDS(file.path(results_wd, "scen_II", 
                             "preds_best1.rds"))$b100.rds$predictions
p_best99 <- readRDS(file.path(results_wd, "scen_II", 
                              "preds_best99.rds"))$b100.rds$predictions
p_jmb <- readRDS(file.path(results_wd, "scen_II", 
                           "preds_jmb.rds"))$jmb100.rds$predictions

# Compare the fitted values to the observed values
mean_dat <- d_sim$data_full %>%
  select(id, marker, obstime, mu) %>%
  left_join(d_sim$data %>% select(id, marker, obstime, y), 
            by = c("id", "marker", "obstime")) %>%
  cbind(p_btru$mu_long %>% mutate("TRUE" = Mean) %>%
          select("TRUE")) %>%
  cbind(p_best1$mu_long %>% mutate(EST = Mean) %>% 
          select(EST)) %>%
  cbind(p_best99$mu_long %>% mutate(TRUNC = Mean) %>%
          select(TRUNC)) %>%
  cbind(p_jmb$mu_long %>% mutate(JMB = Mean) %>%
          select(JMB)) %>%
  pivot_longer(c(4, 6:ncol(.))) %>%
  mutate(name = factor(name, levels = c("mu", "TRUE", "EST", "TRUNC", "JMB")))

ggplot(mean_dat %>% filter(id %in% ids) %>%
         mutate(marker = factor(marker, labels = paste("Outcome", 1:2)),
                id = factor(id, labels = c("Subject 16", "Subject 56",
                                           "Subject 136"))),
       aes(x = obstime, y = value, color = name)) +
  facet_grid(id ~ marker, scales = "free_y")+
  theme_bw() +
  geom_point(aes(y = y), color = "grey") +
  geom_line() +
  scale_x_continuous(limits = c(0, 1), labels = c(0, 0.5, 1),
                     breaks = c(0, 0.5, 1)) +
  scale_color_manual(values = c("black", scales::hue_pal()(4))) +
  labs(y = expression(mu(t)~","~hat(mu)(t)), linetype = NULL, color = NULL,
       x = "Time")
# save 4x8


# Plot baseline hazard ----------------------------------------------------

p_btru <- readRDS(file.path(results_wd, "scen_II", 
                            "preds_btru.rds"))
lambga_dat <- lapply(seq_along(p_btru), function (x) {
  d_rirs <- readRDS(file.path(results_wd, "scen_II", "data", 
                              sub("b", "d", names(p_btru)[x])))
  max_ids <- d_rirs$data_full %>% 
    group_by(x3) %>%
    slice(which.max(obstime)) %>% 
    .$id %>% 
    as.vector
  dat <- d_rirs$data_full %>% 
    select(id, obstime, x3) %>% 
    cbind("lambga" = p_btru[[x]]$simulations$lambga_long) %>%
    cbind("est" = p_btru[[x]]$predictions$lambga_long$Mean) %>%
    filter(id %in% max_ids)
  dat
})
lambga_dat <- 
  do.call(rbind, Map(cbind, it = seq_along(lambga_dat), lambga_dat)) %>%
  mutate(itid = interaction(it, id))

ggplot(lambga_dat, aes(x = obstime, group = itid)) + 
  geom_line(aes(y = est), col = "red", alpha = 0.1) +
  geom_line(aes(y = lambga)) + 
  facet_grid(~x3) +
  theme_bw()
