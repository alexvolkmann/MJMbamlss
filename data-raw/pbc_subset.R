## code to prepare `pbc_subset` dataset goes here
library(JMbayes2)
library(tidyverse)
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

# Subset for examples
pbc_subset <- p_long %>%
  filter(id %in% take) %>%
  filter(id %in% head(levels(id), n = 15))%>%
  droplevels()

usethis::use_data(pbc_subset, overwrite = TRUE)
