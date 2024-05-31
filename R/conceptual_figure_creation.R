library(tidyverse)
library(synchrony)
library(compositions)
library(patchwork)
library(gganimate)
library(magick)
library(ggtext)
library(vegan)

synchrony <- function(df){
  species.sd = apply(df, MARGIN = 2, FUN = sd)
  community.var = var(rowSums(df))
  return(
    tibble(denom = round((sum(species.sd)^2)),
           num = round(community.var),
           sync = community.var/(sum(species.sd, na.rm = TRUE)^2))
  )
}

# metacommunity variability
sim_sync <- function(n = 15, cov = 0.8, return.df = F,
                     m1 = 30, m2 = 30, m3 = 30,
                     rho1 = 0.5,
                     v3 = 1){
  # rho1 <- 0.5
  # rho2 <- 0.5
  # rho3 <- 0.5
  
  S <- matrix(c(1.0,  cov,  cov,
                cov,  1.0,  cov,
                cov,  cov,  v3)
              ,nrow=3, ncol=3)
  mu <- c(m1, m2, m3)
  eps <- rlnorm.rplus(n,log(mu),S)
  
  
  df <- tibble(
    time = 1:n,
    A = as.vector(arima.sim(list(ar=rho1),
                            n,innov=eps[,1],
                            start.innov=eps[,1])),
    
    B = as.vector(arima.sim(list(ar=rho1),
                            n,innov=eps[,2],
                            start.innov=eps[,2])),
    
    C = as.vector(arima.sim(list(ar=rho1),
                            n,innov=eps[,3],
                            start.innov=eps[,3]))
  ) %>%
    tidyr::gather(var, value, -time)
  
  return(df)
  
}


set.seed <- 123
output2 <- NULL
sim_series2 <- NULL
for (i in 1:1000){
  t <- sim_sync(n = 10, 
                cov = -0.499, 
                v3 = 1,
                return.df = T,
                m1 = 80, 
                m2 = 80, 
                m3 = 80) %>% 
    filter(time <= 5) %>% 
    dplyr::rename(species = var)
  
  
  evenness <- t %>% 
    group_by(time) %>% 
    dplyr::summarise(evenness = diversity(value, "shannon")/log(specnumber(value))) %>% 
    pull(evenness) %>% 
    mean
  
  sum_of_sds_sqrt <- t %>% 
    group_by(species) %>% 
    dplyr::summarise(num = sd(value)^2) %>% 
    ungroup() %>% 
    dplyr::summarise(sum_sd = sqrt(sum(num))) %>% 
    pull(sum_sd)
  
  sd_of_sum <- t %>% 
    group_by(time) %>% 
    dplyr::summarise(est = sum(value)) %>% 
    ungroup() %>% 
    dplyr::summarise(sd = sd(est)) %>% 
    pull(sd)
  
  # CPE
  cpe_bio <- sum_of_sds_sqrt/sd_of_sum
  
  sum_of_sds <- t %>% 
    group_by(species) %>% 
    dplyr::summarise(num = sd(value)) %>% 
    ungroup() %>% 
    dplyr::summarise(sum_sd = sum(num)) %>% 
    pull(sum_sd)
  
  sae_bio <- sum_of_sds/sum_of_sds_sqrt
  
  
  assign("output2", rbind(output2, 
                         tibble(sae_bio,
                                cpe_bio,
                                async = sae_bio * cpe_bio,
                                evenness,
                                sim = i)))
  
  assign("sim_series2", rbind(t %>% 
                               mutate(sim = i),
                             sim_series2))
}

# high async
high_async <- output2 %>% 
  filter(cpe_bio == max(cpe_bio))

conc_plt <- 
  sim_series2 %>% 
  filter(sim == high_async$sim) %>% 
  mutate(Month = ifelse(time == 1, 3,
                        ifelse(time == 2,
                               5,
                               ifelse(time == 3,
                                      7,
                                      ifelse(time == 4,
                                             9,
                                             ifelse(time == 5,
                                                    11, NA)))))) %>% 
  ggplot() +
  geom_line(aes(y = value, x = Month, color = species),
            linewidth = 1) + 
  geom_point(aes(y = value, x = Month, color = species)) + 
  ggsci::scale_color_d3() +
  scale_x_continuous(breaks = c(3,5,7,9,11)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom")
conc_plt
ggsave(conc_plt,
       file = here::here("figs/conceptual_fig.svg"),
       width = 3.5,
       height = 2.2)

sim_series2 %>% 
  filter(sim == high_async$sim) %>% 
  group_by(time) %>% 
  dplyr::summarise(total = sum(value)) %>% 
  ggplot() +
    geom_line(aes(y = total, x = time)) + 
    scale_y_continuous(limits = c(700, 800))

high_stab <- sim_series2 %>% 
  filter(sim == high_async$sim) %>% 
  group_by(time) %>% 
  dplyr::summarise(total = sum(value)) 


set.seed <- 123
output <- NULL
sim_series <- NULL
for (i in 1:100){
  t <- sim_sync(n = 10, 
                cov = 0.7, 
                v3 = 1,
                return.df = T,
                m1 = 80, 
                m2 = 80, 
                m3 = 80) %>% 
    filter(time <= 5) %>% 
    dplyr::rename(species = var)
  
  
  evenness <- t %>% 
    group_by(time) %>% 
    dplyr::summarise(evenness = diversity(value, "shannon")/log(specnumber(value))) %>% 
    pull(evenness) %>% 
    mean
  
  sum_of_sds_sqrt <- t %>% 
    group_by(species) %>% 
    dplyr::summarise(num = sd(value)^2) %>% 
    ungroup() %>% 
    dplyr::summarise(sum_sd = sqrt(sum(num))) %>% 
    pull(sum_sd)
  
  sd_of_sum <- t %>% 
    group_by(time) %>% 
    dplyr::summarise(est = sum(value)) %>% 
    ungroup() %>% 
    dplyr::summarise(sd = sd(est)) %>% 
    pull(sd)
  
  # CPE
  cpe_bio <- sum_of_sds_sqrt/sd_of_sum
  
  sum_of_sds <- t %>% 
    group_by(species) %>% 
    dplyr::summarise(num = sd(value)) %>% 
    ungroup() %>% 
    dplyr::summarise(sum_sd = sum(num)) %>% 
    pull(sum_sd)
  
  sae_bio <- sum_of_sds/sum_of_sds_sqrt
  
  
  assign("output", rbind(output, 
                          tibble(sae_bio,
                                 cpe_bio,
                                 async = sae_bio * cpe_bio,
                                 evenness,
                                 sim = i)))
  
  assign("sim_series", rbind(t %>% 
                                mutate(sim = i),
                              sim_series))
}

# high async
low_async <- output %>% 
  filter(cpe_bio == min(cpe_bio))

sync_plt <- 
  sim_series %>% 
  filter(sim == low_async$sim) %>% 
  mutate(Month = ifelse(time == 1, 3,
                        ifelse(time == 2,
                               5,
                               ifelse(time == 3,
                                      7,
                                      ifelse(time == 4,
                                             9,
                                             ifelse(time == 5,
                                                    11, NA)))))) %>% 
  ggplot() +
  geom_line(aes(y = value, x = Month, color = species),
            linewidth = 1) + 
  geom_point(aes(y = value, x = Month, color = species)) + 
  ggsci::scale_color_d3() +
  scale_x_continuous(breaks = c(3,5,7,9,11)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom")
sync_plt


ggsave(sync_plt,
       file = here::here("figs/conceptual_fig1.svg"),
       width = 3.5,
       height = 2.2)

low_stab <- 
  sim_series %>% 
  filter(sim == low_async$sim) %>% 
  group_by(time) %>% 
  dplyr::summarise(total = sum(value)) %>% 
  mutate(Month = ifelse(time == 1, 3,
                        ifelse(time == 2,
                               5,
                               ifelse(time == 3,
                                      7,
                                      ifelse(time == 4,
                                             9,
                                             ifelse(time == 5,
                                                    11, NA))))))


comb_stab <- high_stab %>% 
  mutate(Month = ifelse(time == 1, 3,
                        ifelse(time == 2,
                               5,
                               ifelse(time == 3,
                                      7,
                                      ifelse(time == 4,
                                             9,
                                             ifelse(time == 5,
                                                    11, NA)))))) %>% 
  ggplot() +
  geom_point(aes(y = total, x = Month)) +
  geom_line(aes(y = total, x = Month),
            linewidth = 1) +
  geom_point(data = low_stab, aes(y = total, x = Month)) +
  geom_line(data = low_stab, aes(y = total, x = Month),
            linewidth = 1) +
  scale_x_continuous(breaks = c(3,5,7,9,11)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom") + 
  scale_y_continuous(limits = c(0, 1000))

ggsave(comb_stab,
       file = here::here("figs/conceptual_fig_comb_stab.svg"),
       width = 3.5,
       height = 2.2)

low_stab %>% 
  ggplot() +
  geom_point(aes(y = total, x = time)) +
  geom_line(aes(y = total, x = time),
            linewidth = 1) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom") + 
  scale_y_continuous(limits = c(0, 3000))
