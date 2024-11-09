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
                     m1 = 30, m2 = 30,
                     rho1 = 0,
                     v3 = 1){
  # rho1 <- 0.5
  # rho2 <- 0.5
  # rho3 <- 0.5
  
  S <- matrix(c(1.0,  cov,
                cov,  1.0)
              ,nrow=2, ncol=2)
  mu <- c(m1, m2)
  eps <- rlnorm.rplus(n,log(mu),S)
  
  
  df <- tibble(
    time = 1:n,
    A = as.vector(arima.sim(list(ar=rho1),
                            n,innov=eps[,1],
                            start.innov=eps[,1])),
    
    B = as.vector(arima.sim(list(ar=rho1),
                            n,innov=eps[,2],
                            start.innov=eps[,2]))
  ) %>%
    tidyr::gather(var, value, -time)
  
  return(df)
  
}


set.seed <- 123
output2 <- NULL
sim_series2 <- NULL
for (i in 1:1000){
  t <- sim_sync(n = 10, 
                cov = -0.49, 
                # v3 = 1,
                return.df = T,
                m1 = 80, 
                m2 = 80) %>% 
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

async_df1 <-
sim_series2 %>% 
  filter(sim == high_async$sim)

async_plt <-
  async_df1 %>%
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
  labs(y = "Harvests",
       color = "Species") +
  ggsci::scale_color_cosmic() +
  scale_x_continuous(breaks = c(3,5,7,9,11)) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        legend.position  = c(0.5, 0.86),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "transparent"))

async_plt
ggsave(async_plt,
       file = here::here("async_figure/async_plt1.svg"),
       dpi = 300,
       width = 3,
       height = 2.2)

async1_total <- 
  sim_series2 %>% 
    filter(sim == high_async$sim) %>% 
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
async_df2 <- 
  sim_series2 %>% 
    filter(sim == high_async$sim) %>% 
    mutate(value = ifelse(time %in% 1:2 & species == "A", 0, value))


async_plt2 <-
  async_df2 %>%
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
  labs(y = "Harvests") +
  guides(color = "none")+
  ggsci::scale_color_cosmic() +
  scale_x_continuous(breaks = c(3,5,7,9,11)) +
  theme_bw() +
  theme(axis.ticks = element_blank())

async_plt2
ggsave(async_plt2,
       file = here::here("async_figure/async_plt2.svg"),
       dpi = 300,
       width = 3,
       height = 2.2)

async2_total <- 
  async_df2 %>% 
  group_by(time) %>% 
  dplyr::summarise(total = sum(value))%>%
  mutate(Month = ifelse(time == 1, 3,
                        ifelse(time == 2,
                               5,
                               ifelse(time == 3,
                                      7,
                                      ifelse(time == 4,
                                             9,
                                             ifelse(time == 5,
                                                    11, NA))))))

async_df3 <-
  sim_series2 %>% 
  filter(sim == high_async$sim) %>% 
  group_by(species) %>% 
  mutate(value = ifelse(species == "B", lead(value, 1), value),
         value = ifelse(species == "B" & time == 5, rnorm(n = 1,
                                                          mean = 80, 
                                                          sd = 5), value))

async_plt3 <-
  async_df3 %>%
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
  labs(y = "Harvests") +
  guides(color = "none")+
  ggsci::scale_color_cosmic() +
  scale_x_continuous(breaks = c(3,5,7,9,11)) +
  theme_bw() +
  theme(axis.ticks = element_blank())

async_plt3
ggsave(async_plt3,
       file = here::here("async_figure/async_plt3.svg"),
       dpi = 300,
       width = 3,
       height = 2.2)

async3_total <- 
  async_df3 %>% 
  group_by(time) %>% 
  dplyr::summarise(total = sum(value))%>%
  mutate(Month = ifelse(time == 1, 3,
                        ifelse(time == 2,
                               5,
                               ifelse(time == 3,
                                      7,
                                      ifelse(time == 4,
                                             9,
                                             ifelse(time == 5,
                                                    11, NA))))))

portfolio_plt <- 
  ggplot() +
  geom_line(data = async1_total, 
            aes(y = total, x = Month),
            linewidth = 1) + 
  geom_point(data = async1_total, 
            aes(y = total, x = Month)) + 
  
  geom_line(data = async2_total, 
            aes(y = total, x = Month)) + 
  geom_point(data = async2_total, 
             aes(y = total, x = Month)) + 
  
  geom_line(data = async3_total, 
            aes(y = total, x = Month)) +
  geom_point(data = async3_total, 
             aes(y = total, x = Month)) + 
  
  labs(y = "Harvests") +
  
  scale_y_continuous(limits = c(400, 950)) +
  scale_x_continuous(breaks = c(3,5,7,9,11)) +
  theme_bw() 
  
mean(async1_total$total)/sd(async1_total$total)
mean(async2_total$total)/sd(async2_total$total)
mean(async3_total$total)/sd(async3_total$total)

ggsave(portfolio_plt,
       file = here::here("async_figure/portfolio_plt.svg"),
       dpi = 300,
       width = 3,
       height = 2.2)
