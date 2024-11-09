library(tidyverse)
library(readxl)

df <- read_excel(here::here('async_figure/conc_fig_data.xlsx'))

tot_plt <- 
  df %>% 
  group_by(var, time) %>% 
  dplyr::summarise(total = sum(value)) %>% 
  ggplot() +
    geom_line(aes(y = total, x = time, color = var), show.legend = F) + 
    theme_bw() +
    scale_x_continuous(breaks = c(3,5,7,9,11))

ggsave(tot_plt,
       file = here::here("async_figure/tot_plt.svg"),
       dpi = 300,
       width = 3,
       height = 2.2)


stab_plt <- 
  df %>% 
  group_by(var, time) %>% 
  dplyr::summarise(total = sum(value)) %>% 
  group_by(var) %>% 
  dplyr::summarise(stability = 1/(sd(total)/mean(total))) %>% 
  ggplot() + 
    geom_bar(aes(y = stability, x = var, fill = var), show.legend = F,
             stat = "identity") + 
  theme_bw()

ggsave(stab_plt,
       file = here::here("async_figure/stan_plt.svg"),
       dpi = 300,
       width = 3,
       height = 2.2)
