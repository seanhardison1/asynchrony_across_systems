library(tidyverse)
library(lubridate)
library(tsibble)
library(magrittr)

for (i in c("ac",
            "sb","wp",
            "sp",
            # "sh",
            "gz","cc")){
  load(here::here(paste0("data/model_output/pred_out_seas_",i,"_sim.rdata")))
}


pred_out_ann <- bind_rows(pred_out_seas_ac,
                           pred_out_seas_sb,
                           pred_out_seas_sp,
                           pred_out_seas_wp,
                          # pred_out_seas_sh,
                          pred_out_seas_gz,
                          pred_out_seas_cc) %>% 
  
  # filter(common %in% selected_specs) %>% 
  {. ->> pred_out_seas_complete} %>% 

  mutate(month = month(yearmonth(ymon))) %>% 
  as_tibble() %>% 
  {. ->> pred_out_seas} %>% 
  group_by(metacomm, year, common) %>% 
  dplyr::summarise(est = sum(est)) %>% 
  as_tibble() 

ggplot(pred_out_ann) +
  geom_area(aes(y = est, x = year, fill = common)) +
  facet_wrap(~metacomm)

save(pred_out_ann, 
     pred_out_seas, 
     file = here::here("data/mod_results_423.rdata"))
# file = here::here("data/mod_results_117.rdata"))
