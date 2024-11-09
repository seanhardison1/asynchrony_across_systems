library(tsibble)
library(lubridate)
library(tidyverse)
library(sf)
library(magrittr)
library(ggeffects)
# library(synchrony)
# library(vegan)

# biomass indices----
load(here::here("data/mod_results_423.rdata"))
# load(here::here("data/mod_results_117.rdata"))


# processed landings----
load(here::here("data/processed_landings.rdata"))

# stability-synchrony functions (Lamy et al. 2018)----
source(here::here("R/stability_calculations/part_stab_functions2.r"))

va_specs <- c("Atlantic croaker",
              "spot",
              "striped bass")

md_specs <- c("Atlantic croaker",
              "striped bass",
              "white perch",
              "gizzard shad",
              "channel catfish")


total_revenue_annual <- all_landings_complete %>% 
  filter((species %in% va_specs & metacomm == "VA") |
           (species %in% md_specs & metacomm == "MD")) %>%
  group_by(metacomm, year) %>% 
  dplyr::summarise(total_revenue_annual = sum(total_value),
                   total_landings_annual = sum(total_landings))

# asynchrony-stability partitioning----
sync_out <- tibble()
s_partition <- tibble()
month_vec2 <- 1:12
for (m in c("MD","VA")){
  for (y in 2002:2018){
    
    
    if (m == "VA"){
      specs <- va_specs
    } else if (m == "MD"){
      specs <- md_specs
    }
    
    # COMMUNITY BIOMASS-----
    int <- pred_out_seas %>% 
      dplyr::select(metacomm, year, ymon, est, common, month) %>% 
      filter(year == y,
             metacomm == m,
             common %in% specs) %>%
      dplyr::select(-month) %>% 
      pivot_wider(names_from = "common",
                  values_from = "est",
                  values_fill = 0) %>% 
      ungroup() %>% 
      dplyr::select(-1, -2, -3) %>% 
      dplyr::select(-which(colSums(.) == 0))
    
    ## stability/asynchrony partitioning -----
    int2 <- part_stab_comp(Y = int, 
                           s = 1, t = length(month_vec))
    
    ### synchrony and variability
    phi_SC_k <- int2$patch$phi_SC_k 
    cv_cr <- int2$part$CV_CR # metacommunity variability
    cv_sr <- int2$part$CV_SR # metapopulation variability
    cv_sr_part <- int2$CV_SR_par %>% # partitioned metapopulation variability
      mutate(year = y,
             metacomm = m,
             common = row.names(.))
    phi_SC_k <- int2$patch$phi_SC_k # species synchrony

    ## asynchrony partitioning (adapted Zhao et al. 2022)----
    t <- int %>% 
      mutate(time = c(3,5,7,9,11)) %>% 
      gather(species, value, -time)
    
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
    
    ### CPE BIO----
    cpe_bio <- sum_of_sds_sqrt/sd_of_sum
    
    sum_of_sds <- t %>% 
      group_by(species) %>% 
      dplyr::summarise(num = sd(value)) %>% 
      ungroup() %>% 
      dplyr::summarise(sum_sd = sum(num)) %>% 
      pull(sum_sd)
    
    ### SAE BIO-----
    sae_bio <- sum_of_sds/sum_of_sds_sqrt
    

    # HARVEST VALUE STABILITY (all months)-----
    int_val_sync <- all_landings %>% 
      filter(year == y,
             metacomm == m,
             species %in% specs) %>% 
      dplyr::select(-total_landings, -total_trips, -metacomm) %>% 
      right_join(.,tibble(month = month_vec2,
                          year = y,
                          metacomm = m) %>% 
                   expand_grid(species = specs)) %>%
      mutate(total_value = ifelse(is.na(total_value), 0, total_value)) %>% 
      pivot_wider(names_from = "species",
                  values_from = "total_value",
                  values_fill = 0) %>% 
      ungroup() %>% 
      arrange(year, month) %>% 
      dplyr::select(-1, -2, -3) %>% 
      dplyr::select(-which(colSums(.) == 0))
    
    ## harvest value stability/asynchrony partitioning -----
    int_val_sync2 <- part_stab_comp(Y = int_val_sync, 
                                    s = 1, t = length(month_vec2))
    cv_cr_val <- int_val_sync2$part$CV_CR
    cv_sr_val <- int_val_sync2$part$CV_SR
    
    # HARVEST LANDINGS STABILITY (subset of months)-----
    int_land_sync_short <- all_landings %>% 
      filter(year == y,
             metacomm == m,
             species %in% specs) %>% 
      dplyr::select(-total_trips, -total_value, -metacomm) %>%
      right_join(.,tibble(month = month_vec,
                          year = y,
                          metacomm = m) %>% 
                   expand_grid(species = specs)) %>%
      mutate(total_landings = ifelse(is.na(total_landings), 0, total_landings)) %>% 
      pivot_wider(names_from = "species",
                  values_from = "total_landings",
                  values_fill = 0) %>% 
      ungroup() %>% 
      arrange(year, month) %>% 
      dplyr::select(-1, -2, -3) %>% 
      dplyr::select(-which(colSums(.) == 0))
    
    ## stability/asynchrony partitioning----
    int_land_sync_short2 <- part_stab_comp(Y = int_land_sync_short, 
                                           s = 1, t = length(month_vec))
    cv_cr_land_short <- int_land_sync_short2$part$CV_CR
    sync_land_short <- int_land_sync_short2$patch$phi_SC_k
    cv_sr_land_short <- int_land_sync_short2$part$CV_SR
    
    cv_sr_land_part <- int_land_sync_short2$CV_SR_part %>%
      mutate(year = y,
             metacomm = m,
             common = row.names(.)) %>%
      dplyr::rename(CV_SR_part_land = CV_SR_part,
                    SD_SR_part_land = SD_SR_part,
                    weight_land = weight)
    
    ## asynchrony partitioning (adapted Zhao et al. 2022)----
    t_land <- int_land_sync_short %>% 
      mutate(time = c(3,5,7,9,11)) %>% 
      gather(species, value, -time)
    
    sum_of_sds_sqrt <- t_land %>% 
      group_by(species) %>% 
      dplyr::summarise(num = sd(value)^2) %>% 
      ungroup() %>% 
      dplyr::summarise(sum_sd = sqrt(sum(num))) %>% 
      pull(sum_sd)
    
    sd_of_sum <- t_land %>% 
      group_by(time) %>% 
      dplyr::summarise(est = sum(value)) %>% 
      ungroup() %>% 
      dplyr::summarise(sd = sd(est)) %>% 
      pull(sd)
    
    ### CPE LAND------
    cpe_land <- sum_of_sds_sqrt/sd_of_sum
    
    sum_of_sds <- t_land %>% 
      group_by(species) %>% 
      dplyr::summarise(num = sd(value)) %>% 
      ungroup() %>% 
      dplyr::summarise(sum_sd = sum(num)) %>% 
      pull(sum_sd)
    
    ### SAE LAND------
    sae_land <- sum_of_sds/sum_of_sds_sqrt
    
    # HARVEST LANDINGS STABILITY (all months)-----
    int_land_sync <- all_landings %>% 
      filter(year == y,
             metacomm == m,
             species %in% specs) %>% 
      dplyr::select(-total_trips, -total_value, -metacomm) %>% 
      right_join(.,tibble(month = month_vec2,
                          year = y,
                          metacomm = m) %>% 
                   expand_grid(species = specs)) %>%
      mutate(total_landings = ifelse(is.na(total_landings), 0, total_landings)) %>%
      pivot_wider(names_from = "species",
                  values_from = "total_landings",
                  values_fill = 0) %>% 
      ungroup() %>% 
      arrange(year, month) %>% 
      dplyr::select(-1, -2, -3) %>% 
      dplyr::select(-which(colSums(.) == 0))
    
    ## stability/asynchrony partitioning-----
    int_land_sync2 <- part_stab_comp(Y = int_land_sync, 
                                     s = 1, t = length(month_vec2))
    sync_land_long <- int_land_sync2$patch$phi_SC_k
    cv_sr_spec_long <- int_land_sync2$part$CV_SR
    cv_cr_land_long <- int_land_sync2$part$CV_CR
   
    # cv_sr_land_part <- int_land_sync2$CV_SR_part %>%
    #   mutate(year = y,
    #          metacomm = m,
    #          common = row.names(.)) %>%
    #   dplyr::rename(CV_SR_part_land = CV_SR_part,
    #                 SD_SR_part_land = SD_SR_part,
    #                 weight_land = weight)
    
    assign('s_partition', rbind(cv_sr_part %>% 
                                  left_join(.,cv_sr_land_part,
                                            by = c("common",
                                                   "metacomm",
                                                   "year")), s_partition))
    assign("sync_out", rbind(tibble(sync_land_short,
                                    sync_land_long,
                                    cv_sr_spec_long,
                                    cv_cr_land_long,
                                    phi_SC_k,
                                    cv_cr,
                                    cv_sr,
                                    cv_cr_land_short,
                                    cv_sr_land_short,
                                    cv_cr_val,
                                    cpe_bio,
                                    sae_bio,
                                    cpe_land,
                                    sae_land,
                                    cv_sr_val,
                                    metacomm = m,
                                    year = y),
                             sync_out))
    
  }
}

# Atlantic croaker biomass----
ac_bio <- pred_out_seas %>% 
  dplyr::select(metacomm, year, ymon, est, common, month) %>%
  filter(common %in% c("Atlantic croaker")) %>%
  group_by(metacomm, year) %>% 
  dplyr::summarise(ac_bio = sum(est))


va_df <- 
  sync_out %>% 
  filter(metacomm == "VA") %>% 
  left_join(.,total_revenue_annual %>% 
              filter(metacomm == "VA")) %>% 
  left_join(.,ac_bio %>% 
              filter(metacomm == "VA")) %>% 
  mutate(
    total_revenue_annual2 = as.numeric(scale(total_revenue_annual)),#a
    total_landings_annual2 = as.numeric(scale(total_landings_annual)),
    pe_land = 1/sync_land_short,
    pe_land_long = 1/sync_land_long,
    s_sr = 1/cv_sr,
    s_land = 1/cv_cr_land_long,
    s_land_short = 1/cv_cr_land_short,
    s_land_spec_short = 1/cv_sr_land_short,
    s_land_spec_long = 1/cv_sr_spec_long,
    species_pe = 1/phi_SC_k,
    s_bio = 1/cv_cr,
    s_bio_spec = 1/cv_sr,
    s_val = 1/cv_cr_val,
    s_val_spec = 1/cv_sr_val
  )


md_df <- 
  sync_out %>% 
  filter(metacomm == "MD") %>% 
  left_join(.,total_revenue_annual %>% 
              filter(metacomm == "MD")) %>% 
  mutate(
    total_revenue_annual2 = as.numeric(scale(total_revenue_annual)),#a
    total_landings_annual2 = as.numeric(scale(total_landings_annual)),
  )%>% 
  mutate(
    total_revenue_annual2 = as.numeric(scale(total_revenue_annual)),#a
    total_landings_annual2 = as.numeric(scale(total_landings_annual)),
    pe_land = 1/sync_land_short,
    pe_land_long = 1/sync_land_long,
    s_sr = 1/cv_sr,
    s_land = 1/cv_cr_land_long,
    s_land_short = 1/cv_cr_land_short,
    s_land_spec_short = 1/cv_sr_land_short,
    s_land_spec_long = 1/cv_sr_spec_long,
    species_pe = 1/phi_SC_k,
    s_bio = 1/cv_cr,
    s_bio_spec = 1/cv_sr,
    s_val = 1/cv_cr_val,
    s_val_spec = 1/cv_sr_val
  ) %>% 
  left_join(.,all_landings_complete %>% 
              filter(metacomm == "MD",
                     species %in% md_specs,
                     month %in% month_vec) %>%
              mutate(month_grp = ifelse(month %in% c(3:5),
                                        "spring","summer/fall")) %>% 
              group_by(year, month_grp, month) %>%
              dplyr::summarise(tot_land = sum(total_landings, na.rm = T)) %>% 
              group_by(year, month_grp) %>% 
              dplyr::summarise(mean_land = mean(tot_land, 
                                                na.rm = T)) %>% 
              spread(.,month_grp, mean_land) %>% 
              mutate(land_ratio = spring/`summer/fall`))

