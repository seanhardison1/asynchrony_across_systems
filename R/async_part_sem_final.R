library(piecewiseSEM)
library(nlme)
library(gratia)
library(patchwork)
library(mgcv)
library(ggtext)
library(tidyverse)
library(pals)

source(here::here("data-raw/landings_processing.R"))
source(here::here("R/stability_calculations/stability_final.R"))

ndf_va <- va_df %>% 
  mutate(species_pe2 = log(species_pe),
         s_land_short2 = log(s_land_short),
         s_land_spec_short2 = log(s_land_spec_short),
         s_bio_spec2 = log(s_bio_spec),
         s_bio2 = log(s_bio),
         s_land2 = log(s_land),
         pe_land2 = log(pe_land),
         cpe_bio2 = log(cpe_bio),
         sae_bio2 = log(sae_bio),
         s_val2 = log(s_val),
         cpe_land2 = log(cpe_land),
         sae_land2 = log(sae_land)) %>% 
  left_join(.,all_landings %>%
              filter((species %in% va_specs & metacomm == "VA"),
                     month %in% c(3,5,7,9,11)) %>%
              
              # filter((total_landings > 0 & total_trips > 0)) %>% 
              # ^total trips and total landings must be > 0 - if total_trips > 0 and total_landings == 0,
              # then it's a false zero due to privacy concerns
              
              dplyr::select(month, species, total_landings,
                            total_trips, year, metacomm) %>% 
              group_by(year, species) %>% 
              dplyr::summarise(total_trips = sum(total_trips, na.rm = T)) %>% 
              group_by(year) %>% 
              dplyr::summarise(mean_trips = mean(total_trips, na.rm = T))) %>% 
  left_join(.,all_landings %>%
              filter((species %in% va_specs & metacomm == "VA"),
                     month %in% c(3,5,7,9,11)) %>%
              
              # filter((total_landings > 0 & total_trips > 0)) %>% 
              # ^total trips and total landings must be > 0 - if total_trips > 0 and total_landings == 0,
              # then it's a false zero due to privacy concerns
              
              dplyr::select(month, species, total_landings,
                            total_trips, year, metacomm) %>% 
              group_by(year, species) %>% 
              dplyr::summarise(total_trips = sum(total_trips, na.rm = T)) %>% 
              spread(species, total_trips) %>% 
              dplyr::rename(ac_effort = 2,
                            sp_effort = 3,
                            sb_effort = 4)) #%>% 
  # dplyr::select(cpe_land, cpe_bio, sae_land, sae_bio, 
  #               s_land_spec, mean_trips, pe_land,
  #               s_land)
# VA 1 -------
m1 <- lm(cpe_land2 ~ cpe_bio2 +
           sae_land2 +
           sae_bio2, data = ndf_va)
s <- DHARMa::simulateResiduals(m1, n = 1000);plot(s)
performance::check_collinearity(m1)
acf(residuals(m1))
summary(m1)

m2 <- gls(cpe_bio2 ~ sae_bio2, data = ndf_va,
          correlation = corAR1())
acf(residuals(m2, type = "normalized"))
summary(m2)

m3 <- gls(sae_land2 ~ sae_bio2, data = ndf_va, 
          correlation = corAR1())
acf(residuals(m3, type = "normalized"))
summary(m3)

m4 <- gls(mean_trips ~ sae_bio2, data = ndf_va)
acf(residuals(m4, type = "normalized"))
summary(m4)
 
y1 <- lm(s_land_spec_short2 ~ mean_trips, data = ndf_va)
s <- DHARMa::simulateResiduals(y1, n = 1000);plot(s)
acf(residuals(y1))
summary(y1)

# y2 <- lm(mean_trips ~ sae_bio2, data = ndf_va)
# s <- DHARMa::simulateResiduals(y2, n = 1000);plot(s)
# acf(residuals(y2))
# summary(y2)

p1 <- psem(m1, 
           m2,
           m3,
           m4,
           y1)
p2 <- update(p1, cpe_bio2 %~~% s_land_spec_short2)
# p3 <- update(p2, mean_trips %~~% sae_bio2)
summary(p2)
plot(p2)


# VA 2 --------
m4 <- lm(s_land2 ~ s_land_short2, data = ndf_va)
s <- DHARMa::simulateResiduals(m4, n = 1000);plot(s)
acf(residuals(m4))
summary(m4)

m5 <- lm(total_landings_annual2 ~ s_land2, data = ndf_va)
s <- DHARMa::simulateResiduals(m5, n = 1000);plot(s)
acf(residuals(m5))
summary(m5)

m6 <- lm(s_val2 ~ s_land2,
          data = ndf_va)
s <- DHARMa::simulateResiduals(m6, n = 1000);plot(s)
acf(residuals(m6))
summary(m6)

m7 <- lm(total_revenue_annual2 ~ s_val2 +
           total_landings_annual2,
         data = ndf_va)
s <- DHARMa::simulateResiduals(m7, n = 1000);plot(s)
acf(residuals(m7))
summary(m7)

p2 <- psem(m4, m5, m6, m7)
summary(p2)
plot(p2)

ndf_md <- md_df %>% 
  mutate(species_pe2 = log(species_pe),
         s_land_short2 = log(s_land_short),
         s_land_spec_short2 = log(s_land_spec_short),
         s_land2 = log(s_land),
         pe_land2 = log(pe_land),
         cpe_bio2 = log(cpe_bio),
         sae_bio2 = log(sae_bio),
         cpe_land2 = log(cpe_land),
         sae_land2 = log(sae_land),
         s_val2 = log(s_val),
         s_land_short_sqrd = s_land_short^2)


ggplot(ndf_md) + 
  geom_point(aes(y = cpe_land2, x = year))

ggplot(ndf_md) + 
  geom_point(aes(y = cpe_land2, x = sae_bio2))


ggplot(pred_out_ann) + 
  geom_line(aes(y = est, x = year, color = common)) +
  facet_wrap(~metacomm)

ggplot(md_landings %>% 
         filter(species %in% md_specs) %>% 
         group_by(year, species) %>% 
         dplyr::summarise(total_landings = sum(total_landings))) + 
  geom_line(aes(y = total_landings, x = year, color = species))

# md 1 -------
m8 <- lm(cpe_land2 ~ cpe_bio2 + sae_land2 + sae_bio2 + 
           land_ratio, data = ndf_md)
performance::check_collinearity(m8)
s <- DHARMa::simulateResiduals(m8, n = 1000);plot(s)
acf(residuals(m8))
summary(m8)

m9 <- lm(cpe_bio2 ~ sae_bio2, data = ndf_md)
s <- DHARMa::simulateResiduals(m9, n = 1000);plot(s)
acf(residuals(m9))
summary(m9)

m10 <- lm(sae_land2 ~ sae_bio2, data = ndf_md)
s <- DHARMa::simulateResiduals(m10, n = 1000);plot(s)
acf(residuals(m10))
summary(m10)

p3 <- psem(m8, m9, m10)
summary(p3)
plot(p3) 

# md 2 --------
m11 <- lm(s_land2 ~ s_land_short2, data = ndf_md)
s <- DHARMa::simulateResiduals(m11, n = 1000);plot(s)
acf(residuals(m11))
summary(m11)

m12 <- lm(total_landings_annual2 ~ s_land2, data = ndf_md)
s <- DHARMa::simulateResiduals(m12, n = 1000);plot(s)
acf(residuals(m12))
summary(m12)

m13 <- lm(s_val2 ~ s_land2,
         data = ndf_md)
s <- DHARMa::simulateResiduals(m13, n = 1000);plot(s)
acf(residuals(m13))
summary(m13)

m14 <- gls(total_revenue_annual2 ~ s_val2 +
           total_landings_annual2,
         data = ndf_md,
         correlation = corAR1())
# s <- DHARMa::simulateResiduals(m13, n = 1000);plot(s)
acf(residuals(m14, type = "normalized"))
summary(m14)

p3 <- psem(m11, m12, m13, m14)
summary(p3)
plot(p3, alpha = 0.05)


# asynchrony in communities contribute to fishery asynchrony but only
# after unevenness in species' population variabilities are accounted for. 
# This suggests that fisheries in VA are stabilized by the compensatory dynamics
# of the species that they target. The CPE_land ~ SAE_bio effect occurs 
# because the evenness of population variabilities increased following the decline
# of biomass dominant species like Atlantic croaker, while at the same time harvest
# seasonal compensation has declined. This could be due to a decline in fishing effort
# for Atlanic croaker as the have declined, with more Atlantic croaker being caught
# later in the year, resulting in greater overlap between Atlantic croaker and spot harvests.
# 


fill_vec <- setNames(ggsci::pal_d3("category20")(n = 6),unique(c(va_specs, 
                                                                 md_specs)))

summary(lm(cpe_land ~ land_ratio, data = ndf_md))

cpe_fig_middle <-
  ggplot(ndf_md) + 
  geom_point(aes(y = cpe_land, x = land_ratio, group = year)) +
    geom_smooth(aes(y = cpe_land, x = land_ratio), 
                method = "lm",
                color = "black") +
  labs(y = "Harvest compensation (CPE<sub>Harvest</sub>)",
       x = "Seasonal harvest ratio (SHR)") +
  theme_bw() +
  theme(axis.title.y = element_markdown(size = 15),
        axis.title.x = element_markdown(size = 15))

cpe_fig_left <- md_landings %>% 
  filter(species %in% md_specs,
         month %in% month_vec) %>% 
  right_join(.,tibble(month = month_vec,
                      metacomm = m) %>% 
               expand_grid(species = md_specs,
                           year = c(2002:2018))) %>% 
  mutate(total_landings = ifelse(is.na(total_landings), 0, total_landings)) %>% 
  filter(year %in% c(2009)) %>% 
ggplot() + 
  geom_rect(aes(ymin = 0, ymax = Inf,
                xmin = 3, xmax = 5),
            fill = "#FFCCCC80",
            color = "#FFCCCC80") +
  labs(y = "Total harvests (lbs)",
       x = "Month") +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     labels = function(x) format(x, scientific = TRUE)) +
  geom_line(aes(y = total_landings, x = month, color = species),
            linewidth = 1) +
  geom_point(aes(y = total_landings, x = month, color = species)) +
  scale_color_manual(values = fill_vec) +
  guides(color = "none") +
  theme_bw() +
  theme(axis.title.y = element_markdown(size = 13),
        axis.title.x = element_markdown(size = 12),
        axis.text = element_markdown(size = 12)) 

cpe_fig_right <-
  md_landings %>% 
  filter(species %in% md_specs,
         month %in% month_vec) %>% 
  right_join(.,tibble(month = month_vec,
                      # year = c(2009, 2012),
                      metacomm = m) %>% 
               expand_grid(species = md_specs,
                           year = c(2002:2018))) %>% 
  mutate(total_landings = ifelse(is.na(total_landings), 0, total_landings)) %>% 
  filter(year %in% c(2012)) %>% 
  ggplot() + 
  geom_rect(aes(ymin = 0, ymax = Inf,
                xmin = 3, xmax = 5),
            fill = "#FFCCCC80",
            color = "#FFCCCC80") +
  labs(y = "Total harvests (lbs)",
       x = "Month") +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     labels = function(x) format(x, scientific = TRUE)) +
  geom_line(aes(y = total_landings, x = month, color = species),
            linewidth = 1) +
  geom_point(aes(y = total_landings, x = month, color = species)) +
  scale_color_manual(values = fill_vec) +
  guides(color = "none") +
  theme_bw() +
  theme(axis.title.y = element_markdown(size = 13),
        axis.title.x = element_markdown(size = 12),
        axis.text = element_markdown(size = 12)) 


cpe_fig_guide <-
  md_landings %>% 
  filter(species %in% md_specs,
         month %in% month_vec) %>% 
  right_join(.,tibble(month = month_vec,
                      # year = c(2009, 2012),
                      metacomm = m) %>% 
               expand_grid(species = md_specs,
                           year = c(2002:2018))) %>% 
  mutate(total_landings = ifelse(is.na(total_landings), 0, total_landings)) %>% 
  filter(year %in% c(2012)) %>% 
  ggplot() + 
  geom_rect(aes(ymin = 0, ymax = Inf,
                xmin = 3, xmax = 5),
            fill = "#FFCCCC80",
            color = "#FFCCCC80") +
  labs(y = "Total harvests (lb)",
       x = "Month") +
  scale_y_continuous(expand = c(0.02, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     labels = function(x) format(x, scientific = TRUE)) +
  geom_line(aes(y = total_landings, x = month, color = species),
            linewidth = 1) +
  geom_point(aes(y = total_landings, x = month, color = species),
             size = 1.5) +
  scale_color_manual(values = fill_vec) +
  # guides(color = "none") +
  theme_bw() +
  theme(axis.title.y = element_markdown(size = 13),
        axis.title.x = element_markdown(size = 12),
        axis.text = element_markdown(size = 12),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_markdown(size = 11)) 


ggsave(cpe_fig_middle, filename = here::here("figs/cpe_fig_middle.png"),
       dpi = 300,
       width = 5, height = 4)


ggsave(cpe_fig_left, filename = here::here("figs/cpe_fig_left.png"),
       dpi = 300,
       width = 5, height = 3.5)


ggsave(cpe_fig_right, filename = here::here("figs/cpe_fig_right.png"),
       dpi = 300,
       width = 5, height = 3.5)


ggsave(cpe_fig_guide, filename = here::here("figs/cpe_fig_guide.png"),
       dpi = 300,
       width = 8, height = 4)


# biomass trends ----
ac <- pred_out_ann %>% filter(common == "Atlantic croaker",
                              metacomm == "VA")

mod_ac <- lm(est ~ year, data = ac)
summary(mod_ac)
acf(residuals(mod_ac))

sp <- pred_out_ann %>% filter(common == "spot",
                              metacomm == "VA")

mod_sp <- lm(est ~ year, data = sp)
summary(mod_sp)
acf(residuals(mod_sp))

sb <- pred_out_ann %>% filter(common == "striped bass",
                              metacomm == "VA")

mod_sb <- lm(est ~ year, data = sb)
summary(mod_sb)
acf(residuals(mod_sb))

biomass_trend_stats <- 
  bind_rows(
    broom::tidy(mod_ac) %>% 
      filter(term == "year") %>% 
      mutate(Species = "Atlantic croaker",
             `Error structure` = "iid"),
    broom::tidy(mod_sp) %>% 
      filter(term == "year") %>% 
      mutate(Species = "spot",
             `Error structure` = "iid"),
    broom::tidy(mod_sb) %>% 
      filter(term == "year") %>% 
      mutate(Species = "striped bass",
             `Error structure` = "iid")
  ) %>% 
    dplyr::select(Species, 
                  Trend = estimate, 
                  `Std. error` = std.error, 
                  `T statistic` = statistic,
                  `P value` = p.value,
                  `Error structure`) %>% 
  mutate(`P value` = round(`P value`, 3))

write.csv(biomass_trend_stats, file = here::here("data/biomass_trend_stats.csv"),
          row.names = F)

# biomass time series----
fill_vec <- setNames(ggsci::pal_d3("category20")(n = 6),unique(c(va_specs, 
                                                                 md_specs)))
annual_biomass <- 
  ggplot(data = pred_out_ann) +
  geom_line(aes(y = est, x = year, color = common),
            linewidth = 1, show.legend = F) +
  labs(y = "Biomass (kg)") +
  scale_color_manual(values = fill_vec) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     labels = function(x) format(x, scientific = TRUE)) +
  facet_wrap(~metacomm) +
  theme_bw() +
  theme(axis.title.y = element_markdown(size = 15),
        strip.background = element_rect(fill = 'transparent', color = NA),
        strip.text = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 14))

ggsave(annual_biomass, file = here::here("figs/annual_biomass_time_series.png"),
       dpi = 300,
       width = 7,
       height = 3)

total_harvests <- 
  bind_rows(
    va_landings %>% 
      filter(species %in% c("Atlantic croaker",
                            "spot",
                            "striped bass")) %>% 
      group_by(species, year) %>% 
      dplyr::summarise(total_landings = sum(total_landings)) %>% 
      mutate(metacomm = "VA"),
    
    md_landings %>% 
      filter(species %in% c("Atlantic croaker",
                            "channel catfish",
                            "striped bass",
                            "gizzard shad",
                            "white perch")) %>% 
      group_by(species, year) %>% 
      dplyr::summarise(total_landings = sum(total_landings)) %>% 
      mutate(metacomm = "MD")
  )

annual_harvests <- ggplot(data = total_harvests) +
  geom_line(aes(y = total_landings, x = year, color = species),
            linewidth = 1) +
  labs(y = "Harvests (lbs)") +
  scale_color_manual(values = fill_vec) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     labels = function(x) format(x, scientific = TRUE)) +
  facet_wrap(~metacomm) +
  theme_bw() +
  theme(axis.title.y = element_markdown(size = 15),
        strip.background = element_rect(fill = 'transparent', color = NA),
        strip.text = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 14))

biomass_and_harvests <- annual_biomass +
  annual_harvests + 
    plot_layout(nrow = 2) +
    plot_annotation(tag_levels = "a",
                    tag_prefix = "(",
                    tag_suffix = ")")

ggsave(biomass_and_harvests,
       file = here::here("figs/biomass_harvest_annual_time_series.png"),
       dpi = 300,
       width = 7,
       height = 5.5)

# trends in sae and cpe-----
## community----
sae_bio_mod <- lm(sae_bio ~ year, data = ndf_va)
summary(sae_bio_mod)
acf(residuals(sae_bio_mod))
sae_bio_trend = 
  predict(sae_bio_mod,
          newdata = tibble(year = 2002:2018),
          se.fit = T)

cpe_bio_mod <- gls(cpe_bio ~ year, data = ndf_va,
                   correlation = corARMA(p = 2))
summary(cpe_bio_mod)
acf(residuals(cpe_bio_mod, type = "normalized"))
cpe_bio_trend = 
  predict(cpe_bio_mod,
          newdata = tibble(year = 2002:2018),
          se.fit = T)

## harvests----
sae_land_mod <- lm(sae_land ~ year, data = ndf_va)
summary(sae_land_mod)
acf(residuals(sae_land_mod))
sae_land_trend = 
  predict(sae_land_mod,
          newdata = tibble(year = 2002:2018),
          se.fit = T)

cpe_land_mod <- gls(cpe_land ~ year, data = ndf_va,
                    correlation = corAR1())
summary(cpe_land_mod)
acf(residuals(cpe_bio_mod, type = "normalized"))
cpe_land_trend = 
  predict(cpe_land_mod,
          newdata = tibble(year = 2002:2018),
          se.fit = T)

asynchrony_trend_preds <- 
  tibble(sae_bio_pred = sae_bio_trend$fit,
         sae_bio_se = sae_bio_trend$se.fit,
         
       cpe_bio_pred = cpe_bio_trend$fit,
       cpe_bio_se = cpe_bio_trend$se.fit,
       
       sae_land_pred = sae_land_trend$fit,
       sae_land_se = sae_land_trend$se.fit,
       
       cpe_land_pred = cpe_land_trend$fit,
       cpe_land_se = cpe_land_trend$se.fit) %>% 
  mutate(year = 2002:2018)

asynchrony_trend_stats <- 
  bind_rows(
  broom::tidy(sae_bio_mod) %>% 
    filter(term == "year") %>% 
    mutate(Term = "SAE_species",
           `Error structure` = "iid") %>% 
  dplyr::select(Term,
                Trend = estimate,
                `Std. error` = std.error,
                `T statistic` = statistic,
                `P value` = p.value,
                `Error structure`),

  broom.mixed::tidy(cpe_bio_mod) %>% 
  filter(term == "year") %>% 
  mutate(Term = "CPE_species",
         `Error structure` = "AR(2)",
         System = "Community") %>% 
  dplyr::select(Term,
                Trend = estimate,
                `Std. error` = std.error,
                `T statistic` = statistic,
                `P value` = p.value,
                `Error structure`),

  broom.mixed::tidy(cpe_land_mod) %>% 
    filter(term == "year") %>% 
    mutate(Term = "CPE_harvests",
           `Error structure` = "AR(1)",
           System = "Fishery") %>% 
    dplyr::select(Term,
                  Trend = estimate,
                  `Std. error` = std.error,
                  `T statistic` = statistic,
                  `P value` = p.value,
                  `Error structure`),
  
  broom.mixed::tidy(sae_land_mod) %>% 
    filter(term == "year") %>% 
    mutate(Term = "SAE_harvests",
           `Error structure` = "iid",
           System = "Fishery") %>% 
    dplyr::select(Term,
                  Trend = estimate,
                  `Std. error` = std.error,
                  `T statistic` = statistic,
                  `P value` = p.value,
                  `Error structure`)
)

write.csv(asynchrony_trend_stats,
          file = here::here("data/asynchrony_trend_stats.csv"),
          row.names = F)

# time series of sae and cpe-----
sae_bio_trend_plt <- ndf_va %>% 
  ggplot() +
    geom_line(data = asynchrony_trend_preds,
              aes(y = sae_bio_pred, x = year),
              linewidth = 1) +
    geom_ribbon(data = asynchrony_trend_preds,
                aes(ymin = sae_bio_pred - 1.96 * sae_bio_se,
                    ymax = sae_bio_pred + 1.96 * sae_bio_se,
                    x = year),
                alpha = 0.25) +
  geom_line(aes(y = sae_bio, x = year)) +
    labs(y = "SAE<sub>Species</sub>") +
    geom_text(data = asynchrony_trend_stats %>% 
                filter(Term == "SAE_species"),
              aes(y = 1.8, x = 2003.5, label = paste("P < 0.001")), #hacky here - couldn't get print to work
              size = 4) +
    theme_bw() +
    theme(axis.title.y = element_markdown(size = 14),
          axis.title.x = element_blank()) 

cpe_bio_trend_plt <- ndf_va %>% 
  ggplot() +
  geom_line(data = asynchrony_trend_preds,
            aes(y = cpe_bio_pred, x = year),
            linewidth = 1) +
  geom_ribbon(data = asynchrony_trend_preds,
              aes(ymin = cpe_bio_pred - 1.96 * cpe_bio_se,
                  ymax = cpe_bio_pred + 1.96 * cpe_bio_se,
                  x = year),
              alpha = 0.25) +
  geom_line(aes(y = cpe_bio, x = year)) +
  labs(y = "CPE<sub>Species</sub>") +
  geom_text(data = asynchrony_trend_stats %>% 
              filter(Term == "CPE_species"),
            aes(y = 1.7, x = 2003.5, label = paste("P =",round(`P value`,
                                                             3))),
            size = 4) +
  theme_bw() +
  theme(axis.title.y = element_markdown(size = 14),
        axis.title.x = element_blank()) 

sae_land_trend_plt <- ndf_va %>% 
  ggplot() +
  geom_line(aes(y = sae_land, x = year)) +
  labs(y = "SAE<sub>Harvest</sub>") +
  theme_bw() +
  theme(axis.title.y = element_markdown(size = 14),
        axis.title.x = element_blank()) 

cpe_land_trend_plt <- ndf_va %>% 
  ggplot() +
  geom_line(data = asynchrony_trend_preds,
            aes(y = cpe_land_pred, x = year),
            linewidth = 1) +
  geom_ribbon(data = asynchrony_trend_preds,
              aes(ymin = cpe_land_pred - 1.96 * cpe_land_se,
                  ymax = cpe_land_pred + 1.96 * cpe_land_se,
                  x = year),
              alpha = 0.25) +
  geom_line(aes(y = cpe_land, x = year)) +
  labs(y = "CPE<sub>Harvest</sub>") +
  geom_text(data = asynchrony_trend_stats %>% 
              filter(Term == "CPE_harvests"),
            aes(y = 1.7, x = 2003.5, label = paste("P =",
                                                   round(`P value`,
                                                             3))),
            size = 4) +
  theme_bw() +
  theme(axis.title.y = element_markdown(size = 14),
        axis.title.x = element_blank()) 

async_trend_plts <-
  sae_bio_trend_plt +
  cpe_bio_trend_plt +
   sae_land_trend_plt +
    cpe_land_trend_plt +
  plot_annotation(tag_levels = "a",
                  tag_suffix = ")",
                  tag_prefix = "(")

ggsave(async_trend_plts,
       file = here::here("figs/async_trend_plts.png"),
       dpi = 300,
       width = 7.5,
       height = 5)

# landings shift-----
land_shift_df <- 
va_landings %>% 
  filter(species %in% c("spot", "Atlantic croaker","striped bass")) %>% 
  right_join(.,expand_grid(month = 1:12,
                           species = c("Atlantic croaker", "spot", "striped bass"),
                           year = 2002:2018)) %>% 
  mutate(total_landings = ifelse(is.na(total_landings),
                                 0,
                                 total_landings)) %>% 
  group_by(species, year) %>% 
  dplyr::summarise(weighted_month = weighted.mean(month, total_landings)) 

ac_shift <- lm(weighted_month ~ year, data = land_shift_df %>% 
                 filter(species == "Atlantic croaker"))
summary(ac_shift)
acf(residuals(ac_shift))

sp_shift <- lm(weighted_month ~ year, data = land_shift_df %>% 
                 filter(species == "spot"))
summary(sp_shift)
acf(residuals(sp_shift))

sb_shift <- lm(weighted_month ~ year, data = land_shift_df %>% 
                 filter(species == "striped bass"))
summary(sb_shift)
acf(residuals(sb_shift))

shift_summaries <- 
  bind_rows(
    broom::tidy(ac_shift) %>% 
      mutate(species = "Atlantic croaker"),
    broom::tidy(sp_shift)%>% 
      mutate(species = "spot"),
    broom::tidy(sb_shift)%>% 
      mutate(species = "striped bass")
  ) %>% 
  filter(term == "year")
  
## figure---- 
shift_plt <- ggplot() +
  geom_smooth(data = land_shift_df %>% 
               filter(species %in% c("striped bass", "Atlantic croaker")),
            aes(y = weighted_month, x = year, group = species),
            method = "lm",
            color = "black") + 
  labs(y = "Average month weighted by monthly harvests") +
  geom_line(data = land_shift_df,
             aes(y = weighted_month, x = year, color = species),
            linewidth = 1) + 
  geom_text(data = shift_summaries %>% 
              filter(species %in% c("Atlantic croaker",
                                    "striped bass")),
            aes(y = c(7.5, 4.75),
                x = c(2007.5, 2012.5),
                label = c("P < 0.001", "P < 0.001"))) +
  scale_color_manual(values = fill_vec) +
    theme_bw() +
  theme(axis.title.y = element_markdown(size = 10),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10)) 

ggsave(shift_plt,
       file = here::here("figs/month_shift_plt.png"),
       dpi = 300,
       width = 5,
       height = 4)


## table----
harvest_shift_stats <- shift_summaries %>% 
  filter(term == "year") %>% 
  mutate(Term = "Species",
         `Error structure` = "iid") %>% 
  dplyr::select(Species = species,
                Trend = estimate,
                `Std. error` = std.error,
                `T statistic` = statistic,
                `P value` = p.value,
                `Error structure`)

write.csv(harvest_shift_stats,
          here::here("data/harvest_shift_stats.csv"),
          row.names = F)

# effort trends-----
cpue <-
  all_landings_complete %>% 
  filter((species %in% va_specs & metacomm == "VA") |
           (species %in% md_specs & metacomm == "MD")) %>% 
  group_by(species, year, metacomm) %>%
  dplyr::summarise(cpue = sum(total_landings)/sum(total_trips)) 
# ggplot() +
# geom_line(aes(y = cpue, x = year,
#               color = species)) +
# facet_wrap(~metacomm)

va_cpue <- 
  cpue %>% 
  filter(metacomm == "VA")

ac_cpue_mod <- lm(cpue ~ year, data = va_cpue %>% 
                    filter(species == "Atlantic croaker"))
summary(ac_cpue_mod)
acf(residuals(ac_cpue_mod))
s <- DHARMa::simulateResiduals(ac_cpue_mod, n = 1000);plot(s)

sp_cpue_mod <- lm(cpue ~ year, data = va_cpue %>% 
                    filter(species == "spot"))
summary(sp_cpue_mod)
acf(residuals(sp_cpue_mod))

sb_cpue_mod <- lm(cpue ~ year, data = va_cpue %>% 
                    filter(species == "striped bass"))
summary(sb_cpue_mod)
acf(residuals(sb_cpue_mod))

effort_declines <- 
  ndf_va %>% 
  dplyr::select(`Atlantic croaker` = ac_effort, spot = sp_effort,
                `striped bass` = sb_effort, year) %>% 
  gather(common, value, -year) %>% 
ggplot() + 
  geom_line(aes(y = value, x = year, 
                color = common),
            linewidth = 1) +
  scale_color_manual(values = fill_vec) +
  labs(y = "Total trips") +
  theme_bw() +
  theme(axis.title.y = element_markdown(size = 10),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10)) +
va_cpue %>% 
  ggplot() + 
  geom_line(aes(y = cpue, x = year, 
                color = species),
            linewidth = 1) +
  scale_color_manual(values = fill_vec) +
  labs(y = "Catch per trip (lbs trip<sup>-1</sup>)") +
  theme_bw() +
  theme(axis.title.y = element_markdown(size = 10),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10)) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a",
                  tag_suffix = ")",
                  tag_prefix = "(")  & 
  theme(legend.position = "bottom")

ggsave(effort_declines,
       file = here::here("figs/effort_declines_cpue.png"),
       dpi = 300,
       width = 7.5,
       height = 3.5)

# harvest stability trends-----
stab_trend_mod_va <- lm(s_land ~ year, data = ndf_va)
summary(stab_trend_mod_va)
acf(residuals(stab_trend_mod_va))

stab_trend_mod_md <- lm(s_land ~ year, data = ndf_md)
summary(stab_trend_mod_md)
acf(residuals(stab_trend_mod_md))

stab_summary_md <- broom::tidy(stab_trend_mod_md)
stab_summary_va <- broom::tidy(stab_trend_mod_va)

stab_summary_stats <-
  bind_rows(
    stab_summary_md %>% 
      filter(term == "year") %>% 
      mutate(Term = "S_Portfolio",
             Region = "Maryland",
             `Error structure` = "iid") %>% 
      dplyr::select(Region,
                    Term,
                    Trend = estimate,
                    `Std. error` = std.error,
                    `T statistic` = statistic,
                    `P value` = p.value,
                    `Error structure`),
    stab_summary_va %>% 
      filter(term == "year") %>% 
      mutate(Term = "S_Portfolio",
             Region = "Virginia",
             `Error structure` = "iid") %>% 
      dplyr::select(Region,
                    Term,
                    Trend = estimate,
                    `Std. error` = std.error,
                    `T statistic` = statistic,
                    `P value` = p.value,
                    `Error structure`)
  ) 

write.csv(stab_summary_stats,
          file = here::here("data/stab_summary_stats.csv"),
          row.names = F)

va_stab_plt <- ggplot(data = ndf_va) +
  geom_smooth(aes(y = s_land, x = year), method = "lm",
              color = "black") +
  geom_line(aes(y = s_land, x = year)) +
  geom_text(aes(y = 1.5, x = 2016, label = "P < 0.001")) +
  geom_text(aes(y = 0.5, x = 2003, label = "Virginia")) +
  labs(y = "Portfolio harvest stability (S<sub>Portfolio, L</sub>)") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10)) 

md_stab_plt <- 
  ggplot(data = ndf_md) +
  geom_line(aes(y = s_land, x = year)) +
  labs(y = "Portfolio harvest stability (S<sub>Portfolio, L</sub>)") +
  geom_text(aes(y = 1.25, x = 2003.25, label = "Maryland")) +
  theme_bw() +
  theme(axis.title.y = element_markdown(size = 10),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10)) 

stab_plts <- 
  md_stab_plt +
  va_stab_plt +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")")

ggsave(stab_plts,
       file = here::here("figs/harvest_stability_trends.png"),
       dpi = 300,
       width = 8,
       height = 3)

# within-year biomass and harvest figures----

ylabf <- function(x) {
  parse(text = paste0(x / 10^5))
}

## harvests-----
md_specs <- c("Atlantic croaker",
              "channel catfish",
              "gizzard shad",
              "white perch",
              "striped bass")

va_specs <- c("Atlantic croaker",
              "spot",
              "striped bass")

md_harvests_plt <-
  md_landings %>% 
  filter(species %in% md_specs) %>% 
  right_join(.,expand_grid(month = 1:12,
                           species = md_specs,
                           year = 2002:2018)) %>% 
  mutate(total_landings = ifelse(is.na(total_landings),
                                 0,
                                 total_landings)) %>% 
  mutate(year = factor(year, levels = 2002:2018)) %>% 
  ggplot() +
  geom_segment(data = expand_grid(month = month_vec,
                                  metacomm = "MD","VA"),
               aes(x = month,
                   xend = month,
                   y = 0,
                   yend = Inf),
               alpha = 0.25,
               linetype = 2) +
  geom_line(aes(y = total_landings, x = month, color = year,
                group = year)) +
  facet_wrap(~species) +
  labs(x = "Month",
       y = labs(y = "Harvests (&times;10<sup>5</sup> lbs)"),
       color = "Year") +
  dream::theme_fade()+
  scale_x_continuous(breaks = month_vec) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     labels = ylabf) +
  scale_color_manual(values=as.vector(coolwarm(17))) +
  theme(axis.title.y = element_markdown(size = 8)) 

va_harvests_plt <-
  va_landings %>% 
  filter(species %in% va_specs) %>% 
  right_join(.,expand_grid(month = 1:12,
                           species = va_specs,
                           year = 2002:2018)) %>% 
  mutate(total_landings = ifelse(is.na(total_landings),
                                 0,
                                 total_landings)) %>% 
  mutate(year = factor(year, levels = 2002:2018)) %>% 
  ggplot() +
  geom_segment(data = expand_grid(month = month_vec,
                                  metacomm = "MD","VA"),
               aes(x = month,
                   xend = month,
                   y = 0,
                   yend = Inf),
               alpha = 0.25,
               linetype = 2) +
  geom_line(aes(y = total_landings, x = month, color = year,
                group = year)) +
  facet_wrap(~species) +
  labs(x = "Month",
       y = labs(y = "Harvests (&times;10<sup>5</sup> lbs)"),
       color = "Year") +
  dream::theme_fade()+
  scale_x_continuous(breaks = month_vec) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     labels = ylabf) +
  scale_color_manual(values=as.vector(coolwarm(17)))  +
  theme(axis.title.y = element_markdown(size = 8))

within_yr_harvests <- 
  md_harvests_plt +
  va_harvests_plt +
  plot_layout(nrow = 2, guides = "collect") +
  plot_annotation(tag_levels = "a",
                  tag_suffix = ")",
                  tag_prefix = "(")

ggsave(within_yr_harvests,
       file = here::here("figs/within_year_harvests.png"),
       dpi = 300,
       height = 6.5,
       width = 8)  

## biomass-----
within_yr_biomass <- 
pred_out_seas %>% 
  filter(metacomm == "MD",
         common %in% md_specs) %>% 
  mutate(year = factor(year, levels = 2002:2018)) %>% 
  ggplot() +
  geom_line(aes(y = est, x = month, color = year,
                group = year)) +
  facet_wrap(~common, scales = "free_y") +
  labs(x = "Month",
       y = labs(y = "Biomass (&times;10<sup>5</sup> kg)"),
       color = "Year") +
  dream::theme_fade()+
  scale_x_continuous(breaks = month_vec,
                     limits = c(1,12)) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     labels = ylabf) +
                     # labels = function(x) format(x, scientific = TRUE)) +
  scale_color_manual(values=as.vector(coolwarm(17))) +
  theme(axis.title.y = element_markdown(size = 8)) +
  
pred_out_seas %>% 
  filter(metacomm == "VA",
         common %in% va_specs) %>% 
  mutate(year = factor(year, levels = 2002:2018)) %>% 
ggplot() +
  geom_line(aes(y = est, x = month, color = year,
                group = year)) +
  facet_wrap(~common, scales = "free_y") +
  dream::theme_fade()+
  labs(x = "Month",
       y = labs(y = "Biomass (&times;10<sup>5</sup> kg)"),
       color = "Year") +
  scale_x_continuous(breaks = month_vec,
                     limits = c(1,12)) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     labels = scientific_10_5)+
  scale_color_manual(values=as.vector(coolwarm(17))) +
  plot_layout(nrow = 2, guides = "collect",
              heights = c(1, 0.4)) +
  plot_annotation(tag_levels = "a",
                  tag_suffix = ")",
                  tag_prefix = "(") &
  theme(axis.text = element_markdown(size = 8),
        axis.title.y = element_markdown(size = 8))

ggsave(within_yr_biomass,
       file = here::here("figs/within_year_biomass.png"),
       dpi = 300,
       height = 6.5,
       width = 8)  
