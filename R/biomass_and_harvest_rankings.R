library(tidyverse)
library(readxl)
library(tsibble)
library(lubridate)

load(here::here("data/processed_landings.rdata"))
load(here::here("data/chesmmap.rdata"))

# consumer price index
cpi <-
  read.table(here::here("data/Bureau of Labor Statistics Data.txt"),
           header = T,
           sep = "\t") %>% 
  dplyr::select(1:4) %>% 
  mutate(month = as.numeric(str_extract_all(Period, "\\d{2}"))) %>% 
  dplyr::rename(year = Year,
                value = Value)

ref <- cpi %>% 
  filter(year == 2010,
         month == 1) %>% 
  pull(value)

bio %>% 
  filter(common %in% c("gizzard shad",
                       "channel catfish",
                       "white perch")) %>% 
  group_by(common) %>% 
  dplyr::summarise(tot = sum(biomass > 0)/n()) %>% 
  arrange(desc(tot)) %>% 
  head(20)

md_landings <- read_excel(here::here("data/ListofSelectedSpecies.2002-2018.xlsx")) %>% 
  dplyr::rename(common = SPECNAME,
                year = Year,
                month = Month,
                pounds = POUNDS,
                dollars = Dollars,
                location = NOAACODE,
                days_out = `COUNT OF DAYS`) %>% 
  mutate(ymon = paste(year, month,"01", sep = "-"),
         ymon = yearmonth(as.Date(ymon, "%Y-%m-%d")),
         common = case_when(common == 'BUTTERFISH UNC' ~ "butterfish",
                            common == "BLACK SEA BASS" ~ "black sea bass",
                            common == "BLUE CRAB" ~ "blue crab",
                            common == "CROAKER" ~ "Atlantic croaker",
                            common == "DRUM BLACK" ~ "black drum",
                            common == "FLOUNDER SUMMER" ~ "summer flounder",
                            common == 'CATFISH - BLUE' ~ "blue catfish",
                            common == 'CATFISH - CHANNEL' ~ "channel catfish",
                            common == 'GIZZARD SHAD' ~ "gizzard shad",
                            common == 'HORSESHOE CRAB' ~ "horsehoe crab",
                            common == 'HORSESHOE CRAB (MALE)' ~ "horsehoe crab (male)",
                            common == 'SEA TROUT GRAY' ~ "weakfish",
                            common == "SB - UNCLASSIFIED" ~ "striped bass",
                            common == 'MENHADEN' ~ "Atlantic menhaden",
                            common == 'PORGY UNC' ~ "porgy",
                            common == 'RIBBON' ~ "Atlantic cutlassfish",
                            common == "SHARK - DOGFISH - SMOOTH" ~ "smooth dogfish",
                            common == "WHITING UNC" ~ "kingfishes",
                            common == "SPOT" ~ "spot",
                            common == "WHITE PERCH" ~ "white perch",
                            TRUE ~ "NA"),
         month = month(ymon),
         year = as.numeric(year)) %>%
  # filter(month %in% c(3,5,7,9,11)) %>%
  group_by(year, month, species = common) %>%
  dplyr::summarise(total_value = sum(dollars, na.rm  = T),
                   total_landings = sum(pounds, na.rm  = T),
                   total_trips = sum(days_out, na.rm = T))


va_landings <-
  readxl::read_excel(here::here("data/Hardison 2002-2018_all.xlsx"),
                     sheet = 1) %>% 
  dplyr::rename(year = 1,
                month = 2,
                spec_grp = 3,
                species = 4,
                location = 5,
                location_full = 6,
                usd = 7,
                lbs = 8) %>% 
  filter(!is.na(month)) %>%
  mutate(ymon = yearmonth(paste(year, month,"01",sep="-")),
         year = year(ymon),
         month = month(ymon)) %>% 
  group_by(year, month, species) %>% 
  filter(location != "TS") %>%
  dplyr::summarise(total_landings = sum(lbs, na.rm = T),
                   total_value = sum(usd, na.rm = T)) %>% 
  mutate(species = ifelse(species == "BLUE CRAB",
                          "CRAB, BLUE",
                          ifelse(species == "WHITING,KING",
                                 "WHITING, KING",
                                 ifelse(species == "CATFISH, NK (UNCLASSIFIED)",
                                        "CATFISH, NK",
                                        ifelse(species == "BLUE CTAFISH",
                                               "CATFISH, BLUE", 
                                               ifelse(species == "BASS, STRIPED",
                                                      "striped bass", species)))))) %>% 
  mutate(species = str_to_lower(species),
         species = case_when(species == "croaker, atlantic" ~ "Atlantic croaker",
                             species == "crab, blue" ~ "blue crab",
                             species == "crab, horseshoe" ~ "horseshoe crab",
                             # species == "bass, striped" ~ "striped bass",
                             species == "dogfish, smooth" ~ "smooth dogfish",
                             species == "drum, black" ~ "black drum",
                             species == "flounder, summer" ~ "summer flounder",
                             species == "perch, white" ~ "white perch",
                             species == "seatrout, grey" ~ "weakfish",
                             species == "shad, gizzard" ~ "gizzard shad",
                             species == "whiting, king" ~ "kingfishes",
                             species == "ray, cownose" ~ "cownose ray",
                             species == "catfish, blue" ~ "blue catfish",
                             species == "catfish, nk" ~ "channel catfish",
                             TRUE ~ species)) %>% 
  {. ->> species_df} #%>% 

md_landings %>% 
  group_by(species) %>% 
  dplyr::summarise(tot = sum(total_landings)) %>% 
  arrange(desc(tot)) %>% 
  filter(!species %in% c("blue crab",
                         "Atlantic menhaden")) %>% 
  mutate(rank = 1:nrow(.),
         harv_cum = cumsum(tot)/sum(tot)) 


va_landings %>% 
  group_by(species) %>% 
  dplyr::summarise(tot = sum(total_landings)) %>% 
  arrange(desc(tot)) %>% 
  filter(!species %in% c("blue crab",
                         "menhaden")) %>% 
  mutate(rank = 1:nrow(.),
         harv_cum = cumsum(tot)/sum(tot)) 

# average annual value
va_landings %>% 
  left_join(., cpi) %>% 
  dplyr::select(-Series.Id, -Period) %>% 
  rowwise() %>% 
  mutate(real_total_value = (total_value * ref)/value) %>% 
  ungroup() %>% 
  group_by(species, year) %>% 
  dplyr::summarise(real_total_value = sum(real_total_value)) %>% 
  group_by(species) %>% 
  dplyr::summarise(average_value = mean(real_total_value)) %>% 
  arrange(desc(average_value))

md_landings %>% 
  left_join(., cpi) %>% 
  dplyr::select(-Series.Id, -Period) %>% 
  rowwise() %>% 
  mutate(real_total_value = (total_value * ref)/value) %>% 
  ungroup() %>% 
  group_by(species, year) %>% 
  dplyr::summarise(real_total_value = sum(real_total_value)) %>% 
  group_by(species) %>% 
  dplyr::summarise(average_value = mean(real_total_value)) %>% 
  arrange(desc(average_value))


# average annual value


md_landings %>% 
  group_by(species, year) %>% 
  dplyr::summarise(total_value = sum(total_value, na.rm = T)) %>% 
  group_by(species) %>% 
  dplyr::summarise(average_value = mean(total_value)) %>% 
  arrange(desc(average_value))

# average annual harvests
va_landings %>% 
  group_by(species, year) %>% 
  dplyr::summarise(total_landings = sum(total_landings)) %>% 
  group_by(species) %>% 
  dplyr::summarise(average_landings = mean(total_landings)) %>% 
  arrange(desc(average_landings))

md_landings %>% 
  group_by(species, year) %>% 
  dplyr::summarise(total_landings = sum(total_landings)) %>% 
  group_by(species) %>% 
  dplyr::summarise(average_landings = mean(total_landings)) %>% 
  arrange(desc(average_landings))

# total value range by region
va_specs <- c("Atlantic croaker", "spot", "striped bass")

va_val <- va_landings %>% 
  filter(species %in% va_specs) %>% 
  group_by(year, species) %>% 
  dplyr::summarise(total_value = sum(total_value)) %>% 
  group_by(year) %>% 
  dplyr::summarise(total = sum(total_value)) %>% 
  mutate(state = "va")


md_specs <- c("Atlantic croaker","striped bass",
              "white perch", "gizzard shad", "channel catfish")

md_landings %>% 
  filter(species %in% md_specs) %>% 
  group_by(year, species) %>% 
  dplyr::summarise(total_value = sum(total_value)) %>% 
  group_by(year) %>% 
  dplyr::summarise(total = sum(total_value))  %>% 
  mutate(state = "md") %>% 
  bind_rows(., 
            va_val) %>% 
  group_by(year) %>% 
  dplyr::summarise(total = sum(total)) %>% 
  pull(total) %>% 
  range()
