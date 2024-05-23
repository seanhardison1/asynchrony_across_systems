library(tidyverse)
library(readxl)
library(tsibble)
library(lubridate)

load(here::here("data/processed_landings.rdata"))
load(here::here("data/chesmmap.rdata"))


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
  group_by(species, year) %>% 
  dplyr::summarise(total_value = sum(total_value)) %>% 
  group_by(species) %>% 
  dplyr::summarise(average_value = mean(total_value)) %>% 
  arrange(desc(average_value))

md_landings %>% 
  group_by(species, year) %>% 
  dplyr::summarise(total_value = sum(total_value)) %>% 
  group_by(species) %>% 
  dplyr::summarise(average_value = mean(total_value)) %>% 
  arrange(desc(average_value))
