library(tidyverse)
library(sf)
library(tsibble)
library(lubridate)
library(readxl)

ncrs <- "+proj=utm +zone=18 +datum=WGS84 +units=km +no_defs"


source(here::here("R/stability_calculations/part_stab_functions2.r"))

# landings-----
total_trips <-
  readxl::read_excel(here::here("data/Hardison 2002-2018_all.xlsx"),
                     sheet = 2) %>% 
  slice(1:6557) %>% 
  mutate(ym = yearmonth(as.Date(paste(MONTH,"01",sep="-"))),
         year = year(ym),
         month = month(ym)) %>%
  dplyr::select(-MONTH) %>% 
  rename_with(str_to_lower, SPEC_GROUP:month) %>% 
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
  group_by(year, month, species) %>% 
  dplyr::summarise(total_trips = sum(trips, na.rm = T))

total_trips_combined <- read_xlsx(here::here("data/2023 Hardison Request 0905 without Species.xlsx"),
                          sheet = 1) %>% 
  filter(WATER_Name != "TANGIER SOUND") %>% 
  dplyr::rename(year = YEAR,
                month = MONTH,
                trips = TRIPS) %>% 
  mutate(ymon = yearmonth(paste(year, month,"01",sep="-")),
         year = year(ymon),
         month = month(ymon)) %>% 
  group_by(year, month, ymon) %>% 
  filter(month %in% c(3,5,7,9,11)) %>% 
  dplyr::summarise(combined_trips = sum(trips))

total_trips_annual <- total_trips_combined %>% 
  group_by(year) %>% 
  dplyr::summarise(total_trips = sum(combined_trips))

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
  left_join(.,total_trips) %>%
  filter(!is.na(total_trips)) %>% 
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
  filter(!str_detect(species, "crab")) %>%
  {. ->> species_df} #%>% 
# mutate(metacomm = "VA",
# lpt = total_landings/total_trips)# %>%
# dplyr::rename(metacomm = location)


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

# md_specs <- unique(md_landings$species)
# va_specs <- unique(va_landings$species)
# specs <- unique(c(md_specs, va_specs))
# specs <- specs[!specs %in% c("blue crab","Atlantic menhaden","horseshoe crab",
#                              "horsehoe crab (male)", "menhaden","channel catfish",
#                              "blue catfish")]
# specs[str_detect(specs, "striped bass")] <- "striped bass"

selected_specs <-
  c("Atlantic croaker",
    "weakfish",
    "summer flounder",
    "spot",
    "white perch",
    "channel catfish",
    "striped bass",
    "sheepshead",
    "kingfishes",
    "gizzard shad")

# Combine landings data from MD and VA
month_vec <- c(3,5,7,9,11)
all_landings <- bind_rows(va_landings %>% 
                            mutate(metacomm = "VA"),
                          md_landings %>% 
                            mutate(metacomm = "MD")) %>% 
  mutate(total_landings = total_landings / 2.20462262) %>% 
                            # filter(year > 2004)) %>% # trip data is only valid in MD after 2004
  {. ->> all_landings_complete} %>% 
  filter(species %in% selected_specs) %>% 
  # mutate(month = ifelse(month %in% c(1,2,3),
  #                        3,# "winter-spring",
  #                        ifelse(month %in% c(4,5),
  #                              5,# 'spring',
  #                        ifelse(month %in% c(6,7),
  #                               7,#"summer",
  #                               ifelse(month %in% c(8,9),
  #                                      9,#'fall',
  #                                      ifelse(month %in% c(10,11,12),
  #                                             11)))))) %>% #'fall-winter'
  group_by(year, metacomm, species, month) %>%
  dplyr::summarise(total_value = sum(total_value),
                   total_landings = sum(total_landings),
                   total_trips = sum(total_trips))

# ggplot(all_landings %>% filter(metacomm == "MD")) +
#   geom_line(aes(y = total_landings, x = month, color = species)) +
#   facet_wrap(~year)
         # month %in% month_vec) #%>%  # Removing months where no data from ChesMMAP
# filter((species %in% c("Atlantic croaker","spot","striped bass") & metacomm == "VA") |
#          (species %in% c("striped bass","white perch","Atlantic croaker",
#                         "spot") & metacomm == "MD"))
# management area polygons for biomass indices
va_mgmt <- NULL
for (i in c("CBLE.kml", "CBLW.kml", "CBUE.kml", 
            "CBUW.kml","CB_LW_MD.kml","CB_MID_MD.kml",
            "CB_UP_MD.kml")){ 
  init <- st_read(here::here('data',i)) %>% 
    st_zm() %>% 
    st_transform(.,ncrs)
  
  assign("va_mgmt", rbind(init,
                          va_mgmt))
}

save(selected_specs, month_vec, all_landings,total_trips_combined,
     total_trips_annual,
     ncrs, chesmmap_poly, va_mgmt,all_landings_complete,
     file = here::here("data/processed_landings.rdata"))

