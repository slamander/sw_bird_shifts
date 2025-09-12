pacman::p_load(
  tidyverse
  , here
  , sf
  , bbsBayes
  # , 
  # , 
)

fetch_bbs_data()

bbs_strat <- bbsBayes::stratify(
  by = "latlong")
bbs_strat <- pija_strat

pija_bbs_raw <- bbsBayes::prepare_data(
  bbs_strat,
  species_to_run = "Pinyon Jay",
  model = "slope")

pija_route_bbs <- bbs_strat$bird_strat %>% 
  left_join(bbs_strat$route_strat) %>%
  filter(!is.na(Longitude)) %>%
  filter(!is.na(Latitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(4326))

ggplot(pija_route_bbs) +
  geom_sf()

strat <- bbsBayes::load_map("latlong") %>%
  st_transform(st_crs(region)) %>%
  rename(strat_name = ST_12)

pija_bbs <- pija_route_bbs %>%
  left_join(strat)

plot(strat["AREA_SQ_KM"])