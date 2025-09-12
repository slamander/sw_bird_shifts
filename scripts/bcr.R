pacman::p_load(
  tidyverse
  , sf
  , here
  , janitor
  # ,
  # ,
)

imbcr <- read_csv("data/bcr/sw_birds.csv") %>%
  clean_names() %>%
  rename(
    transect = transect_num, 
    observer = transect_visit_observer, 
    start_time = point_visit_start_time, 
    lat = point_latitude, 
    lon = point_longitude, 
    alpha = bird_code,
    en_common_name = species) %>%
  mutate(date = as.Date(date)) %>% 
  st_as_sf(coords = c("lon", "lat"),
           crs = "+proj=longlat +datum=WGS84") %>%
  bind_cols(
    st_coordinates(.)) %>%
  st_transform(st_crs(region)) %>%
  rename(lat = "Y", lon = "X") %>%
  left_join(bird_key) %>%
  filter(!is.na(alpha)) %>%
  select(transect, year, county, state, lat, lon, date, 
         aou, alpha, en_common_name, 
         binomial, order, family, genus, species) 

pija_imbcr <- imbcr %>%
  filter(en_common_name == "Pinyon Jay")

write_rds(pija_imbcr, "data/bcr/pija_imbcr.rds")
