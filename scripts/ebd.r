pacman::p_load(
  tidyverse
  , here
  , auk
  # , 
  # , 
)

obs_file <- here("data/ebd/ebd/observations.txt")
chk_file <- here("data/ebd/ebd/checklists.txt")

pija_ebd_raw <- obs_file %>%
  auk_ebd() %>%
  auk_species(species = "Pinyon Jay") %>%
  auk_country(country = c("US")) %>%
  auk_bbox(bbox = c(-130, 25, -100, 50)) %>%
  # auk_date(date = c("2012-01-01", "2012-12-31")) %>%
  # auk_time(start_time = c("06:00", "09:00")) %>%
  # auk_duration(duration = c(0, 60)) %>%
  auk_complete() %>% 
  auk_filter(file = "data/ebd/filters/ebd_pija.txt", overwrite = TRUE) %>%
  read_ebd()

pija_ebd <- pija_ebd_raw %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(4326)) %>%
  st_transform(st_crs(region))

write_rds(pija_ebd, "data/ebd/pija_ebd.rds")
