pacman::p_load(
  tidyverse
  , here
  , sf
  , bbsBayes
  # , 
  # , 
)

fetch_bbs_data()

strat <- bbsBayes::load_map("latlong") %>%
  st_transform(st_crs(region)) %>%
  rename(strat_name = ST_12)

# sf_use_s2(TRUE)

strat_region <- strat %>%
  st_intersection(st_union(region))

bbs_strat <- bbsBayes::stratify(
  by = "latlong")

pija_bbs_raw <- bbsBayes::prepare_data(
  bbs_strat,
  species_to_run = "Pinyon Jay",
  model = "gam")

pija_bbs_raw <- data.frame(
  count = pija_bbs_raw$count,
  strat_name = pija_bbs_raw$strat_name, 
  strat = pija_bbs_raw$strat_name,
  obser = pija_bbs_raw$obser, 
  firstyr = pija_bbs_raw$firstyr, 
  route = pija_bbs_raw$route, 
  year = pija_bbs_raw$year, 
  month = pija_bbs_raw$month, 
  day = pija_bbs_raw$day
)

pija_bbs_strat <- pija_bbs_raw %>% 
  group_by(strat_name) %>%
  summarize(count = sum(count, na.rm = T)) %>%
  left_join(strat_region) %>%
  st_set_geometry("geometry")

pija_strat_bbs <- strat_region %>%
  left_join(
    pija_bbs_raw %>% 
      group_by(strat_name) %>%
      summarize(count = sum(count, na.rm = T))
  )

write_rds(pija_bbs_strat, "data/bbs/pija_bbs_strat.rds")
write_rds(pija_strat_bbs, "data/bbs/pija_strat_bbs.rds")

