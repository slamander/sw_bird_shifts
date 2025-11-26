pacman::p_load(
  tidyverse
  , janitor
  , sf
  , terra
  , rnaturalearth
  , elevatr
  , rmapshaper
  #, 
  #, 
)

crs <- read_rds("data/region/crs.rds")

na <- bind_rows(
    ne_states("United States of America"),
    ne_states("Mexico"),
    ne_states("Canada")) %>%
  clean_names() %>%
  filter(!name %in% c("Hawaii")) %>%
  st_transform(crs = crs) %>%
  rename(
    country = admin, 
    state = name ) %>%
  select(country, state, fips) %>%
  st_union()
st_crs(na) == crs
write_rds(na, "data/region/na.rds")

region <- na %>%
  st_crop(
    xmin = -2350000, 
    ymin = 330000, 
    xmax = -240000,
    ymax = 3500000) %>%
  ms_filter_islands(min_area = 1.5e9); plot(region)
write_rds(region, "data/region/region.rds")

strat <- region %>%
  st_make_grid(
    cellsize = c(20000, 20000),
    square = T) %>%
  st_intersection(region) %>%
  st_as_sf() %>%
  mutate(cell_id = row_number()); nrow(strat); plot(strat["cell_id"])
write_rds(strat, "data/region/strat.rds")

routes <- read_csv("data/bbs/routes/routes.csv") %>%
  clean_names() %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(4326)) %>%
  st_transform(crs = crs) %>%
  st_intersection(na) %>%
  mutate(route_id = paste(state_num, route, sep = "-")) %>%
  bind_cols(st_coordinates(.)) %>%
  select(country, country_num, state, state_num, fips, 
  route_name, route, route_id, active, stratum, bcr, X, Y, geometry) 
write_rds(routes, "data/bbs/routes/routes.rds")

dem <- get_elev_raster(
    data.frame(
        x = c(-2350000, -240000), 
        y = c(330000, 3600000)),
    prj = crs, 
    # z = 7) %>%
    z = 6) %>%
  rast() %>%
  crop(ext(-2420000, -240000, 330000, 3500000))

dem_hs <- shade(
  terrain(dem, "slope", unit = "radians"), 
  terrain(dem, "aspect", unit = "radians"),
  angle = 45,
  direction = 300,
  normalize = TRUE
)

dem_df <- c(
  dem, 
  dem_hs) %>%
  as.data.frame(xy = T) 
names(dem_df) <- c("x", "y", "elev", "hillshade")
write_rds(dem_df, "data/region/dem_6.rds")
