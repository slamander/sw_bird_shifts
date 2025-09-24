pacman::p_load(
  tidyverse
  , janitor
  , sf
  , terra
  , rnaturalearth
  , elevatr
  #, 
  #, 
)

crs <- read_rds("data/region/crs.rds")

region <- bind_rows(
    ne_states("United States of America"),
    ne_states("Mexico"),
    ne_states("Canada")) %>%
  clean_names() %>%
  filter(!name %in% c("Alaska", "Hawaii")) %>%
  st_transform(crs) %>%
  # st_crop(xmin = -130, ymin = 28, xmax = -100, ymax = 50) %>%
  st_crop(
    xmin = -2236887.4, 
    ymin = -1700000, 
    xmax = -276637.4, 
    ymax = 1650000); plot(region["admin"])
write_rds(region, "data/region/poly.rds")

dem <- get_elev_raster(
    data.frame(
        x = c(-2336887.4, -276637.4), 
        y = c(-1700000, 1650000)),
    prj = crs, 
    z = 7) %>%
  rast() %>%
  crop(ext(-2350000, -276637.4, -1700000, 1650000))

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
write_rds(dem_df, "data/region/dem.rds")

sf_use_s2(FALSE)

strat <- bbsBayes::load_map("latlong") %>%
  rename(strat_name = ST_12) %>%
  st_intersection(
    bind_rows(
    ne_states("United States of America"),
    ne_states("Mexico"),
    ne_states("Canada")) %>%
      clean_names() %>%
      filter(!name %in% c("Alaska", "Hawaii")) %>%
      st_union() %>%
  st_transform(crs)
    ) %>%
  st_crop(region)
write_rds(strat, "data/region/strat.rds")
