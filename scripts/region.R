pacman::p_load(
  tidyverse
  , sf
  , terra
  , tidyterra
  , tigris
  , rnaturalearth
  , elevatr
  #, 
  #, 
)

region <- bind_rows(
    ne_states("United States of America"),
    ne_states("Mexico")) %>%
  clean_names() %>%
  st_crop(xmin = -130, ymin = 28, xmax = -100, ymax = 50); 
plot(region["admin"])
write_rds(region, "data/region.rds")

dem <- get_elev_raster(
  data.frame(x = c(-130, -100), y = c(28, 51)),
  prj = st_crs(region), 
  z = 7) %>%
  rast() %>%
  crop(region, mask = T) 

dem_hs <- shade(
  terrain(dem, "slope", unit = "radians"), 
  terrain(dem, "aspect", unit = "radians"),
  angle = 45,
  direction = 300,
  normalize = TRUE
)

plot(dem_hs, col = grey(1:100 / 100))

dem_df <- c(
  dem, 
  dem_hs) %>%
  as.data.frame(xy = T) 
names(dem_df) <- c("x", "y", "elev", "hillshade")
