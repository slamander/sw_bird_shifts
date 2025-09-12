pacman::p_load(
  tidyverse
  , sf
  , terra
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
  st_crop(xmin = -130, ymin = 28, xmax = -100, ymax = 50); plot(region["admin"])
write_rds(region, "data/region.rds")

dem <- get_elev_raster(
  data.frame(x = c(-130, -100), y = c(28, 51)),
  prj = st_crs(region), 
  z = 7) %>%
  crop(region) %>%
  mask(region) %>%
  as.data.frame(xy = T) 
names(dem) <- c("x", "y", "elev")
write_rds(dem, "data/dem.rds")
