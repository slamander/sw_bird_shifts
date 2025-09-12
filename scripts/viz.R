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

region <- read_rds("data/region.rds")
dem <- read_rds("data/dem.rds")

dem <- get_elev_raster(
  data.frame(x = c(-130, -100), y = c(28, 51)),
  prj = st_crs(region), 
  z = 7) %>%
  crop(region) %>%
  mask(region) %>%
  as.data.frame(xy = T) 
names(dem) <- c("x", "y", "elev")

pija_ebd <- pija %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326) %>%
  st_transform(st_crs(region))

pija_imbcr <- imbcr %>%
  filter(en_common_name == "Pinyon Jay")  %>%
  st_transform(st_crs(region))

cols <- c("ebd" = "#f04546", "imbcr" = "#3591d1","bbs" = "#62c76b")

ggplot() +
  geom_raster(data = dem, aes(x = x, y = y, fill = elev)) +
  geom_sf(data = region, fill = NA) +
  geom_sf(data = pija_ebd, shape = "+", aes(col = "ebd")) +
  geom_sf(data = pija_imbcr, shape = "+", aes(col = "imbcr")) +
  scale_colour_manual("data", values = cols) +
  scale_fill_gradient(low = "grey40", high = "white", na.value = NA) +
  theme_void() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.05, 0.25))
ggsave("figures/pija_ebd_imbcr.jpeg", dpi = 600)

