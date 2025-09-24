pacman::p_load(
  tidyverse
  , sf
  , terra
  , tidyterra
  , elevatr
  , ggnewscale
  , ggpubr
  #, 
  #, 
)

crs <- read_rds("data/region/crs.rds")
region <- read_rds("data/region/poly.rds")
dem <- read_rds("data/region/dem.rds")

pija_abund <- rast("data/sat/pija/abundance_2023.tif") %>%
  project(crs(region)) %>%
  crop(region, mask = T)
pija_abund_filt <- pija_abund; pija_abund_filt[pija_abund_filt == 0] <- NA

pija_ebd <- read_rds("data/ebd/pija_ebd_strat.rds")
pija_imbcr <- read_rds("data/bcr/pija_imbcr.rds")
pija_bbs <- read_rds("data/bbs/pija_bbs_strat.rds")

base_map <- ggplot() + 
  geom_raster(
    data = dem,
    aes(x, y, fill = elev),
    show.legend = FALSE) +
  scale_fill_hypso_tint_c(
    palette = "etopo1",
    labels = scales::label_number(),
    breaks = c(-5000, 0, 2000, 4500)) +
  geom_sf(data = region, fill = NA) 

gg_data_4panel <- ggarrange(
  # eBird status and trends
 base_map + 
  new_scale_fill() +
  geom_raster(
    data = as.data.frame(pija_abund_filt, xy = T),
    aes(x = x, 
        y = y,
        alpha = resident,
        fill = resident)) +
  scale_fill_gradient(
    "Abundance",
    low = "#63beffff", 
    high = "#006cb9ff", 
    na.value = NA) +
  scale_alpha_continuous(
    range = c(0.5, 1)) +
    guides(alpha = "none") +
  theme_void() +
  theme(legend.position = "bottom"),

  # IMBCR map
  base_map + 
    geom_sf(data = pija_imbcr, 
      shape = "+", size = 5, alpha = 0.75,
      aes(col = en_common_name)) +
    scale_color_manual(
      "Presences", 
      values = "#3590d1ff") +
  theme(legend.position = "bottom"), 

  # eBird map
  base_map + 
    geom_sf(data = pija_ebd, 
      aes(alpha = count), 
      col = NA, 
      fill = "#3590d1ff") +
    scale_alpha_continuous(
      "Observations",
      limits = c(1, 92462),
      range = c(0.5, 1), 
      na.value = NA) +
  theme(legend.position = "bottom"), 

  # BBS map
  base_map + 
    geom_sf(data = pija_bbs, 
      aes(alpha = count), 
      col = NA, 
      fill = "#3590d1ff") +
    scale_alpha_continuous(
      "Observations",
      range = c(0.5, 1), 
      limits = c(24, 2043),
      na.value = NA) +
  theme(legend.position = "bottom"),

  # panel params
  ncol = 4, 
  nrow = 1,
  align = "h"
)

ggsave(
  "figures/pija_ebd_imbcr_bbs_4panel.jpeg", 
  height = 7, 
  width = 16, 
  dpi = 600)
