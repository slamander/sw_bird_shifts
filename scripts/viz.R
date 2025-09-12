pacman::p_load(
  tidyverse
  , sf
  , terra
  , tigris
  , rnaturalearth
  , elevatr
  , ggnewscale
  , ggpubr
  #, 
  #, 
)

region <- read_rds("data/region.rds")
dem <- read_rds("data/dem.rds")

pija_ebd <- read_rds("data/ebd/pija_ebd.rds")
pija_imbcr <- read_rds("data/bcr/pija_imbcr.rds")
pija_bbs <- read_rds("data/bbs/pija_bbs_strat.rds")

# cols <- c(
#   "eBird" = "#f04546"
#   , "IMBCR" = "#3591d1"
#   , "BBS" = "#62c76b"
#   # , "bbs" = "#62c76b"
# )

base_map <- ggplot() +
  geom_raster(
    data = dem_df,
    aes(x, y, fill = elev),
    guide = "none") +
  scale_fill_hypso_tint_c() +
  geom_sf(data = region, fill = NA) +
  theme_void() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.05, 0.25))

gg_data_3panel <- ggarrange(
  # eBird data
  base_map + 
    geom_sf(data = pija_ebd, shape = "+", col = "#3590d1ff"),
  # IMBCR data
  base_map + 
    geom_sf(data = pija_imbcr, shape = "+", col = "#3590d1ff"),
  # BBS data
  base_map + 
    geom_sf(data = pija_bbs, aes(alpha = count), col = NA, fill = "#3590d1ff") +
    scale_alpha_continuous(
      "counts",
      range = c(0.25, 1), 
      na.value = NA), 
  ncol = 3, 
  nrow = 1
)

ggsave(
  "figures/pija_ebd_imbcr_bbs_3panel.jpeg", 
  height = 4, 
  width = 12, 
  dpi = 600)
