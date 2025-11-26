pacman::p_load(
  tidyverse
  , sf
  , terra
  , tidyterra
  , ggnewscale
  , ggpubr
  , ggspatial
  #, 
  #, 
)

region <- read_rds("data/region/poly.rds")
strat <- read_rds("data/region/strat.rds")
dem <- read_rds("data/region/dem_6.rds")

species <- read_rds("data/species/prop_prelim.rds")

base_map <- ggplot() + 
  geom_raster(
    data = dem,
    aes(x, y, fill = elev),
    show.legend = FALSE) +
  scale_fill_hypso_tint_c(
    # palette = "etopo1", 
    palette = "wiki-2.0", # "etopo1"
    labels = scales::label_number(),
    breaks = c(-5000, 0, 2000, 4500)) +
  geom_sf(data = region, fill = NA)  +
  ggnewscale::new_scale_fill() +
  annotation_scale(
    bar_cols = c("grey30", "white"), 
    line_width = 0.1, 
    pad_y = unit(1.2, "cm"), 
    pad_x = unit(0.75, "cm")) +
  annotation_north_arrow(
    location = "bl", 
    which_north = "true", 
    style = north_arrow_orienteering(fill = c("grey30", "white"), 
    line_col = "grey20"),
    pad_y = unit(1.8, "cm"), 
    pad_x = unit(3.4, "cm"),
    height = unit(1, "cm"), 
    width = unit(0.75, "cm")) +
  annotate(
    'rect',
    xmin = -2350000,
    xmax = -276637.4,
    ymin = -1700000.0,
    ymax = 1650000.0,
    alpha = 0.5, # This was put back to 0.5
    fill = NA,
    col = 'black') +
  theme_void() +
  theme(
    legend.box.margin = margin(t = 7, r = 5, b = 7, l = 5),
    legend.box.background = element_rect(
      color = "black", fill = "#ffffff87", linewidth = 0.1),
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.2),
    # legend.direction = "horizontal",
    legend.title.position = "left",
    legend.title = element_text(hjust = 0.5, angle = 90))

######## BBS #########

bbs_strat <- read_rds("data/bbs/bbs_strat.rds")

bbs_strat_rich <- bbs_strat %>%
  group_by(strat_name, en_common_name) %>%
  summarize(mean_count = mean(count, na.rm = T)) %>%
  group_by(strat_name) %>%
  summarize(rich = sum(mean_count > 0, na.rm = T))

gg_bbs_rich <- base_map +
  geom_sf(data = bbs_strat_rich, aes(fill = rich), alpha = 0.6) +
  scale_fill_gradientn(
    "SGCN richness",
    colours = wesanderson::wes_palette("Zissou1", 100, "continuous")) + 
  guides(fill = guide_colourbar(barwidth = 1, barheight = 8)) +
      annotate(
    'text', 
    x = -900000,
    y = 1570000.0,
    size = 5,
    label = "Breeding Bird Survey data:",
    fontface = "bold"
  ) +
  annotate(
    'text', 
    x = -800000,
    y = 1480000.0,
    size = 3.5,
    label = "N = 36 species of greatest concern"
  ); # ggsave("figures/bbs_strat_richness.jpeg", height = 8, width = 5, dpi = 600)

###### EBD #######

ebd_strat <- read_rds("data/ebd/ebd_strat.rds")

ebd_strat_rich <- ebd_strat %>%
  group_by(strat_name, en_common_name) %>%
  summarize(mean_count = mean(count, na.rm = T)) %>%
  group_by(strat_name) %>%
  summarize(rich = sum(mean_count > 0, na.rm = T))

gg_ebd_rich <- base_map +
  geom_sf(data = ebd_strat_rich, aes(fill = rich), alpha = 0.6) +
  scale_fill_gradientn(
    "SGCN richness",
    colours = wesanderson::wes_palette("Zissou1", 100, "continuous")) + 
  guides(fill = guide_colourbar(barwidth = 1, barheight = 8)) +
      annotate(
    'text', 
    x = -550000, #-800000
    y = 1570000.0,
    size = 5,
    label = "eBird data:",
    fontface = "bold"
  ) +
  annotate(
    'text', 
    x = -800000, #-800000
    y = 1480000.0,
    size = 3.5,
    label = "N = 36 species of greatest concern"
  ); # ggsave("figures/bbs_strat_richness.jpeg", height = 8, width = 5, dpi = 600)

ggarrange(
  gg_bbs_rich, 
  gg_ebd_rich, 
  ncol = 2
); # ggsave("figures/2panel_strat_richness.jpeg", height = 8, width = 10, dpi = 600)
