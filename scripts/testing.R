pacman::p_load(
  tidyverse
  , sf
  , terra
  , tidyterra
  , ggnewscale
  , ggpubr
  #, 
  #, 
)

species <- read_rds("data/species/prop_prelim.rds")

region <- read_rds("data/region/poly.rds")
strat <- read_rds("data/region/strat.rds")

# Fucked up with the year
# bbs_sp_strat <- read_rds("data/bbs/species/strat/alhu.rds")
# bbs_sp_point <- read_rds("data/bbs/species/point/alhu.rds")

bbs_strat <- read_rds("data/bbs/bbs_strat.rds") # preload it instead??

# gather bbs point data
bbs_sp <- bbsBayes::prepare_data(
  bbs_strat,
  species_to_run = "Allen's Hummingbird",
  model = "slope")

# clean it up
bbs_sp_df <- data.frame(
  count = bbs_sp$count,
  strat_name = bbs_sp$strat_name, 
  strat = bbs_sp$strat_name,
  obser = bbs_sp$obser,
  firstyr = bbs_sp$firstyr, 
  route = bbs_sp$route, 
  year = bbs_sp$r_year, 
  month = bbs_sp$month, 
  day = bbs_sp$day
)

bbs_sp_strat <- bbs_sp_df %>%
  group_by(strat_name, year) %>%
  summarize(count = sum(count, na.rm = T)) %>%
  left_join(strat) %>%
  st_set_geometry("geometry")

bbs_sp_strat_sum <- bbs_sp_strat %>%
  group_by(strat_name) %>%
  summarize(count = sum(count, na.rm = T))

ggplot() + 
  geom_sf(data = region, fill = "grey80") +
  geom_sf(data = bbs_sp_strat_sum, aes(fill = count)) +
  # facet_wrap(~year) +
  theme_void()

bbs_strat <- lapply(
  species$alpha, 
    FUN = function(x){ 
      read_rds(paste0("data/bbs/species/strat/", x, ".rds")) %>%
        mutate(alpha = x)
    }
  ) %>%
  bind_rows() %>%
  left_join(species)

bbs_strat_sum <- bbs_strat %>%
  group_by(strat_name, en_common_name) %>%
  summarize(count = sum(count, na.rm = T))

ggplot() + 
  geom_sf(data = region, fill = "grey80") +
  geom_sf(data = bbs_strat_sum, aes(fill = log(count + 1))) +
  scale_fill_gradientn(
    "log(counts)",
    colours = wesanderson::wes_palette("Zissou1", 100, "continuous")) + 
  facet_wrap(~en_common_name, nrow = 3, ncol = 12) +
  theme_void() +
  theme(legend.position = "bottom")
ggsave("figures/prelim_sp_strat.jpeg", height = 6, width = 12, dpi = 600)

