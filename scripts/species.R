pacman::p_load(
  tidyverse
  , sf
  , janitor
  , birdnames
  , wildlifeR
  , rebird
)

########################################################################################
# 2) load species taxonomic information to create a key
########################################################################################

ebird_key <- rebird::ebirdtaxonomy() %>%
  rename(
    binomial = sciName,
    en_common_name = comName, 
    alpha = comNameCodes, 
    family = familySciName) %>%
  mutate(
    genus = as.character(noquote(str_split_fixed(binomial, " ", 5)[,1])),
    species = as.character(noquote(str_split_fixed(binomial, " ", 5)[,2])),
    binomial = paste(genus, species, sep = " "),
    en_common_name = str_remove(en_common_name, "\\(.*\\)")) %>%
  separate_longer_delim(alpha, ",") %>%
  filter(is.na(extinctYear)) %>%
  filter(str_count(alpha) == 4) %>%
  filter(!is.na(alpha)) %>%
  distinct(alpha, binomial, .keep_all = TRUE) %>%
  select(alpha, en_common_name, binomial, order, family, genus, species)

aou_key <- wildlifeR::AOU_species_codes %>%
  as_tibble() %>%
  rename(
    aou = spp.num, 
    alpha = alpha.code, 
    en_common_name = name, 
    binomial = spp) %>%
  mutate(
    aou = ifelse(str_count(aou) == 5, aou,
                 sapply(aou, FUN = function(x){
                   paste0(rep("0", as.integer(5 - str_count(x))), x)[1]})),
    genus = as.character(noquote(str_split_fixed(binomial, " ", 2)[,1])), 
    species = as.character(noquote(str_split_fixed(binomial, " ", 2)[,2])), 
    hybrid = ifelse(binomial == "hybrid", T, F),
    binomial = ifelse(hybrid == TRUE, NA, binomial), 
    genus = ifelse(hybrid == TRUE, NA, genus),
    species = ifelse(hybrid == TRUE, NA, species)
  ) %>%
  select(aou, alpha, en_common_name, hybrid, binomial, genus, species) 

alpha_key <- birdnames::bird_list %>%
  as_tibble() %>%
  rename(alpha = alpha.code,
         en_common_name = common.name, 
         binomial = species, 
         species = specific.name) %>%
  mutate(
    binomial = ifelse(is.na(species), paste(genus, "sp.", sep = " "), binomial)) %>%
  select(alpha, en_common_name, binomial, order, family, genus, species) 

bird_key <- full_join(alpha_key, aou_key) %>%
  full_join(ebird_key)

sum(is.na(bird_key$aou))/nrow(bird_key)
sum(is.na(bird_key$alpha))/nrow(bird_key)

write_rds(bird_key, "data/species/bird_key.rds")

########################################################################################





########################################################################################
# 3) load proposal species list
########################################################################################

prop_prelim <- read_csv("data/species/proposal_primary_prelim_sp.csv") %>%
  left_join(bird_key) %>%
  mutate(
    proposal_focal = T,
    alpha = tolower(alpha)) %>%
  select(order, genus, species, 
    binomial, en_common_name,
    alpha, intersect_area_pcnt, proposal_focal) %>%
  distinct(en_common_name, .keep_all = T) %>% 
  arrange(en_common_name)
write_rds(prop_prelim, "data/species/prop_prelim.rds")

prop_second <- read_csv("data/species/proposal_secondary_prelim_sp.csv") %>%
  left_join(bird_key) %>%
  mutate(
    proposal_second = T,
    alpha = tolower(alpha)) %>%
  select(order, genus, species, 
    binomial, en_common_name,
    alpha, intersect_area_pcnt, 
    proposal_second) %>%
  distinct(en_common_name, .keep_all = T) %>% 
  arrange(en_common_name)
write_rds(prop_second, "data/species/prop_second.rds")

prop_species <- left_join(prop_second,  prop_prelim) %>%
  select(order, genus, species, 
    binomial, en_common_name,
    alpha, intersect_area_pcnt, 
    proposal_focal, proposal_second) %>% 
  arrange(en_common_name)
write_rds(prop_species, "data/species/prop_species.rds")

########################################################################################



########################################################################################
# 4) load individual database species lists
########################################################################################

imbcr_species <- read_csv("data/IMBCR/sw_birds.csv") %>%
  clean_names() %>%
  rename(
    transect = transect_num, 
    observer = transect_visit_observer, 
    start_time = point_visit_start_time, 
    lat = point_latitude, 
    lon = point_longitude, 
    alpha = bird_code,
    en_common_name = species) %>%
  mutate(date = as.Date(date)) %>% 
  st_as_sf(coords = c("lon", "lat"),
           crs = "+proj=longlat +datum=WGS84") %>%
  bind_cols(
    st_coordinates(.)) %>%
  rename(lat = "Y", lon = "X") %>%
  left_join(bird_key) %>%
  filter(!is.na(alpha)) %>%
  select(transect, year, county, state, lat, lon, date, 
         aou, alpha, en_common_name, 
         binomial, order, family, genus, species) %>%
  st_drop_geometry() %>%
  left_join(bird_key) %>%
  group_by(aou, alpha, en_common_name, 
         binomial, order, family, genus, species) %>%
  summarize(IMBCR = n()) %>%
  filter(!is.na(species))

sum(is.na(imbcr_species$aou))/nrow(imbcr_species)
sum(is.na(imbcr_species$alpha))/nrow(imbcr_species)

bbs_species <- read_csv("data/species/species_codes.csv") %>%
  mutate(binomial = paste(genus, species, sep = " ")) %>%
  left_join(bird_key) %>%
  select(aou, alpha, en_common_name, 
         binomial, order, family, genus, species) %>%
  arrange(order, family, genus, species); 

sum(is.na(bbs_species$aou))/nrow(bbs_species)
sum(is.na(bbs_species$alpha))/nrow(bbs_species)

# ebird <- read_table("data/ebd/xxx")

########################################################################################



########################################################################################
# 5) synchronizing species lists
########################################################################################

imbcr_species
sum(is.na(imbcr_species$aou))/nrow(imbcr_species)
sum(is.na(imbcr_species$alpha))/nrow(imbcr_species)

bbs_species
sum(is.na(bbs_species$alpha))/nrow(bbs_species)
sum(is.na(bbs_species$aou))/nrow(bbs_species)

western_species <- read_csv("data/species_lists/western_species.csv")
names(bbs_species)
names(imbcr_species)
names(ebird_key)

masterlist_species <- full_join(
  bbs_species %>%
    mutate(BBS = T),
  imbcr_species) %>%
  left_join(ebird_key %>%
              mutate(ebird = T), 
            by = c("en_common_name", "species"), 
            relationship = "many-to-many") %>%
  left_join(prop_species) %>%
  select(-c(alpha.y, binomial.y, order.y, family.y, genus.y))

masterlist_species %>% 
  filter(BBS == T & ebird == T & IMBCR > 0) %>%
  arrange(ebird, desc(IMBCR), BBS) %>%
  distinct() %>%
  View()

View(masterlist_species)

nrow(prop_species)
########################################################################################




########################################################################################
# 6) spatial filtering of IMBCR
########################################################################################

imbcr_sites <- imbcr %>%
  group_by(transect) %>%
  summarize(start_date = min(date, na.rm = T), 
            num_species = n_distinct(species),
            num_obs = n())

### filtering sites:
hist(imbcr_sites$num_species, breaks = 100)
hist(imbcr_sites$num_obs, breaks = 100)
hist(imbcr_sites$start_date, breaks = 100)

min_species <- 2
min_obs <- 10

sum(imbcr_sites$num_species < min_species)
sum(imbcr_sites$num_obs < min_obs)
sum(imbcr_sites$start_date > as.Date("2023-01-01"))

imbcr_fltr <- imbcr %>%
  filter(!transect %in% 
    imbcr$transect[imbcr_sites$num_species < min_species]) %>%
  filter(!transect %in% 
    imbcr$transect[imbcr_sites$num_obs < min_obs]) %>%
  filter(!transect %in% 
    imbcr$transect[imbcr_sites$start_date > as.Date("2023-01-01")]) 

imbcr_sites_fltr <- imbcr_sites %>%
  filter(!transect %in% unique(imbcr_fltr$transect))
imbcr_fltr
### BBS data

read_csv("data/BBS/BBS_sw/Arizona.csv")

bbs_raw <- lapply(list.files("data/BBS/BBS_sw"), FUN = function(x){
  read_csv(paste0("data/BBS/BBS_sw/", x)) %>%
    as.data.frame() %>%
    mutate(StateNum = as.double(StateNum))}) %>%
  bind_rows() %>%
  clean_names() %>%
  left_join(
    data.frame(
      state = c("Arizona", "California", "Colorado", "Idaho", "Iowa",
                "Kansas", "Minnesota", "Missouri", "Montana", "North Dakota",
                "Nebraska", "Nevada", "New Mexico", "Oklahoma", "Oregon",
                "South Dakota", "Texas", "Utah", "Washington", "Wyoming"), 
      state_num = c(6, 14, 17, 33, 36, 38, 50, 52, 53, 64, 54, 55, 60, 67, 69,
                    81, 83, 85, 89, 92))) 
bbs_raw$aou %>% unique()

bbs_raw %>%
  left_join(bbs_species %>%
              select(aou, en_common_name, binomial) %>%
              mutate(aou = as.character(aou)))
  
bbs_routes <- read_csv("data/BBS/routes.csv") %>%
  clean_names() %>%
  mutate(state_num = as.double(state_num))

region %>%
  filter(admin == "United States of America") %>%
  names()

### visualizing data
pal <- wes_palette("Zissou1", 5, type = "discrete")

gg_start_date <- ggplot() +
  geom_sf(data = region) +
  geom_sf(data = imbcr_sites_fltr, shape = 3, aes(col = as.Date(start_date))) +
  scale_colour_gradientn("start date", colors = pal, trans = "date",
                         breaks = scales::breaks_width("2 years")) +
  theme_void() +
  ggtitle("IMBCR: transect establishment date") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "inside",
        legend.position.inside = c(0.1, 0.2),
        legend.key.height=unit(0.75, "cm"))

gg_species_count <- ggplot() +
  geom_sf(data = region, fill = "grey90") +
  geom_sf(data = imbcr_sites_fltr, shape = 3, aes(col = num_species)) +
  scale_color_gradientn("species counts", 
                        colors = pal,
                        limits = c(0, 50),
                        breaks = seq(0, 50, 10)) +
  theme_void() +
  ggtitle("IMBCR: transect species counts") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "inside",
        legend.position.inside = c(0.1, 0.2),
        legend.key.height=unit(0.75, "cm"))

ggpubr::ggarrange(
  gg_start_date, 
  gg_species_count, 
  nrow = 2); ggsave("figures/2_panel_sites.jpeg", width = 12, height = 6, unit = "in")

ggplot() +
  geom_sf(data = region) +
  geom_point(data = imbcr %>% 
               filter(species %in% head(imbcr_summ$species, 10)), 
             aes(x = lon, y = lat, col = species),
             shape = 3) +
  ggtitle("IMBCR: Top 10 must observed bird species") +
  guides(col = guide_legend(
    direction = "horizontal", 
    nrow = 3)) +
  facet_wrap(~year) +
  theme_void() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.7, 0.1),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

ggsave("figures/top_10_species.jpeg", 
       width = 10, height = 12, unit = "in")
