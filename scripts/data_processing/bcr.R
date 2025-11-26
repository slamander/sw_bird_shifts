pacman::p_load(
  tidyverse
  , sf
  , here
  , janitor
  # ,
  # ,
)

crs <- read_rds("data/region/crs.rds")
taxa <- read_rds("data/species/focal_species.rds")
bird_key <- read_rds("data/species/bird_key.rds")
strat <- read_rds("data/region/strat.rds")

# Process and clean IMBCR data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bcr <- read_csv("data/bcr/sw_birds.csv") %>%
  clean_names() %>%
  rename(
    transect = transect_num, 
    observer = transect_visit_observer, 
    distance = radial_distance,
    count = cl_count,
    start_time = point_visit_start_time, 
    latitude = point_latitude, 
    longitude = point_longitude, 
    en_common_name = species) %>%
  mutate(
    date = as.Date(date),
    alpha = tolower(bird_code)
  ) %>% 
  left_join(bird_key) %>%
  filter(!is.na(alpha)) %>%
  filter(year %in% 2010:2015) %>%
  filter(!distance == -1) %>%
  filter(!is.na(distance)) %>%
  filter(!distance > 250) %>%
  filter(!time_period == -1) %>%
  arrange(alpha, year, transect, point, time_period) %>%
  select(
    year, county, state, latitude, longitude, date, 
    transect, point, distance, time_period, count,
    aou, alpha, en_common_name) 

# Select species to focus on ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bcr_alphas <- lapply(
  taxa$alpha, 
  FUN = function(x){
    grepl(
      x, 
      bcr %>%
        pull(alpha) %>% 
        tolower() %>%
        unique()) %>%
      sum(isTRUE(.))
  }) %>%
  unlist()
bcr_species <- taxa[which(bcr_alphas == 1),]

go_list <- c(
  lapply(bcr_species$alpha, FUN = function(x){
    file.exists(paste0("data/bcr/species/point/raw/", x, ".rds") ) %>%
      isFALSE(.) == T}) %>%
    unlist() %>%
    which(), 
  lapply(bcr_species$alpha, FUN = function(x){
    file.exists(paste0("data/bcr/species/point/dist/", x, ".rds") ) %>%
      isFALSE(.) == T}) %>%
    unlist() %>%
    which(), 
  lapply(bcr_species$alpha, FUN = function(x){
    file.exists(paste0("data/bcr/species/strat/raw/", x, ".rds") ) %>%
      isFALSE(.) == T}) %>%
    unlist() %>%
    which(), 
  lapply(bcr_species$alpha, FUN = function(x){
    file.exists(paste0("data/bcr/species/strat/dist/", x, ".rds") ) %>%
      isFALSE(.) == T}) %>%
    unlist() %>%
    which()
) %>% 
  unique(); taxa$en_common_name[go_list]
  
go_list <- go_list[-c(which(go_list %in% c(
    19   # this one only has 1 point in Alaska
  , 98   # this one is an eastern bird
  , 107  # this one is an Alaskan bird
  # ,  
)))]; taxa$en_common_name[go_list]

ggplot() + 
  geom_histogram(data = bcr, aes(distance))

# Create distance parameters for reshaping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
distance_breaks <- seq(0, 250, by = 50)

distance_bin_labels <- paste0(
  head(distance_breaks, -1), sep = "_", 
  tail(distance_breaks, -1))

time_periods <- 1:6

all_bin_columns <- expand.grid(
  distance_bin = distance_bin_labels,
  time_period = 1:6
) %>%
  mutate(col_name = paste0(
    "count", "_", 
    distance_bin, "_", 
    time_period)) %>%
  pull(col_name)

############ loop for 

for(i in go_list){

  bcr_point <- bcr %>%
    filter(alpha == bcr_species$alpha[i])
  
  if(nrow(filter(bcr_point, count > 0)) >= 30){
      
    cat(paste0(
      "Iteration #: ", i, 
      "| Species: ", 
      bcr_species$en_common_name[i], 
      " has sufficient observations | continuing processing \n"))

    write_rds(
      bcr_point,
        paste0(
        "data/bcr/species/point/raw/", 
        bcr_species$alpha[i], ".rds"))
  
    bcr_point_dist_bins <- bcr_point %>%
      mutate(
        distance_bin = cut(
          distance,
          breaks = distance_breaks, 
          labels = distance_bin_labels,
          include.lowest = TRUE, 
          right = FALSE)
      ) %>%
      mutate(
        bin_label = paste0("count", "_", distance_bin, "_", time_period)
      ) %>%
      group_by(year, point, transect, bin_label) %>%
      summarise(
        n_observations = sum(count), 
        .groups = "drop") %>%
      pivot_wider(
        names_from = bin_label,
        values_from = n_observations,
        values_fill = 0
      )

    missing_cols <- setdiff(all_bin_columns, names(bcr_point_dist_bins))
    bcr_point_dist_bins[missing_cols] <- 0

    bcr_point_dist <- bcr_point_dist_bins %>%
      select(
        year, transect, point, all_of(all_bin_columns), 
      ); write_rds(
        bcr_point_dist,
        paste0(
          "data/bcr/species/point/dist/", 
          bcr_species$alpha[i], 
          ".rds")
        )

    bcr_strat <- bcr_point %>%
      st_as_sf(
        coords = c("longitude", "latitude"),
        crs = "+proj=longlat +datum=WGS84"
      ) %>%
      bind_cols(st_coordinates(.)) %>%
      st_transform(st_crs(region)) %>%
      rename(latitude = "Y", longitude = "X") %>%
      st_intersection(strat) %>%
      st_drop_geometry() %>%
      full_join(strat) %>%
      st_set_geometry("x") %>%
      rename(geometry = x) %>%
      arrange(cell_id, year) %>%
      group_by(cell_id) %>%
      ungroup(); write_rds(
        bcr_strat,
          paste0(
            "data/bcr/species/strat/raw/", 
            bcr_species$alpha[i], 
            ".rds")
        )
    
    bcr_strat_dist_bins <- bcr_strat %>%
      mutate(
        distance_bin = cut(
          distance,
          breaks = distance_breaks, 
          labels = distance_bin_labels,
          include.lowest = TRUE, 
          right = FALSE)
      ) %>%
      mutate(
        bin_label = paste0("count", "_", distance_bin, "_", time_period)
      ) %>%
      group_by(year, point, transect, bin_label) %>%
      summarise(
        n_observations = sum(count), 
        .groups = "drop") %>%
      pivot_wider(
        names_from = bin_label,
        values_from = n_observations,
        values_fill = 0
      );
    
    missing_cols <- setdiff(all_bin_columns, names(bcr_strat_dist_bins))
    bcr_strat_dist_bins[missing_cols] <- 0

    bcr_strat_dist <- bcr_strat_dist_bins %>%
      select(
        year, transect, point, all_of(all_bin_columns), 
      ); write_rds(
        bcr_strat_dist,
          paste0(
            "data/bcr/species/strat/dist/", 
            bcr_species$alpha[i], 
            ".rds")
        )
      }

  cat(paste0(
    "Iteration #: ", i, 
    "| Species: ", 
    bcr_species$en_common_name[i], 
    " has insufficient observations \n"))
}

