pacman::p_load(
  tidyverse
  , here
  , auk
  , lubridate
  , sf
  # , 
  # , 
)

crs <- read_rds("data/region/crs.rds")
strat <- read_rds("data/region/strat.rds")
taxa <- read_rds("data/species/focal_species.rds")

# paths to ebd / checklist data
obs_file <- here("data/ebd/ebd/observations_clean.txt")
chk_file <- here("data/ebd/ebd/checklists_clean.txt")

# Reduce global data to western usa
ebd_region <- auk_ebd(
    file = obs_file,
    file_sampling = chk_file) %>%
  auk_species(species = taxa$en_common_name) %>%
  auk_bbox(bbox = c(-128.41723, 24.75365, -98.37560, 54.61227)) %>%
  auk_duration(c(0, 300)) %>%
  auk_protocol(c("Stationary")) %>% 
  auk_complete() %>%
  auk_filter(
    file = "data/ebd/filters/ebd_observations_prelim_sp_region.txt", 
    file_sampling = "data/ebd/filters/ebd_checklists_prelim_sp_region.txt",
    overwrite = TRUE)

ebd_region_observations <- read_ebd("data/ebd/filters/ebd_observations_prelim_sp_region.txt")
checklist_ids <- unique(ebd_region_observations$sampling_event_identifier)

ebd_region_checklists <- read_sampling("data/ebd/filters/ebd_checklists_prelim_sp_region.txt") %>%
  distinct(sampling_event_identifier, .keep_all = TRUE) %>%
  filter(sampling_event_identifier %in% 
    unique(ebd_region_observations$sampling_event_identifier))

ebd_region_observations <- ebd_region_observations %>%
  filter(sampling_event_identifier %in% 
    unique(ebd_region_checklists$sampling_event_identifier))

# zerofill and compile species counts
# ebd_region_zerofill <- auk_zerofill(
#   x = ebd_region_observations, 
#   sampling_events = ebd_region_checklists, 
#   collapse = TRUE
# ); write_rds(ebd_region_zerofill, "data/ebd/ebd_region_prelim_zerofill.rds")
ebd <- read_rds("data/ebd/ebd_region_prelim_zerofill.rds")

#######################################################################
################## checking which species to process ##################
#######################################################################

ebd_species <- taxa %>%
  filter(binomial %in% 
    ebd$scientific_name %>% 
    unique())

go_list <- c(
  lapply(ebd_species$alpha, FUN = function(x){
    file.exists(paste0("data/ebd/species/point/", x, ".rds") ) %>%
      isFALSE(.) == T}) %>%
    unlist() %>%
    which(), 
  lapply(ebd_species$alpha, FUN = function(x){
    file.exists(paste0("data/ebd/species/strat/", x, ".rds") ) %>%
      isFALSE(.) == T}) %>%
    unlist() %>%
    which()) %>%
  unique(); ebd_species$en_common_name[go_list]
  
# go_list <- go_list[-c(which(go_list %in% c(
#     19   # this one only has 1 point in Alaska
#   , 98   # this one is an eastern bird
#   , 107  # this one is an Alaskan bird
#   # ,  
# )))]; taxa$en_common_name[go_list]

##############################################################################

for(i in go_list){

  ebd_point <- ebd %>%
      filter(scientific_name == taxa$binomial[i]) %>%
      filter(!observation_count == "X")
  
  if(nrow(filter(ebd_point, count > 0)) >= 30){
      
    cat(paste0(
      "Iteration #: ", i, 
      "| Species: ", 
      ebd_species$en_common_name[i], 
      " has sufficient observations | continuing processing \n"))

    point_data <- ebd_point %>%
      st_as_sf(
        coords = c("longitude", "latitude"), 
        crs = st_crs(4326)) %>%
      st_transform(crs) %>%
      select(
        checklist_id, country, state, county, locality_id, 
        observation_date, observer_id, group_identifier, sampling_event_identifier, 
        observation_type, observation_count, duration_minutes, 
        effort_distance_km, effort_area_ha, number_observers, 
        breeding_category, age_sex) %>%
      rename(
        breeding = breeding_category,
        date = observation_date,
        observer = observer_id,
        group_id = group_identifier,
        observation = observation_type,
        distance = effort_distance_km,
        area = effort_area_ha,
        nobservers = number_observers) %>%
      mutate(
        count = as.integer(observation_count),
        year = year(date),
        month = lubridate::month(date),
        day = lubridate::day(date)
      ); write_rds(
          point_data, 
          paste0(
            "data/ebd/species/point/", 
            taxa$alpha[i], 
            ".rds")
          )
        
    # ex_strat_sum <- point_data %>%
    #   st_join(strat, join = st_within) %>%  
    #   group_by(cell_id, year) %>%
    #   summarize(count = sum(count, na.rm = T))
    
    # ex_ebd_strat <- strat %>%
    #   left_join(st_drop_geometry(ex_strat_sum), by = "cell_id") %>%
    #   filter(!count == 0)
    # write_rds(ex_ebd_strat, paste0("data/ebd/species/strat/", taxa$alpha[i], ".rds"))

      }

  cat(paste0(
    "Iteration #: ", i, 
    "| Species: ", 
    ebd_species$en_common_name[i], 
    " has insufficient observations \n"))
}

