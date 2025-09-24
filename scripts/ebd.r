pacman::p_load(
  tidyverse
  , here
  , auk
  , lubridate
  # , 
  # , 
)

crs <- read_rds("data/region/crs.rds")
strat <- read_rds("data/region/strat.rds")

obs_file <- here("data/ebd/ebd/observations.txt")
chk_file <- here("data/ebd/ebd/checklists.txt")

species <- read_rds("data/species/prop_prelim.rds")

for(i in 1:nrow(species)){

  tryCatch({
      ex_raw <- obs_file %>%
        auk_ebd() %>%
        auk_species(species = paste(species$en_common_name[i])) %>%
        auk_country(country = "US") %>%
        auk_filter(
          file = paste("data/ebd/filters/", species$alpha[i], ".txt"), 
          overwrite = TRUE) %>%
        read_ebd()
      
      ex <- ex_raw %>%
        st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(4326)) %>%
        st_transform(crs) %>%
        select(
          taxonomic_order, category, common_name, scientific_name,
          observation_count, breeding_category, behavior_code, 
          observation_date, observer_id, observation_type, effort_distance_km, 
          effort_area_ha, number_observers, group_identifier) %>%
        rename(
          order = taxonomic_order, 
          en_common_name = common_name,
          binomial_name = scientific_name,
          count = observation_count,
          breeding = breeding_category,
          behavior = behavior_code,
          date = observation_date,
          observer = observer_id,
          observation = observation_type,
          distance = effort_distance_km,
          area = effort_area_ha,
          nobservers = number_observers,
          group_id = group_identifier) %>%
        mutate(
          year = lubridate::year(date),
          month = lubridate::month(date),
          day = lubridate::day(date))
        write_rds(ex, paste("data/ebd/species/point/", species$alpha[i], ".rds"))
        
        ex_strat_sum <- ex %>%
          mutate(count = as.numeric(count)) %>%
          st_join(strat, join = st_within) %>%  
          group_by(strat_name, year) %>%
          summarize(count = sum(count, na.rm = T))
        
        ex_ebd_strat <- strat %>%
          left_join(st_drop_geometry(ex_strat_sum), by = "strat_name") %>%
          filter(!count == 0)
        write_rds(ex_ebd_strat, paste("data/ebd/species/strat/", species$alpha[i], ".rds"))

  }, error = function(e){
    # Baby, what's wrong?
    cat(paste("Error at species", species$en_common_name[i], ":", conditionMessage(e), "\n"))
     # Moving on...
    next
  })
    cat(paste("All done with species: ", species$en_common_name[i]))
}
