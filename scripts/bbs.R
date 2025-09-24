pacman::p_load(
  tidyverse
  , here
  , sf
  , bbsBayes
  # , 
  # , 
)

# fetch_bbs_data()

crs <- read_rds("data/region/crs.rds")
strat <- read_rds("data/region/strat.rds")

# bbs_strat <- bbsBayes::stratify(by = "latlong")
write_rds(bbs_strat, "data/bbs/bbs_strat.rds")

species <- read_rds("data/species/prop_prelim.rds")

for(i in 1:nrow(species)){
  
  tryCatch({
    # I guess this needs to be within the loop
    # bbs_strat <- bbsBayes::stratify(by = "latlong")
    bbs_strat <- read_rds(data/bbs/bbs_strat.rds) # preload it instead??

    # gather bbs point data
    bbs_raw <- bbsBayes::prepare_data(
      bbs_strat,
      species_to_run = paste0(species$en_common_name[i]),
      model = "slope")
    
    # clean it up
    bbs_df <- data.frame(
      count = bbs_raw$count,
      strat_name = bbs_raw$strat_name, 
      strat = bbs_raw$strat_name,
      obser = bbs_raw$obser,
      firstyr = bbs_raw$firstyr, 
      route = bbs_raw$route, 
      year = bbs_raw$r_year, 
      month = bbs_raw$month, 
      day = bbs_raw$day
    ); write_rds(bbs_df, paste0("data/bbs/species/point/", species$alpha[i], ".rds"))
    
    # stratify
    bbs_strat <- bbs_df %>%
      group_by(strat_name, year) %>%
      summarize(count = sum(count, na.rm = T)) %>%
      left_join(strat) %>%
      st_set_geometry("geometry")
    write_rds(bbs_strat, paste0("data/bbs/species/strat/", species$alpha[i], ".rds"))

  }, error = function(e){
    # Baby, what's wrong?
    cat(paste("Error at species", species$en_common_name[i], ":", conditionMessage(e), "\n"))
     # Moving on...
    next
  })
    cat(paste("All done with species: ", species$en_common_name[i]))
}
