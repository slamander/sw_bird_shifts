rm(list = ls())
pacman::p_load(
  tidyverse
  , here
  , sf
  , bbsBayes
  , janitor
  , rnaturalearth
  # , 
  # , 
)

crs <- read_rds("data/region/crs.rds")
taxa <- read_rds("data/species/focal_species.rds")
region <- read_rds("data/region/region.rds") 
strat <- read_rds("data/region/strat.rds")

bbs_point <- bbs_strat <- list()

bbs_raw <- load_bbs_data(level = "stop")

no_file_check <- function(files, path, ext){
    lapply(taxa$alpha, FUN = function(x)
    file.exists(paste0(path), x, ext)) %>%
      isFALSE(.) == T %>%
    unlist() %>%
    which()
    }

go_list <- c(
  lapply(taxa$alpha, FUN = function(x){
    file.exists(paste0("data/bbs/species/point/", x, ".rds") ) %>%
      isFALSE(.) == T}) %>%
    unlist() %>%
    which(), 
  lapply(taxa$alpha, FUN = function(x){
    file.exists(paste0("data/bbs/species/strat/", x, ".rds") ) %>%
      isFALSE(.) == T}) %>%
    unlist() %>%
    which()) %>%
  unique(); taxa$en_common_name[go_list]
  
go_list <- go_list[-c(which(go_list %in% c(
    19   # this one only has 1 point in Alaska
  , 98   # this one is an eastern bird
  , 107  # this one is an Alaskan bird
  # ,  
)))]; taxa$en_common_name[go_list]

bbs_point <- bbs_raw$bird %>%
  filter(AOU == as.numeric(taxa$aou[go_list[1]])) %>%
  left_join(bbs_raw$route) %>%
  clean_names() %>%
  filter(!is.na(longitude)) %>%
  filter(!is.na(latitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(4326)) %>%
  st_transform(crs = crs); ggplot() + 
  geom_sf(data = bbs_point)

for(i in go_list){

  # tryCatch({

  bbs_point <- bbs_raw$bird %>%
      filter(AOU == as.numeric(taxa$aou[i])) %>%
      left_join(bbs_raw$route) %>%
      clean_names() %>%
      filter(!is.na(longitude)) %>%
      filter(!is.na(latitude)) %>%
      st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(4326)) %>%
      st_transform(crs = crs) %>%
      st_intersection(strat) %>%
      bind_cols(st_coordinates(.)) %>%
      rowwise() %>%
      mutate(
        count1 = sum(c_across(stop1:stop5)), 
        count2 = sum(c_across(stop6:stop10)), 
        count3 = sum(c_across(stop11:stop15)), 
        count4 = sum(c_across(stop16:stop20)), 
        count5 = sum(c_across(stop21:stop25)), 
        count6 = sum(c_across(stop26:stop30)), 
        count7 = sum(c_across(stop31:stop35)), 
        count8 = sum(c_across(stop36:stop40)), 
        count9 = sum(c_across(stop41:stop45)), 
        count10 = sum(c_across(stop46:stop50)), 
        total_count = sum(c_across(count1:count10))
      ) %>%
      ungroup() %>%
      mutate(
        route_id = as.numeric(paste0(cell_id, route)), 
        year = as.numeric(year)) %>%
      arrange(route_id, year) %>%
      group_by(route_id, obs_n) %>%
      mutate(new = ifelse(row_number() == 1, 1, 0)) %>%
      ungroup() %>%
      arrange(route_id, year) %>%
      mutate(
        obs_id = as.numeric(factor(obs_n, levels = unique(obs_n)))
      ) %>%
      rename(
        latitude = Y, 
        longitude = X
      ) %>%
      select(
        state, route_id, cell_id, year, new, 
        count1, count2, count3, count4, count5, count6, 
        count7, count8, count9, count10, total_count,
        obs_id, geometry
      ); write_rds(
          bbs_point, 
          paste0("data/bbs/species/point/", taxa$alpha[i], ".rds")
      )

    bbs_strat <- bbs_point %>%
      distinct(
        state, route_id, cell_id, year, new, 
        count1, count2, count3, count4, count5, count6, 
        count7, count8, count9, count10, total_count,
        obs_id
      ) %>%
      st_drop_geometry() %>%
      full_join(strat) %>%
      st_set_geometry("x") %>%
      rename(geometry = x) %>%
      mutate(is_observed = !is.na(total_count) & total_count > 0) %>%  
      group_by(cell_id) %>%
      mutate(
        time_step = row_number(),
        nyears_site = n()) %>%
      ungroup(); write_rds(
        bbs_strat, 
        paste0("data/bbs/species/strat/", taxa$alpha[i], ".rds")
      )

  # }, error = function(e){

  #   # Baby, what's wrong?
  #   cat(paste("Error at species", taxa$en_common_name[i], ":", conditionMessage(e), "\n"))
  #    # Moving on...
  #   next
  # })
    cat(paste0("All done with species # ", i, ": ", taxa$en_common_name[i]))
}

bbs_point
