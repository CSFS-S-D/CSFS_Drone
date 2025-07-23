# Project     gappyClumpyFunctions.R
# Purpose     Gappy-clumpy thinning functions for use with gappyClumpyThin.R
# Person      Andy Whelan
# Date        December 9, 2024
# Modified    March 4, 2025
################################################################################


################################################################################
######                 Functions                                         #######
################################################################################

#---> Create a forest/non-forest raster from a point cloud
forested = function(las, w=5, res=5, n=5) {
  # las: las point cloud with forest
  # w:   focal window size for smooting
  # res: voxel and pixel resolution for point thinning and rasterization.
  # n:   Number of points to select per voxel
  forested_rast = las %>% filter_poi(Classification!=2) %>% 
    decimate_points(random_per_voxel(res=res, n=n)) %>% 
    rasterize_density(res=res) %>% 
    focal(w=w, "median")
  return(forested_rast)
}


poly_split = function(poly, ttops) {
  # split prod into smaller polygons if there are more than 30000 trees in ttops_temp
  ttops_temp =  st_intersection(ttops,poly)
  
  if(nrow(ttops_temp)<=30000) {
    polys = poly
  }else{
    grps = ceiling(nrow(ttops_temp)/30000)
    polys = divider(vect(poly), grps) %>% 
      st_as_sf() %>% 
      st_cast("MULTIPOLYGON") %>% 
      st_cast("POLYGON")
  }
  return(polys)
}

#---> Thin the forest controlled with some groupy-gappy parameters
forestThin_gg = function(prod, 
                         ttops, 
                         tree_clump_dist_m, 
                         ostory_dbh_cm,
                         target_ba,
                         target_qmd,
                         target_pcts) {
  
  
  thinned=list()
  for(i in 1:nrow(prod)) {
    
    # define ttops as the trees in one polygon from ttops files above.
    ttops_temp =  st_intersection(ttops,prod[i,])
    
    if(nrow(ttops_temp)==0) next
    if(!any(ttops_temp$dbh_cm>=ostory_dbh_cm)) next
    
    stand_sf = prod[i,]
    
    ################################################################################
    #####                     Define Clumps                                    #####
    ################################################################################
    
    # call get_tree_clumps
    ttops_temp = get_tree_clumps(
      ttops
      , stand_sf
      , tree_clump_dist_m = tree_clump_dist_m
      , ostory_dbh_cm = ostory_dbh_cm # 6 in dbh
    )
    
    ################################################################################
    #####                  Clump polygons and metrics                          #####
    ################################################################################
    
    # Get the clump summary data
    clump_polys_temp = get_clump_summary(ttops_temp)
    
    
    
    ################################################################################
    #####                        Clump Spacing                                 #####
    ################################################################################
    
    # Get a clump distance raster
    dist_rast_temp = get_clump_dist_rast(clump_polys_temp,stand_sf)
    
    
    ################################################################################
    #####                     Clump metrics                                    #####
    ################################################################################
    
    # Summarize by number of tree clump grouping variable
    clump_n_trees_grp_summary_temp = get_clump_n_trees_grp_summary(
      trees=ttops_temp, 
      stand_sf,
      clumps = get_clump_summary(ttops_temp)
    )
    
    
    ################################################################################
    #####                           ICO implementation                         #####
    ################################################################################
    
    if(nrow(clump_n_trees_grp_summary_temp)==1 & clump_n_trees_grp_summary_temp$n_trees[1]<=1) next
    if(max(clump_n_trees_grp_summary_temp$basal_area_ft2_per_ac)<=target_ba) next
    #if(all(get_tpa(target_ba,target_qmd)>clump_n_trees_grp_summary_temp$trees_per_ac)) next
    
    # Get targets for clump proportions
    target_data_temp = get_target_check_prescription(
      clump_n_trees_grp_summary_temp
      , target_ba = target_ba
      , target_qmd = target_qmd
      , target_pcts = target_pcts
    )  
    
    if(all(target_data_temp$stand_n_clumps_target==0)) next
    
    #--- Combine clump and opening targets with leave tree criteria into marking guidelines
    
    # what is the smallest group with a target?
    sm_grp_temp = target_data_temp %>% 
      dplyr::filter(pct_stand_n_trees_target>0) %>% 
      dplyr::pull(clump_n_trees_grp) %>% 
      min()
    
    # start with the largest group in target data currently with trees
    grp_temp = target_data_temp %>% 
      dplyr::arrange(desc(clump_n_trees_grp)) %>% 
      dplyr::filter(pct_stand_n_trees_current>0) %>% 
      dplyr::slice(1) %>% 
      dplyr::pull(clump_n_trees_grp)
    
    message(
      paste("started cutting", grp_temp, "at", Sys.time())
    )
    
    # identify clumps to leave untouched
    # clump_polys_temp = get_clump_summary(ttops_temp)
    keep_clumps_temp = 
      clump_polys_temp %>% 
      sf::st_drop_geometry() %>% 
      #keep only trees in current group
      dplyr::inner_join(
        target_data_temp %>% 
          dplyr::filter(clump_n_trees_grp == grp_temp) %>% 
          dplyr::select(
            clump_n_trees_grp, mean_clump_n_trees, min_clump_n_trees, max_clump_n_trees
            , stand_n_clumps_target
          )
        , by = "clump_n_trees_grp"
      ) %>% 
      # keep only clumps that meet criteria
      dplyr::filter(
        n_trees >= min_clump_n_trees
        & n_trees <= max_clump_n_trees
      ) %>% 
      # keep clumps closest the mean number in target
      dplyr::mutate(
        pct_to_target = abs( (n_trees-mean_clump_n_trees)/mean_clump_n_trees )
      ) %>% 
      dplyr::arrange(pct_to_target, desc(qmd_cm), n_trees, clump_id) %>% 
      dplyr::filter(dplyr::row_number()<=stand_n_clumps_target) %>% 
      dplyr::select(clump_id)
    
    # !!!!!!!!!!FIX: WHAT IF WE HAVE CLUMPS OF THIS SIZE BUT THE N_TREES>MAX_TREES AND NEED TO GET CLUMPS OF THIS SIZE?
    #   ... UPDATE CUTTING ALG TO CUT CLUMP DOWN TO MIN-MAX TREE RANGE FOR CLUMPS WITH INF UPPER LIMIT
    
    # start building tree list with keep/cut flag
    # ttops_temp = get_tree_clumps(
    #   my_suid = my_suid
    #   , tree_clump_dist_m = tree_clump_dist_m
    #   , ostory_dbh_cm = ostory_dbh_cm
    # )
    # build tree list 
    keep_trees =
      ttops_temp %>% 
      dplyr::ungroup() %>% 
      sf::st_drop_geometry() %>% 
      dplyr::inner_join(keep_clumps_temp, by = "clump_id") %>% 
      dplyr::select(treeID) %>% 
      dplyr::mutate(is_keep_tree = 1)
    
    # did we keep all of the clumps?
    if( 
      nrow(keep_clumps_temp) == 
      ( clump_polys_temp %>% dplyr::filter(clump_n_trees_grp == grp_temp) %>% nrow() )
    ){
      remaining_trees = dplyr::tibble(
        treeID = character(0)
        , is_keep_tree = numeric(0)
      )
    }else{
      # determine keep/cut for the remaining trees in that group type 
      remaining_trees = 
        ttops_temp %>%
        dplyr::ungroup() %>% 
        # keep only trees in current grp
        dplyr::filter(clump_n_trees_grp == grp_temp) %>% 
        # remove keep trees
        dplyr::anti_join(keep_trees, by = "treeID")
      
      # apply the cutting algorithm
      # get the flag
      remaining_trees$is_keep_tree = get_keep_tree_flag(
        x = remaining_trees
        , tgt = target_data_temp %>%  
          dplyr::ungroup() %>%
          dplyr::arrange(clump_n_trees_grp) %>%
          dplyr::mutate(l = dplyr::lag(clump_n_trees_grp)) %>%
          dplyr::filter(clump_n_trees_grp == grp_temp) %>%
          dplyr::pull(l)
        , clump_dist_m = tree_clump_dist_m
      ) 
      
      # select relevant columns
      remaining_trees = remaining_trees %>% 
        sf::st_drop_geometry() %>% 
        dplyr::ungroup() %>% 
        dplyr::select(treeID, is_keep_tree)
    }
    message(
      paste("done cutting", grp_temp, "at", Sys.time())
    )
    
    ###############################################
    # now process to go on to the next groups
    ###############################################
    # get the next group
    grp_temp = target_data_temp %>%
      dplyr::filter(clump_n_trees_grp<grp_temp) %>% 
      dplyr::arrange(desc(clump_n_trees_grp)) %>% 
      dplyr::slice(1) %>% 
      dplyr::pull(clump_n_trees_grp)
    
    while(grp_temp>=sm_grp_temp & grp_temp != "Individual"){
      message(
        paste("started cutting", grp_temp, "at", Sys.time())
      )
      # identify clumps to leave untouched
      # clump_polys_temp = get_clump_summary(ttops_temp)
      keep_clumps_temp =
        clump_polys_temp %>% 
        sf::st_drop_geometry() %>% 
        #keep only trees in current group
        dplyr::inner_join(
          target_data_temp %>% 
            dplyr::filter(clump_n_trees_grp == grp_temp) %>% 
            dplyr::select(
              clump_n_trees_grp, mean_clump_n_trees, min_clump_n_trees, max_clump_n_trees
              , stand_n_clumps_target
            )
          , by = "clump_n_trees_grp"
        ) %>% 
        # keep only clumps that meet criteria
        dplyr::filter(
          n_trees >= min_clump_n_trees
          & n_trees <= max_clump_n_trees
        ) %>% 
        # keep clumps closest the mean number in target group
        dplyr::mutate(
          pct_to_target = abs( (n_trees-mean_clump_n_trees)/mean_clump_n_trees )
        ) %>% 
        dplyr::arrange(pct_to_target, desc(qmd_cm), n_trees, clump_id) %>% 
        dplyr::filter(dplyr::row_number()<=stand_n_clumps_target) %>% 
        dplyr::select(clump_id)
      
      # add to tree list with keep/cut flag
      keep_trees = keep_trees %>% 
        dplyr::bind_rows(
          ttops_temp %>% 
            dplyr::ungroup() %>% 
            sf::st_drop_geometry() %>% 
            dplyr::inner_join(keep_clumps_temp, by = "clump_id") %>% 
            dplyr::select(treeID) %>% 
            dplyr::mutate(is_keep_tree = 1)
        )
      
      # check if the desired clump number was met and add the remaining trees from previous group if needed
      more_clumps_target = 0
      if(
        nrow(keep_clumps_temp) <
        ( 
          target_data_temp %>% 
          dplyr::filter(clump_n_trees_grp == grp_temp) %>% 
          dplyr::pull(stand_n_clumps_target) 
        )
      ){
        # how many more clumps are needed?
        more_clumps_target = 
          ( 
            target_data_temp %>% 
              dplyr::filter(clump_n_trees_grp == grp_temp) %>% 
              dplyr::pull(stand_n_clumps_target) 
          ) - nrow(keep_clumps_temp)
        
        # determine group size of remaining trees in the group prior that were not in a group selected and were not cut
        potential_trees = 
          ttops_temp %>% 
          dplyr::ungroup() %>% 
          dplyr::inner_join(
            remaining_trees %>% dplyr::filter(is_keep_tree == 1) %>% dplyr::select(treeID)
            , by = "treeID"
          ) %>% 
          st_clump_points(clump_dist_m = tree_clump_dist_m) %>% 
          # ggplot() + geom_sf(aes(fill = clump_n_trees_grp)) + theme_void()
          # get the original clump id to prioritize new clump groups in the same area
          dplyr::inner_join(
            ttops_temp %>% 
              sf::st_drop_geometry() %>% 
              dplyr::ungroup() %>% 
              dplyr::select(treeID, clump_id) %>% 
              dplyr::rename(orig_clump_id = clump_id)
            , by = "treeID"
          ) %>% 
          dplyr::group_by(orig_clump_id) %>% 
          dplyr::mutate(
            pct_desired_grp = 
              sum(ifelse(clump_n_trees_grp == grp_temp, 1, 0)) / dplyr::n()
          ) %>% 
          dplyr::ungroup() %>% 
          # keep only the current group
          dplyr::filter(clump_n_trees_grp == grp_temp)
        
        # pick trees from potential trees based on clump summary
        keep_trees = keep_trees %>% 
          dplyr::bind_rows(
            potential_trees %>% 
              sf::st_drop_geometry() %>% 
              # filter trees based on clumps needed
              dplyr::inner_join(
                get_clump_summary(potential_trees) %>% 
                  sf::st_drop_geometry() %>% 
                  # join with original clump id metrics to prioritize 
                  # keeping clumps in the same area and minimize cutting time
                  dplyr::left_join(
                    potential_trees %>% 
                      sf::st_drop_geometry() %>% 
                      dplyr::group_by(clump_id, orig_clump_id) %>% 
                      dplyr::summarise(pct_desired_grp = max(pct_desired_grp))
                    , by = "clump_id"
                  ) %>% 
                  # keep the number of clumps needed
                  dplyr::arrange(desc(pct_desired_grp), orig_clump_id, desc(qmd_cm), n_trees, clump_id) %>% 
                  dplyr::filter(dplyr::row_number()<=more_clumps_target) %>%
                  dplyr::select(clump_id)
                , by = "clump_id"
              ) %>% 
              dplyr::select(treeID) %>% 
              dplyr::mutate(is_keep_tree = 1)
          )
      } # end if don't have enough clumps
      ####################################################
      # determine keep/cut for the remaining trees in that group type and higher groups 
      ####################################################
      remaining_trees =
        ttops_temp %>%
        dplyr::ungroup() %>% 
        # REMOVE TREES FROM PREVIOUS TREES REMAINING that got cut to make this group size
        dplyr::anti_join(
          remaining_trees %>% dplyr::filter(is_keep_tree == 0)
          , by = "treeID"
        ) %>%
        # keep only trees in current grp or prior grp
        dplyr::filter(clump_n_trees_grp >= grp_temp) %>% 
        # remove keep trees
        dplyr::anti_join(keep_trees, by = "treeID")
      
      if(nrow(remaining_trees) > 0){
        # apply the cutting algorithm
        # get the flag
        remaining_trees$is_keep_tree = get_keep_tree_flag(
          x = remaining_trees
          , tgt = target_data_temp %>%  
            dplyr::ungroup() %>%
            dplyr::arrange(clump_n_trees_grp) %>%
            dplyr::mutate(l = dplyr::lag(clump_n_trees_grp)) %>%
            dplyr::filter(clump_n_trees_grp == grp_temp) %>%
            dplyr::pull(l)
          , clump_dist_m = tree_clump_dist_m
        ) 
        
        # select relevant columns
        remaining_trees = remaining_trees %>% 
          sf::st_drop_geometry() %>% 
          dplyr::ungroup() %>% 
          dplyr::select(treeID, is_keep_tree)
      }
      
      # increment
      message(
        paste("done cutting", grp_temp, "at", Sys.time())
      )
      # get the next group
      grp_temp = target_data_temp %>%
        dplyr::filter(clump_n_trees_grp<grp_temp) %>% 
        dplyr::arrange(desc(clump_n_trees_grp)) %>% 
        dplyr::slice(1) %>% 
        dplyr::pull(clump_n_trees_grp) %>% 
        dplyr::coalesce("Individual")
    }
    
    
    # now individuals
    if(nrow(remaining_trees) > 0) {
      if(sm_grp_temp == "Individual"){
        message(
          paste("started cutting", sm_grp_temp, "at", Sys.time())
        )
        potential_trees = 
          # original data
          ttops_temp %>% 
          dplyr::filter(clump_n_trees_grp == "Individual") %>% 
          # remaining trees
          dplyr::bind_rows(
            ttops_temp %>% 
              dplyr::ungroup() %>% 
              dplyr::inner_join(
                remaining_trees %>% dplyr::filter(is_keep_tree == 1) %>% dplyr::select(treeID)
                , by = "treeID"
              )
          ) %>% 
          # make sure that these are all individuals
          dplyr::mutate(is_orig = ifelse(clump_n_trees_grp=="Individual", 1, 0)) %>% 
          st_clump_points(clump_dist_m = tree_clump_dist_m) %>% 
          dplyr::filter(clump_n_trees_grp == "Individual")
        
        # pick trees from potential trees based on target
        keep_trees = keep_trees %>% 
          dplyr::select(treeID) %>% 
          dplyr::bind_rows(
            potential_trees %>% 
              sf::st_drop_geometry() %>% 
              # keep the number of clumps needed
              dplyr::arrange(desc(is_orig), desc(dbh_cm), desc(tree_height_m)) %>% 
              dplyr::filter(
                dplyr::row_number() <=
                  target_data_temp %>% 
                  dplyr::filter(clump_n_trees_grp == "Individual") %>% 
                  dplyr::pull(stand_n_clumps_target)
              ) %>% 
              dplyr::select(treeID)
          ) %>% 
          dplyr::mutate(is_keep_tree = 1)
        
        message(
          paste("done cutting", sm_grp_temp, "at", Sys.time())
        )
      }
    }
    
    
    # get the final prescription
    prescription_trees = ttops_temp %>% 
      dplyr::ungroup() %>% 
      dplyr::left_join(
        keep_trees
        , by = "treeID"
      ) %>% 
      # tracking vars
      dplyr::mutate(
        # fill in keep tree flag
        is_keep_tree = dplyr::coalesce(is_keep_tree, 0)
        , orig_clump_n_trees_grp = clump_n_trees_grp
      ) %>% 
      dplyr::select(-c(clump_n_trees_grp))
    
    # attach the new clumping to the trees
    # first check if there are any keep trees
    if(all(prescription_trees$is_keep_tree==0)) next 
    
    prescription_trees = prescription_trees %>% 
      dplyr::left_join(
        # reclump
        prescription_trees %>% 
          dplyr::filter(is_keep_tree == 1) %>% 
          st_clump_points(clump_dist_m = tree_clump_dist_m) %>% 
          sf::st_drop_geometry() %>% 
          dplyr::select(treeID, clump_n_trees_grp) %>% 
          dplyr::rename(new_clump_n_trees_grp = clump_n_trees_grp)
        , by = "treeID"
      ) %>% 
      dplyr::mutate(
        new_clump_n_trees_grp = forcats::fct_na_value_to_level(new_clump_n_trees_grp, level = "Cut tree")
      )
    thinned[[length(thinned)+1]]=prescription_trees
  }
  
  return(do.call(rbind, thinned))
}



#---> Convert metric to imperial
calc_imperial_units_fn <- function(df) {
  df %>% 
    # convert to imperial units
    dplyr::mutate(
      dplyr::across(
        .cols = tidyselect::ends_with("_cm")
        , ~ .x * 0.394
        , .names = "{.col}_in"
      )
      , dplyr::across(
        .cols = tidyselect::ends_with("_m")
        , ~ .x * 3.28
        , .names = "{.col}_ft"
      )
      , dplyr::across(
        .cols = tidyselect::ends_with("_m2_per_ha")
        , ~ .x * 4.359
        , .names = "{.col}_ftac"
      )
      , dplyr::across(
        .cols = tidyselect::ends_with("_per_ha") & !tidyselect::ends_with("_m2_per_ha")
        , ~ .x * 0.405
        , .names = "{.col}_ac"
      )
      , dplyr::across(
        .cols = tidyselect::ends_with("_area_ha")
        , ~ .x * 2.471
        , .names = "{.col}_ac"
      )
      , dplyr::across(
        .cols = tidyselect::ends_with("_m2")
        , ~ .x * 10.764
        , .names = "{.col}_ft2"
      )
    ) %>%
    dplyr::rename_with(
      .fn = function(x){dplyr::case_when(
        stringr::str_ends(x,"_cm_in") ~ stringr::str_replace(x,"_cm_in","_in")
        , stringr::str_ends(x,"_m_ft") ~ stringr::str_replace(x,"_m_ft","_ft")
        , stringr::str_ends(x,"_m2_per_ha_ftac") ~ stringr::str_replace(x,"_m2_per_ha_ftac","_ft2_per_ac")
        , stringr::str_ends(x,"_per_ha_ac") ~ stringr::str_replace(x,"_per_ha_ac","_per_ac")
        , stringr::str_ends(x,"_area_ha_ac") ~ stringr::str_replace(x,"_area_ha_ac","_area_ac")
        , stringr::str_ends(x,"_m2_ft2") ~ stringr::str_replace(x,"_m2_ft2","_ft2")
        , TRUE ~ x
      )}
    )
}


#---> function to clump data that are sf points
st_clump_points <- function(
    x # point data of class `sf` 
    , clump_dist_m = 6 # size (radius) of the epsilon neighborhood = maximum distance between points to add to clump
    , clump_breaks = c(0,1,4,9,15,25,Inf) # where to break the clump size groups
    , clump_labels = c("Individual","2-4 trees","5-9 trees","10-15 trees","16-25 trees",">25 trees") # what to name the clump size groups
) {
  # get points as x,y
  point_clusters = x %>% 
    dplyr::mutate(
      X = sf::st_coordinates(.)[,1] %>% as.numeric()
      , Y = sf::st_coordinates(.)[,2] %>% as.numeric()
    )
  #############################################################################
  ##### Identify clusters in each stem map plot                           #####
  #############################################################################
  ### Place trees into clusters using an inter-tree distance of 6 m
  my_dbscan_temp =  point_clusters %>% 
    sf::st_drop_geometry() %>% 
    dplyr::select(X,Y) %>% 
    dbscan::dbscan(eps = clump_dist_m, minPts = 2)
  
  # my_dbscan_temp %>% str()
  
  ### append cluster ID to tree points
  point_clusters$dbscan_cluster = my_dbscan_temp$cluster
  # point_clusters$cluster %>% summary()
  # point_clusters %>% sf::st_drop_geometry() %>% dplyr::count(cluster) %>% dplyr::arrange(desc(n)) %>% dplyr::slice_head(n=11)
  
  ### cluster metrics
  point_clusters = point_clusters %>% 
    dplyr::group_by(dbscan_cluster) %>% 
    dplyr::mutate(
      # unique dbscan_cluster for individuals
      clump_id = dplyr::case_when(
        dbscan_cluster == 0 ~ max(my_dbscan_temp$cluster)+dplyr::row_number()
        , T ~ dbscan_cluster
      ) %>% 
        factor()
    ) %>% 
    dplyr::group_by(clump_id) %>% 
    dplyr::mutate(
      dbscan_cluster = factor(dbscan_cluster)
      , clump_n_trees = dplyr::n()
      , clump_n_trees_grp = cut(
        clump_n_trees
        , breaks = clump_breaks
        , labels = clump_labels
      ) %>% 
        factor(
          ordered = T
          , levels = clump_labels
        )
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(tree_clump_dist_m = clump_dist_m)
  # return
  return(point_clusters)
}


#---> Function to pass a unit id and return list of trees with clump groupings
get_tree_clumps = function(
    ttops
    , stand_sf
    , tree_clump_dist_m=6
    , ostory_ht_m = as.numeric(NA)
    , ostory_dbh_cm = as.numeric(NA)
){
  # check ostory definition
  if(is.na(as.numeric(ostory_dbh_cm)) & is.na(as.numeric(ostory_ht_m))){
    warning("`ostory_dbh_cm` and `ostory_ht_m` are not set...using `ostory_dbh_cm` = 12.7")
    ostory_dbh_cm = 5*2.54
    # filter data
    ttops_temp = ttops %>% sf::st_intersection(stand_sf) %>% 
      dplyr::filter(
        dbh_cm>=as.numeric(ostory_dbh_cm)
      )
  }else if(!is.na(as.numeric(ostory_dbh_cm))){
    # filter data
    ttops_temp = ttops %>% sf::st_intersection(stand_sf) %>% 
      dplyr::filter(
        dbh_cm>=as.numeric(ostory_dbh_cm)
      )
  }else{
    # filter data
    ttops_temp = ttops %>% sf::st_intersection(stand_sf) %>% 
      dplyr::filter(
        tree_height_m>=as.numeric(ostory_ht_m)
      )
  }
  
  ### Place trees into clusters using an inter-tree distance of 6 m
  ttops_temp = st_clump_points(
    x = ttops_temp
    , clump_dist_m = tree_clump_dist_m
  )
  
  # add distance to nearest within clump
  ttops_temp =
    ttops_temp %>% 
    dplyr::group_by(clump_id) %>%
    tidyr::nest() %>% 
    dplyr::mutate(
      distance_clump_nn_m = purrr::map(data, function(x){
        # get index of nearest neighbor
        i = sf::st_nearest_feature(x)
        # get dist
        d = sf::st_distance(x, x[i,], by_element=TRUE) %>% as.numeric()
        return(d)
      })
    ) %>% 
    tidyr::unnest(cols = c(data, distance_clump_nn_m)) %>% 
    sf::st_set_geometry("geom") %>% # set it cuz it got lost 
    dplyr::ungroup() %>% 
    dplyr::mutate(
      tree_clump_dist_m = tree_clump_dist_m
      , clump_id_duplicate = clump_id # can use this even after nesting data by clump_id
      # , ostory_ht_m = ifelse(is.na(ostory_ht_m), as.numeric(NA), as.numeric(ostory_ht_m))
      # , ostory_dbh_cm = ifelse(is.na(ostory_dbh_cm), as.numeric(NA), as.numeric(ostory_dbh_cm))
    )
  # return
  return(ttops_temp)
}


#---> Function to pass a return from st_clump_points/get_tree_clumps and create clump polygons with summary stats
get_clump_summary = function(dta){
  # get tree_clump_dist_m
  tree_clump_dist_m = min(dta$tree_clump_dist_m, na.rm = T)
  # create clump polys and summary
  clump_polys_temp = 
    dta %>% 
    dplyr::ungroup() %>% 
    sf::st_set_geometry("geometry") %>% 
    sf::st_buffer(tree_clump_dist_m/2) %>% 
    dplyr::group_by(clump_id, dbscan_cluster, clump_n_trees_grp) %>%
    dplyr::summarise(
      # union buffered tree points
      geometry = sf::st_union(geometry)
      # summary metrics
      , n_trees = dplyr::n_distinct(treeID)
      , mean_dbh_cm = mean(dbh_cm, na.rm = T)
      , mean_tree_height_m = mean(tree_height_m, na.rm = T)
      , loreys_height_m = sum(basal_area_m2*tree_height_m, na.rm = T) / sum(basal_area_m2, na.rm = T)
      , basal_area_m2 = sum(basal_area_m2, na.rm = T)
      , sum_dbh_cm_sq = sum(dbh_cm^2, na.rm = T)
      , .groups = "drop_last"
    ) %>%
    dplyr::ungroup() %>% 
    sf::st_make_valid() %>% 
    dplyr::mutate(
      clump_area_ha = sf::st_area(geometry) %>% as.numeric() %>% `/`(10000)
      , trees_per_ha = (n_trees/clump_area_ha)
      , basal_area_m2_per_ha = (basal_area_m2/clump_area_ha)
      , pct_stand_basal_area = basal_area_m2/sum(basal_area_m2)
      , pct_stand_n_trees = n_trees/sum(n_trees)
      , qmd_cm = sqrt(sum_dbh_cm_sq/n_trees)
    ) %>%
    dplyr::select(-c(sum_dbh_cm_sq)) %>% 
    # convert to imperial units
    calc_imperial_units_fn() %>% 
    dplyr::mutate(tree_clump_dist_m = tree_clump_dist_m)
  # calculate distance between clumps
  clump_polys_temp = clump_polys_temp %>% 
    dplyr::mutate(
      nearest = sf::st_nearest_feature(clump_polys_temp)
      , distance_nearest_clump_m = sf::st_distance(
        clump_polys_temp
        , clump_polys_temp[nearest,]
        , by_element=TRUE
      ) %>% 
        as.numeric()
    ) %>% 
    dplyr::select(-c(nearest))
  # return
  return(clump_polys_temp)
}


#---> create function to pass a return from get_clump_summary() and get a distance raster
get_clump_dist_rast = function(dta, stand_sf){
  # get tree_clump_dist_m
  tree_clump_dist_m = min(dta$tree_clump_dist_m, na.rm = T)
  # rasterize the clump polygons and then calculate distance between clumps as raster
  dist_rast = 
    terra::rasterize(
      x = dta %>% terra::vect()
      , y = dta %>% 
        terra::vect() %>% 
        terra::rast(res = 0.2)
    ) %>% 
    terra::distance() %>% 
    # crop it to stand extent
    terra::crop(
      stand_sf %>% 
        terra::vect()
    ) %>% 
    terra::mask(
      stand_sf %>% 
        terra::vect()
    )
  ######### part 2
  # now create openings vector data
  openings_vect = 
    dist_rast %>% 
    terra::classify(rcl = c(tree_clump_dist_m/2,Inf), others = NA, include.lowest = T) %>% 
    terra::as.polygons(na.rm = T) %>% 
    sf::st_as_sf() %>% 
    sf::st_cast("POLYGON") %>% 
    dplyr::mutate(layer = dplyr::row_number()) %>% 
    dplyr::mutate(
      openining_area_m2 = sf::st_area(geometry) %>% as.numeric()
      , tree_clump_dist_m = tree_clump_dist_m
    )
  
  # return
  return(list(dist_rast = dist_rast, openings_vect = openings_vect))
}

#---> create a function to summarize by number of tree clump grouping
get_clump_n_trees_grp_summary = function(trees, clumps, stand_sf){
  # get area of harvest unit
  #...will use this area in the area calculations such that...
  #...TPA = trees in a certain group size across the whole stand area
  harvest_area_ha = stand_sf %>% 
    sf::st_as_sf() %>% 
    dplyr::mutate(harvest_area_ha = sf::st_area(.) %>% as.numeric() %>% `/`(10000)) %>% 
    dplyr::pull(harvest_area_ha) %>% 
    .[1]
  # collapse and calculate silv metrics
  dta = trees %>% 
    sf::st_drop_geometry() %>% 
    dplyr::mutate(stand_area_ha = harvest_area_ha) %>% 
    dplyr::group_by(stand_area_ha,clump_n_trees_grp) %>%
    dplyr::summarise(
      # summary metrics
      n_trees = dplyr::n_distinct(treeID)
      , mean_dbh_cm = mean(dbh_cm, na.rm = T)
      , mean_tree_height_m = mean(tree_height_m, na.rm = T)
      , loreys_height_m = sum(basal_area_m2*tree_height_m, na.rm = T) / sum(basal_area_m2, na.rm = T)
      , basal_area_m2 = sum(basal_area_m2, na.rm = T)
      , sum_dbh_cm_sq = sum(dbh_cm^2, na.rm = T)
      , .groups = "drop_last"
    ) %>%
    dplyr::ungroup() %>% 
    # attach clump area
    dplyr::left_join(
      clumps %>% 
        sf::st_drop_geometry() %>% 
        dplyr::group_by(clump_n_trees_grp) %>%
        dplyr::summarise(
          clump_area_ha = sum(clump_area_ha)
          , stand_n_clumps = dplyr::n()
          , .groups = "drop_last"
        ) %>% 
        dplyr::ungroup()
      , by = dplyr::join_by(clump_n_trees_grp)
    ) %>% 
    dplyr::mutate(
      trees_per_ha = (n_trees/stand_area_ha) # (n_trees/clump_area_ha) ... this was not right
      , basal_area_m2_per_ha = (basal_area_m2/stand_area_ha) # (basal_area_m2/clump_area_ha) ... this was not right
      , qmd_cm = sqrt(sum_dbh_cm_sq/n_trees)
      # stand calcs
      , stand_trees_per_ha = sum(n_trees)/stand_area_ha
      , stand_basal_area_m2 = sum(basal_area_m2)
      , stand_basal_area_m2_per_ha = sum(basal_area_m2)/stand_area_ha
      , pct_stand_basal_area = basal_area_m2/stand_basal_area_m2
      , pct_stand_n_trees = n_trees/sum(n_trees)
      , stand_qmd_cm = sqrt(sum(trees$dbh_cm^2, na.rm = T)/sum(n_trees))
    ) %>%
    dplyr::select(-c(sum_dbh_cm_sq)) %>% 
    # convert to imperial units
    calc_imperial_units_fn()
  # return
  return(dta)
}


#---> function to get tpa from ba and qmd
get_tpa = function(ba_ft2_ac, qmd_in){
  tpa = round(ba_ft2_ac/((qmd_in^2)*0.005454))
  return(tpa)
} 


#---> Function to check, setup, and define data with targets based on Churchill et al. 2016
get_target_check_prescription = function(
    clump_n_trees_grp_summary_dta
    , target_ba = as.numeric(NA)
    , target_qmd = as.numeric(NA)
    , target_pcts = as.numeric(NA)
){
  if(
    is.na(target_ba) | is.na(target_qmd) | max(is.na(target_pcts))==1
  ){
    stop("must set all of the function parameters:\n`target_ba`, `target_qmd`, and `target_pcts`")
  }
  #############################################
  # check target BA and TPA
  #############################################
  if(as.numeric(target_ba)>clump_n_trees_grp_summary_dta$stand_basal_area_ft2_per_ac[1]){
    stop(
      "target BA in `target_ba` of "
      , round(as.numeric(target_ba),1), " is greater than current BA of "
      , clump_n_trees_grp_summary_dta$stand_basal_area_ft2_per_ac[1] %>% round(1)
    )
  }
  if(
    get_tpa(target_ba, target_qmd)>clump_n_trees_grp_summary_dta$stand_trees_per_ac[1]
  ){
    stop(
      "target TPA in of "
      , round(as.numeric(get_tpa(target_ba, target_qmd)),1), " is greater than current TPA of "
      , clump_n_trees_grp_summary_dta$stand_trees_per_ac[1] %>% round(1)
      , "\n adjust `target_ba` and/or `target_qmd` to get valid TPA"
    )
  }
  #############################################
  # define data with current and target
  # ... this is "smart" in that percentages are adj based on:
  # ... 0) are there missing targets?
  # ... ... if < 6 numbers provided in `target_pcts` then the largest tree groups get targets of 0
  # ... 1) do targets sum to 1? 
  # ... ... if not trees are distributed proportionally based on targets provided and trees available
  # ... 2) is target in largest clump size > current conditions?
  # ... ... if yes, target is set to current condition
  # ... 3) is target listed in clump size > current largest clump with trees?
  # ... ... if yes, target for largest clump size is shifted to current largest clump with trees
  #############################################
  target_data = 
    # create data for joining if missing clump groups
    dplyr::tibble(
      stand_area_ac = rep(clump_n_trees_grp_summary_dta$stand_area_ac[1],6)
      , clump_n_trees_grp = factor(
        c(1:6)
        , labels = c("Individual", "2-4 trees",   "5-9 trees",    "10-15 trees","16-25 trees",">25 trees")
        , ordered = T
      )
      , mean_clump_n_trees = c(1,3,7,12,20,30)
      , min_clump_n_trees = c(1,2,5,10,16,26)
      , max_clump_n_trees = c(1,4,9,15,25,99999)
    ) %>% 
    dplyr::left_join(
      clump_n_trees_grp_summary_dta %>% 
        dplyr::ungroup() %>% 
        dplyr::select(clump_n_trees_grp, pct_stand_n_trees, stand_n_clumps)
    ) %>% 
    dplyr::mutate(
      pct_stand_n_trees = dplyr::coalesce(pct_stand_n_trees,0)
      , stand_n_clumps = dplyr::coalesce(stand_n_clumps,0)
    ) %>% 
    # add targets
    dplyr::bind_cols(
      pct_stand_n_trees_target = c(as.numeric(target_pcts), rep(0,6))[1:6] # pad target with 0's
    ) %>% 
    # adjust target based on difference from 1
    dplyr::mutate(
      pct_stand_n_trees_target = pct_stand_n_trees_target*(1/sum(pct_stand_n_trees_target))
      # largest clump size with trees
      , largest_w_trees = max(ifelse(dplyr::coalesce(pct_stand_n_trees)>0,clump_n_trees_grp,NA),na.rm = T)
      , largest_w_trees_target = max(ifelse(dplyr::coalesce(pct_stand_n_trees_target)>0,clump_n_trees_grp,NA),na.rm = T)
    ) %>% 
    # move target for largest clump size to the largest current clump size 
    dplyr::mutate(
      pct_stand_n_trees_target = dplyr::case_when(
        as.numeric(clump_n_trees_grp)==largest_w_trees &
          largest_w_trees_target>largest_w_trees ~ max(
            ifelse(as.numeric(clump_n_trees_grp)==largest_w_trees_target,pct_stand_n_trees_target,0)
          )
        , T ~ pct_stand_n_trees_target
      )
    ) %>% 
    # adjust target based on current conditions
    dplyr::mutate(
      pct_stand_n_trees_target = dplyr::case_when(
        as.numeric(clump_n_trees_grp)>largest_w_trees &
          pct_stand_n_trees_target > 0 ~ 0
        , as.numeric(clump_n_trees_grp)==largest_w_trees &
          pct_stand_n_trees_target > pct_stand_n_trees ~ pct_stand_n_trees
        , T ~ pct_stand_n_trees_target
      )
    ) %>% 
    # finally, re-scale again based on adjustments
    dplyr::mutate(
      pct_stand_n_trees_target = dplyr::case_when(
        as.numeric(clump_n_trees_grp)==largest_w_trees ~ pct_stand_n_trees_target
        , T ~ pct_stand_n_trees_target * (
          # pct remaining to scale to
          (1-max(ifelse(as.numeric(clump_n_trees_grp)==largest_w_trees,pct_stand_n_trees_target,0))) /
            # current pct remaining total allocated
            ifelse(sum(
              ifelse(as.numeric(clump_n_trees_grp)!=largest_w_trees,pct_stand_n_trees_target,0))==0,
              max(pct_stand_n_trees_target),sum(
                ifelse(as.numeric(clump_n_trees_grp)!=largest_w_trees,pct_stand_n_trees_target,0))
            )
        )
      )
    ) %>% 
    # add other targets
    dplyr::rename(pct_stand_n_trees_current = pct_stand_n_trees) %>% 
    dplyr::mutate(
      stand_trees_per_ac_current = clump_n_trees_grp_summary_dta$stand_trees_per_ac[1]
      , stand_trees_per_ac_target = dplyr::coalesce(get_tpa(target_ba, target_qmd),0)
      , trees_per_acre_current = stand_trees_per_ac_current*pct_stand_n_trees_current
      , trees_per_acre_target = stand_trees_per_ac_target*pct_stand_n_trees_target
      , clumps_per_acre_current = trees_per_acre_current/mean_clump_n_trees
      , clumps_per_acre_target = trees_per_acre_target/mean_clump_n_trees
      , stand_n_clumps_current = stand_n_clumps
      , stand_n_clumps_target = (clumps_per_acre_target*stand_area_ac) %>% round(0)
    ) %>% 
    dplyr::select(-c(tidyselect::starts_with("largest_w_trees"), stand_n_clumps))
  # ????  
  # target_data %>% glimpse()
  
  # get rid of nans
  # target_data$pct_stand_n_trees_target[is.na(target_data$pct_stand_n_trees_target)]=0
  # target_data$stand_n_clumps_target[is.na(target_data$stand_n_clumps_target)]=0
  
  # issue warning about targets
  if(min(target_data$pct_stand_n_trees_target == c(as.numeric(target_pcts), rep(0,6))[1:6])==0){
    warning(
      "proportion of trees in each clump size target `target_pcts` adjusted!!!"
      , "\nfrom : ", paste(round(target_pcts,2),collapse = ",")
      , "\nto : ", paste(round(target_data$pct_stand_n_trees_target,2),collapse = ",")
    )
  }
  # return
  return(target_data)
}


# Next two functions borrowed from https://github.com/metafor-ulaval/ALSroads/blob/main/R/line_tools.R


#---> Get heading of both ends of a line
st_ends_heading <- function(line){
  M <- sf::st_coordinates(line)
  i <- c(2, nrow(M) - 1)
  j <- c(1, -1)
  
  headings <- mapply(i, j, FUN = function(i, j) {
    Ax = M[i-j,1]
    Ay = M[i-j,2]
    Bx = M[i,1]
    By = M[i,2]
    atan2(Ay-By, Ax-Bx)*180/pi
  })
  names(headings) <- c("head", "tail")
  return(headings)
}


#---> extend the line on both ends
st_extend_line <- function(line, distance, end = "BOTH"){
  if (!(end %in% c("BOTH", "HEAD", "TAIL")) | length(end) != 1) stop("'end' must be 'BOTH', 'HEAD' or 'TAIL'")
  
  M <- sf::st_coordinates(line)[,-3]
  keep <- !(end == c("TAIL", "HEAD"))
  
  ends <- c(1, nrow(M))[keep]
  headings <- st_ends_heading(line)[keep] / 180 * pi
  distances <- if (length(distance) == 1) rep(distance, 2) else distance[1:2]
  
  M[ends, 1:2] <- M[ends, 1:2] + distances[keep] * c(cos(headings), sin(headings))
  newline <- sf::st_linestring(M)
  
  # If input is sfc_LINESTRING and not sfg_LINESTRING
  if (is.list(line)) newline <- sf::st_sfc(newline, crs = sf::st_crs(line))
  
  return(newline)
}


#---> pass an sf dataframe of points and return a line between the farthest points
st_points_to_line <- function(pts, line_ext=0) {
  if(max(class(pts) %in% c("sf"))!=1){
    stop("must provide an object of class `sf`")
  }
  # find farthest distance between points
  dist_temp = sf::st_distance(pts)
  
  # get the points
  f_pts_temp = 
    pts %>% 
    dplyr::ungroup() %>% 
    dplyr::slice(
      # get the farthest points from distance matrix
      which(dist_temp == max(dist_temp), arr.ind = TRUE)[1,]
    )
  
  # draw a line between the farthest two points
  f_line_temp = f_pts_temp %>% 
    # convert to linestring
    dplyr::ungroup() %>% 
    dplyr::summarise(n=dplyr::n()) %>%
    sf::st_cast("LINESTRING") %>%
    dplyr::ungroup() %>% 
    dplyr::select(-c(n))
  
  # and apply the line extension
  farthest_line = st_extend_line(f_line_temp, distance = line_ext)
  # return
  return(farthest_line)
}


#---> Function to calculate Euclidean distance between 2 points
st_euclidean_distance <- function(p1,p2) {
  return(sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2))
}


#---> return a line perpendicular to current line
### https://stackoverflow.com/questions/56771058/perpendicular-lines-at-regular-intervals-along-lines-with-multiple-nodes
# Function to calculate 2 points on a line perpendicular to another defined by 2 points p1,p2
# For point at interval, which can be a proportion of the segment length, or a constant
st_perp_line <- function(interval=0.5, my_line, proportion=TRUE) {
  # get end points of line
  p1 = my_line %>% sf::st_cast("POINT") %>% sf::st_coordinates() %>% .[1,]
  p2 = my_line %>% sf::st_cast("POINT") %>% sf::st_coordinates() %>% .[2,]
  # get length of line to return equal length line
  line_len = sf::st_length(my_line) %>% as.numeric() %>% `/`(2)      
  # get crs of line
  my_crs = sf::st_crs(my_line)
  
  # Calculate x and y distances
  x_len <- p2[1] - p1[1]
  y_len <- p2[2] - p1[2]
  
  # If proportion calculate reference point from tot_length
  if (proportion) {
    point <- c(p1[1]+x_len*interval,p1[2]+y_len*interval)
  }
  # Else use the constant value
  else {
    tot_len <- st_euclidean_distance(p1,p2)
    point <- c(p1[1]+x_len/tot_len*interval,p1[2]+y_len/tot_len*interval)
  }    
  
  # Calculate the x and y distances from reference point to point on line line_len distance away    
  ref_len <- st_euclidean_distance(point,p2)
  xn_len <- (line_len / ref_len) * (p2[1] - point[1])
  yn_len <- (line_len / ref_len) * (p2[2] - point[2])
  
  # fix for identical 
  if(identical(point,p2[1:2]) & x_len>y_len){ # this works for horizontal line
    xn_len <- line_len/2
    yn_len <- 0
  }else if(identical(point,p2[1:2]) & x_len<y_len){ # this works for vertical line
    xn_len <- 0
    yn_len <- line_len/2
  }
  
  # Invert the x and y lengths and add/subtract from the refrence point
  # ref_points <- rbind(point,c(point[1] + yn_len,point[2] - xn_len),c(point[1] - yn_len,point[2] + xn_len))
  ref_points <- rbind(c(point[1] + yn_len,point[2] - xn_len),c(point[1] - yn_len,point[2] + xn_len))
  
  # use the reference points to return a line
  return(
    ref_points %>% 
      dplyr::as_tibble() %>% 
      dplyr::rename_with(tolower) %>% 
      sf::st_as_sf(coords = c("x", "y"), crs = my_crs, remove = F) %>%
      dplyr::summarise(n=dplyr::n()) %>%
      sf::st_cast("LINESTRING") %>%
      dplyr::ungroup() %>% 
      dplyr::select(-c(n))
  )
}


#---> function to cut a single clump
cut_clump_fn <- function(
    c # clump id from need_cut_trees data 
    , need_cut_trees # data with clumps already defined returned by st_clump_points()
    , tgt = tgt # select one level of input passed to st_clump_points() ...
    # ... clump_labels = c("Individual","2-4 trees","5-9 trees","10-15 trees","16-25 trees",">25 trees")
){
  # use only the current clump
  curr_dta = need_cut_trees %>% 
    dplyr::filter(clump_id == c)
  
  if(nrow(curr_dta)<1){
    stop("cannot find data with the clump_id == `c` in `need_cut_trees` data")
  }
  
  # get distance used to create clumps in st_clump_points()
  dist_temp = curr_dta$tree_clump_dist_m[1]
  
  # do we even need to cut?
  reclump = st_clump_points(x = curr_dta, clump_dist_m = dist_temp)
  if(
    (
      reclump %>% 
      dplyr::pull(clump_id) %>% 
      unique() %>% 
      length()
    ) != 1
  ){
    stop("this is not a clump...send the `need_cut_trees` data through the st_clump_points function again")
  }
  # do we even need to cut?
  if(
    unique(reclump$clump_n_trees_grp) == tgt
  ){
    return(
      reclump %>% 
        sf::st_drop_geometry() %>% 
        dplyr::select(treeID, is_keep_tree)
    )
  }
  
  # create clump polygon
  clumps = get_clump_summary(curr_dta)
  
  # get the farthest line between points
  f_line_temp = st_points_to_line(curr_dta, line_ext = dist_temp)
  
  # get perpendicular lines in dataset which we can iterate over to make cuts
  perp_line_sf_temp = 
    # for every 1 m along line length, get a new perp line
    seq(
      from = 0
      , to = sf::st_length(f_line_temp) %>% as.numeric() %>% floor()
      , by = 1
    ) %>% 
    purrr::map(
      st_perp_line
      , my_line = f_line_temp
      , proportion = F
    ) %>% 
    dplyr::bind_rows() %>% 
    dplyr::mutate(line_n = dplyr::row_number())
  
  # find intersection of lines with the polygon and add length of intersection to perp line data
  perp_line_sf_temp = perp_line_sf_temp %>% 
    dplyr::inner_join(
      # intersect and calc len
      perp_line_sf_temp %>% 
        sf::st_intersection(
          clumps %>% 
            dplyr::ungroup() %>% 
            dplyr::select(clump_id) %>% 
            dplyr::filter(clump_id == c)
        ) %>% 
        dplyr::mutate(len_m = sf::st_length(geometry) %>% as.numeric()) %>% 
        sf::st_drop_geometry()
      , by = "line_n"
    )
  
  # make cuts at the points where there is the least overlap with the clump polygon
  # list of potential line cut + tree combinations
  # aggregate the mean length of the intersecting cut lines to the tree level
  # prioritize trees for removal that have smallest length
  cut_tree_lines_temp = 
    curr_dta %>% 
    sf::st_buffer(dist_temp/2) %>% 
    sf::st_join(perp_line_sf_temp) %>% 
    sf::st_drop_geometry() %>% 
    dplyr::group_by(treeID) %>% 
    dplyr::summarise(len_m = mean(len_m, na.rm = T)) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(len_m, treeID) %>% 
    dplyr::mutate(n = dplyr::row_number())
  ###############################################################
  # while
  ###############################################################
  # cut until the desired clump size is achieved based on these cut lines
  while_temp = 1
  i = 1
  while(while_temp==1) {
    # cut trees
    cut_trees_temp = cut_tree_lines_temp %>% 
      dplyr::slice(1:i) %>% 
      dplyr::distinct(treeID)
    
    # get the remaining trees not cut
    trees_remain_temp = curr_dta %>% 
      dplyr::anti_join(cut_trees_temp, by = "treeID")
    
    # ensure that there are trees
    if(nrow(trees_remain_temp)==0){
      if(best_desired_grps_n_temp==0){
        # there are no more possible cuts :/
        # ... so we're going to add trees until the desired clump size is reached
        # start with biggest tree and keep adding trees until the desired clump size is reached
        start_tree = curr_dta %>% 
          dplyr::arrange(desc(dbh_cm)) %>% 
          dplyr::filter(dplyr::row_number() == 1)
        
        increment_m = 0.5
        k = 1
        keep_trees = dplyr::tibble(treeID = character(0))
        while(nrow(keep_trees)==0){
          # create a polygon to intersect with tree points and keep trees until number of trees met
          poly_keep = start_tree %>% 
            sf::st_buffer( (dist_temp/2) + increment_m*k ) %>% 
            dplyr::mutate(dummy = 1) %>% 
            dplyr::select(dummy)
          # intersect and clump
          i_trees = curr_dta %>% 
            sf::st_intersection(poly_keep) %>% 
            st_clump_points(clump_dist_m = dist_temp)
          
          # Did we skip from multiple clumps smaller than the tgt to a clump larger than the tgt?
          # if so, select new start tree and start tree additions over.
          clump_temp_summary = i_trees %>% sf::st_drop_geometry() %>% dplyr::count(clump_n_trees_grp, clump_id)
          if(any(clump_temp_summary$clump_n_trees_grp>tgt)) {
            start_tree = curr_dta %>% 
              dplyr::arrange(desc(dbh_cm)) %>% 
              dplyr::filter(dplyr::row_number() == sample(1:nrow(curr_dta), 1))
            k=1
            next
          }
          
          # keep only the desired clump
          keep_trees = i_trees %>% 
            sf::st_drop_geometry() %>% 
            dplyr::filter(
              clump_id == 
                i_trees %>% 
                dplyr::filter(clump_n_trees_grp == tgt) %>% 
                dplyr::pull(clump_id) %>% 
                .[1] %>% 
                dplyr::coalesce("nope")
            ) %>% 
            dplyr::select(treeID)
          # increment
          k = k+1
        }
        # the remaining trees are cuts
        best_cuts = curr_dta %>% 
          sf::st_drop_geometry() %>% 
          dplyr::anti_join(keep_trees, by = "treeID") %>% 
          dplyr::select(treeID)
      }else if(
        best_desired_grps_n_temp>0 
      ){ 
        ###############################
        # add trees back in until gets worse
        ###############################
        j = 1
        while_add = 1
        while(while_add==1){
          # cut trees
          cut_trees_temp = best_cuts %>% 
            # add trees back in (i.e. remove them from the cut trees)
            dplyr::anti_join(
              cut_tree_lines_temp %>% 
                dplyr::slice(1:j) %>% 
                dplyr::distinct(treeID)
              , by = "treeID"
            )
          # get the remaining trees not cut
          trees_remain_temp = curr_dta %>% 
            dplyr::anti_join(cut_trees_temp, by = "treeID")
          
          # count the groups remaining after cuts
          desired_grps_n_temp = trees_remain_temp %>% 
            st_clump_points(clump_dist_m = dist_temp) %>% 
            sf::st_drop_geometry() %>% 
            # do we have group sizes we want?
            dplyr::filter(clump_n_trees_grp == tgt) %>% 
            nrow()
          
          if(desired_grps_n_temp>=best_desired_grps_n_temp){ # is this better than the best
            best_cuts = cut_trees_temp
            best_desired_grps_n_temp = desired_grps_n_temp
          }else{
            # stop it
            while_add = 0  
          }
          # increment
          j = j + 1
        } # while(while_add==1)
      } # best_desired_grps_n_temp>0 
      # done so stop the whole stop it
      while_temp = 0
    }else{ # if(nrow(trees_remain_temp)==0)
      # count the groups remaining after cuts
      desired_grps_n_temp = trees_remain_temp %>% 
        st_clump_points(clump_dist_m = dist_temp) %>% 
        sf::st_drop_geometry() %>% 
        # do we have group sizes we want?
        dplyr::filter(clump_n_trees_grp == tgt) %>% 
        nrow()
      ### store best cut list
      if(i==1){
        best_cuts = cut_trees_temp
        best_desired_grps_n_temp = desired_grps_n_temp
      }else if(desired_grps_n_temp>best_desired_grps_n_temp){ # is this better than the best
        best_cuts = cut_trees_temp
        best_desired_grps_n_temp = desired_grps_n_temp
      }else if(
        desired_grps_n_temp==best_desired_grps_n_temp
        & i!=nrow(cut_tree_lines_temp)
      ){
        best_cuts = best_cuts
        best_desired_grps_n_temp = best_desired_grps_n_temp
      }else if(
        best_desired_grps_n_temp>0 
        & desired_grps_n_temp<best_desired_grps_n_temp
      ){ # is this worse than the best which was successful?
        ###############################
        # add trees back in until gets worse
        ###############################
        j = 1
        while_add = 1
        while(while_add==1){
          # cut trees
          cut_trees_temp = best_cuts %>% 
            # add trees back in (i.e. remove them from the cut trees)
            dplyr::anti_join(
              cut_tree_lines_temp %>% 
                dplyr::slice(1:j) %>% 
                dplyr::distinct(treeID)
              , by = "treeID"
            )
          # get the remaining trees not cut
          trees_remain_temp = curr_dta %>% 
            dplyr::anti_join(cut_trees_temp, by = "treeID")
          
          # count the groups remaining after cuts
          desired_grps_n_temp = trees_remain_temp %>% 
            st_clump_points(clump_dist_m = dist_temp) %>% 
            sf::st_drop_geometry() %>% 
            # do we have group sizes we want?
            dplyr::filter(clump_n_trees_grp == tgt) %>% 
            nrow()
          
          if(desired_grps_n_temp>=best_desired_grps_n_temp){ # is this better than the best
            best_cuts = cut_trees_temp
            best_desired_grps_n_temp = desired_grps_n_temp
          }else{
            # stop it
            while_add = 0  
          }
          # increment
          j = j + 1
        } # while(while_add==1)
        # done so stop the whole stop it
        while_temp = 0
      }else if( i==nrow(cut_tree_lines_temp) ){ # is this the end?
        # stop it
        while_temp = 0
      }
      ### increment
      i = i+1
    } # else if(nrow(trees_remain_temp)==0)
  } # while(while_temp==1)
  
  # return it
  # return treelist with cut/keep
  # join to original data and pull
  d_temp = curr_dta %>%
    sf::st_drop_geometry() %>%
    dplyr::mutate(is_keep_tree = as.numeric(NA)) %>% 
    dplyr::select(-c(is_keep_tree)) %>% 
    dplyr::left_join(
      best_cuts %>% dplyr::mutate(is_keep_tree = 0)
      , by = dplyr::join_by("treeID")
    ) %>%
    dplyr::mutate(is_keep_tree = dplyr::coalesce(is_keep_tree, 1)) %>%
    dplyr::select(treeID, is_keep_tree)
  # returns treeID and keep tree flag data frame
  return(d_temp)
}


#---> function to if we met the clump cutting target, I guess.
get_cut_clump <- function(
    c # clump id from need_cut_trees data 
    , need_cut_trees # data with clumps already defined returned by st_clump_points()
    , tgt = tgt # select one level of input passed to st_clump_points() ...
    # ... clump_labels = c("Individual","2-4 trees","5-9 trees","10-15 trees","16-25 trees",">25 trees")
) {
  # what if we need to keep cutting because could not find a solution based on initial cut lines?
  cuts_first = cut_clump_fn(
    c = c
    , need_cut_trees = need_cut_trees
    , tgt = tgt
  )
  # did we get the target?
  new_clumps = 
    need_cut_trees %>% 
    dplyr::inner_join(
      cuts_first %>% dplyr::filter(is_keep_tree == 1) %>% dplyr::select(treeID)
      , by = "treeID"
    ) %>%
    st_clump_points()
  # new_clumps %>% sf::st_drop_geometry() %>% dplyr::count(clump_n_trees_grp)
  
  while(
    (
      new_clumps %>% dplyr::filter(clump_n_trees_grp > tgt) %>% nrow()
    ) > 0
  ){
    # redo the cut clump
    # redo the cut clump
    cuts_again = new_clumps %>% 
      dplyr::filter(clump_n_trees_grp > tgt) %>%
      dplyr::pull(clump_id) %>% 
      unique() %>% 
      purrr::map(cut_clump_fn, need_cut_trees = new_clumps, tgt = tgt) %>% 
      dplyr::bind_rows() %>% 
      dplyr::rename(updt = is_keep_tree)
    
    # update the original cut data
    cuts_first = cuts_first %>%
      dplyr::left_join(
        cuts_again
        , by = "treeID"
      ) %>% 
      dplyr::mutate(is_keep_tree = dplyr::coalesce(updt, is_keep_tree)) %>% 
      dplyr::select(-c(updt))
    
    # reset the new clumps
    new_clumps = 
      need_cut_trees %>% 
      dplyr::inner_join(
        cuts_first %>% dplyr::filter(is_keep_tree == 1) %>% dplyr::select(treeID)
        , by = "treeID"
      ) %>%
      st_clump_points()
    # new_clumps %>% sf::st_drop_geometry() %>% dplyr::count(clump_n_trees_grp)
  } # end while
  return(cuts_first) 
}


#---> get cut keep tree flag
get_keep_tree_flag <- function(
    x # x = sf point data
    , tgt # tgt = clump_n_trees_grp level defined in st_clump_points: ...
    ## ... clump_labels = c("Individual","2-4 trees","5-9 trees","10-15 trees","16-25 trees",">25 trees")
    , clump_dist_m = 6 # size (radius) of the epsilon neighborhood = maximum distance between points to add to clump
) {
  # MAKE CUTS INDEPENDENT OF CURRENT CLUMP GROUPING...
  # what if we pass a tree list with trees already cut? and thus, potentially different clump group sizes than is 
  # currently defined in clump_n_trees_grp?
  # 1) apply the clump grouping 
  # 2) only go through the cut alg for clumps that need cutting to the target
  # 3) append the trees in target clump or lower to the list at the end
  
  #############################
  # 1) apply the clump grouping 
  #############################
  new_clump_trees = x %>% 
    st_clump_points(clump_dist_m = clump_dist_m)
  
  # need cutting still
  need_cut_trees = new_clump_trees %>% 
    dplyr::filter(clump_n_trees_grp > tgt)
  
  #############################
  # 2) only go through the cut alg for clumps that need cutting to the target
  #############################
  # for each clump that still needs cutting, iterate over and return keep/cut flag
  compl_need_cut_trees = need_cut_trees %>% 
    dplyr::pull(clump_id) %>% 
    unique() %>% 
    purrr::map(get_cut_clump, need_cut_trees = need_cut_trees, tgt = tgt) %>%  # purrr:map fn
    dplyr::bind_rows()
  
  #############################
  # 4) append the trees in target clump or lower to the list at the end
  #############################
  compl_need_cut_trees = compl_need_cut_trees %>% 
    dplyr::bind_rows(
      new_clump_trees %>% 
        sf::st_drop_geometry() %>% 
        dplyr::filter(clump_n_trees_grp <= tgt) %>% 
        dplyr::mutate(is_keep_tree = 1) %>%
        dplyr::select(treeID, is_keep_tree)
    )
  
  # return
  return(
    # original data so that order is preserved
    x %>% 
      # add on trees that got cut with trees already good
      dplyr::inner_join(
        compl_need_cut_trees
        , by = "treeID"
      ) %>% 
      dplyr::pull(is_keep_tree)
  )
  
} # get_keep_tree_flag


################################################################################
######                          End                                       ######
################################################################################