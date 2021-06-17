
##** Functions to create slabR schemes for common types of model summary

#' Use slabR Schemes with netcdf Files
#'
#' These functions streamline the process of using slabR to summarise arrays contained within netcdf files 
#' 
#' ncdf4 uses a convention of 2 vectors, start and count, to read in a subset of a file. `sceheme_to_start()` and
#' `scheme_to_count` create the vectors to read in the minimum amount of data that contains the points needed in
#' the summary scheme. The `fill` argument allows you to add additional elements to the end of these vectors if necessary.
#' NEMO_MEUDSA model outputs have a useless 4th dimension of length 1. The `fill` default accounts for this.
#' 
#' Reading in a subset of data may change the indices needed in the summary scheme. If the X, Y and layer vectors used
#' to generate the scheme were from the whole netcdf file, then using it on a smaller subset will result in index out 
#' of bound errors. `scheme_reframe()` fixes this by rescaling the indices to start from the smallest in each dimension,
#' now matching the netcdf subset. 
#'
#' @param slabr_scheme a dataframe containing a slabR summary scheme.
#' @param fill a vector of values to add to the end of start or count objects to control extra dimensions.
#' @family NEMO-MEDUSA spatial tools
#' @name netcdf_scheme_helpers
NULL

#' @rdname netcdf_scheme_helpers
#' @export
scheme_to_start <- function(slabr_scheme = scheme, fill = 1){ 
  
  warning("If working with a subsetted netcdf file remember to use 'reframe_scheme()'.")
  string <- c(min(slabr_scheme$x), min(slabr_scheme$y), min(slabr_scheme$layer), fill)
  
  return(string)
}

#' @rdname netcdf_scheme_helpers
#' @export
scheme_to_count <- function(slabr_scheme = scheme, fill = -1){
  
  warning("If working with a subsetted netcdf file remember to use 'reframe_scheme()'.")
  string <- c(diff(range(slabr_scheme$x))+1, diff(range(slabr_scheme$y))+1, diff(range(slabr_scheme$layer))+1, fill) 
  
  return(string)
}

#' @rdname netcdf_scheme_helpers
#' @export
scheme_reframe <- function(slabr_scheme) {
  
  scheme <- mutate(scheme, 
                   x = (x+1) - min(x),
                   y = (y+1) - min(y),
                   layer = (layer+1) - min(layer))
}

#' Return a slabR Scheme for Linear Interpolation Between Array Slices
#'
#' Use this scheme to interpolate a layer at a new depth from an array.
#' 
#' The function takes the spatial information for a whole NEMO-MEDUSA file and determines which slice is immediately
#' above and below a target depth. These two slices are used to interpolate the new layer. The function calculates 
#' the weights needed to interpolate between the values using `calculate_proximity_weight`. If an sf polygon is provided,
#' only points which fall within the polygon will be interpolated. 
#'
#' @param space a list of latitudes, longitudes, and depths for NEMO-MEDUSA model output, as returned by `get_spatial()`
#' @param target_depth a depth (m) to interpolate a slice for.
#' @param sf an optional sf polygon to define the horizontal area of interest.
#' @return A dataframe containing a slabR scheme.
#' @family NEMO-MEDUSA spatial tools
#' @export
scheme_interp_slice <- function(space, target_depth, sf = NULL){
  
  bracket_depths <- c(max(space$nc_depth[space$nc_depth <= target_depth]), # Which depth in our layers is immediately below the target depth
                      min(space$nc_depth[space$nc_depth >= target_depth])) # and which is above
  
  if(bracket_depths[1] == bracket_depths[2]) {
    stop("Your target depth matches a provided depth value, there's no need to interpolate a new layer") }
  
  grid <- reshape2::melt(space$nc_lon) %>% 
    rename(x = "Var1", y = "Var2", longitude = "value") %>% 
    cbind(latitude = reshape2::melt(space$nc_lat)$value) %>%  # Convert space object into a dataframe of lat/lons 
    mutate(depth = list(bracket_depths))                      # Create a list column with the depths straddling our target depth.
  
  if(!is.null(sf)) {                                        # If an SF polygon was provided
    
    points <- st_as_sf(grid, coords = c("longitude", "latitude"), crs = 4326, remove = F) # Convert to SF points
    
    if(st_crs(points) != st_crs(sf)) points <- st_transform(points, st_crs(sf)) # Ensure points and sf share a projection
    
    grid <- st_join(points, sf) %>%                           # Test wihch grid points are in the polygon
      st_drop_geometry() %>%                                  # Drop heavy geometry column
      drop_na()                                               # Drop points outside of the polygon
  }  
  
  scheme <- grid %>%                                         
    group_by(y, x) %>%                                        # Group by each horizontal pixel
    mutate(group = cur_group_id(),                            # Create a grouping column for the summary scheme
           weight = map(depth, calculate_proximity_weight, target = target_depth), # Calculate the weights to interpolate between layers
           layer = list(which(space$nc_depth %in% depth[[1]]))) %>% # Which array layers match the depths we're interpolating?
    unnest(c(layer, weight, depth)) %>%                       # Open up the layer indices, weights, and the real world depths
    ungroup() 
  return(scheme)
}

#' Return a slabR Scheme for a StrathE2E Summary
#'
#' Use this scheme to create a depth-averaged shallow and deep grid for use with StrathE2E.
#' 
#' The function takes the spatial information for a whole NEMO-MEDUSA file and determines how to get a depth-averaged 
#' shallow and deep grid. The function uses `calculate_depth_share()` per pixel to determine how much of the water column
#' is represented by each slice, within either the shallow or deep grid. The depth limits of each layer and the 
#' bathymetry are used when calculating the weights for each point. If an sf polygon is provided, only points which fall 
#' within the polygon will be interpolated. 
#'
#' @param space a list of latitudes, longitudes, and depths for NEMO-MEDUSA model output, as returned by `get_spatial()`
#' @param bathymetry a dataframe containg the NEMO_MEDUSA bathymetry.
#' @param shallow a vector of a min and max depth (m) for the shallow layer to be averaged over.
#' @param deep a vector of a min and max depth (m) for the deep layer to be averaged over.
#' @param sf an optional sf polygon to define the horizontal area of interest.
#' @return A dataframe containing a slabR scheme.
#' @family NEMO-MEDUSA spatial tools
#' @export
scheme_strathE2E <- function(space, bathymetry, shallow, deep, sf = NULL){
  
  grid <- reshape2::melt(space$nc_lon) %>% 
    rename(x = "Var1", y = "Var2", longitude = "value") %>% 
    cbind(latitude = reshape2::melt(space$nc_lat)$value) %>%  # Convert space object into a dataframe of lat/lons 
    mutate(depth = list(space$nc_depth))                      # Create a list column with the depths straddling our target depth.
  
  if(!is.null(sf)) {                                        # If an SF polygon was provided
    
    points <- st_as_sf(grid, coords = c("longitude", "latitude"), crs = 4326, remove = F) # Convert to SF points
    
    if(st_crs(points) != st_crs(sf)) points <- st_transform(points, st_crs(sf)) # Ensure points and sf share a projection
    
    grid <- st_join(points, sf) %>%                           # Test wihch grid points are in the polygon
      st_drop_geometry() %>%                                  # Drop heavy geometry column
      drop_na()                                               # Drop points outside of the polygon
  }  
  
  grid <- left_join(grid, bathymetry) %>% 
    drop_na()
  
  # Shallow slab   
  shallow_scheme <- grid %>% 
    mutate(slab_layer = "1",
           Max_depth = ifelse(Bathymetry > shallow, shallow, Bathymetry),
           Min_depth = 0,
           weight = map2(Min_depth, Max_depth, .f = calculate_depth_share, depths = space$nc_depth), # Calculate the weights to interpolate between layers
           layer = list(which(space$nc_depth %in% depth[[1]]))) %>% # Which array layers match the depths we're interpolating?
    unnest(c(layer, weight, depth))                           # Open up the layer indices, weights, and the real world depths
  
  
  # Deep slab
  scheme <- grid %>% 
    mutate(slab_layer = "2",
           Max_depth = ifelse(Bathymetry > deep, deep, Bathymetry),
           Min_depth = shallow,
           weight = map2(Min_depth, Max_depth, .f = calculate_depth_share, depths = space$nc_depth), # Calculate the weights to interpolate between layers
           layer = list(which(space$nc_depth %in% depth[[1]]))) %>% # Which array layers match the depths we're interpolating?
    unnest(c(layer, weight, depth)) %>%                       # Open up the layer indices, weights, and the real world depths
    rbind(shallow_scheme) %>% 
    filter(weight != 0) %>%                                   # Drop points below the seafloor
    group_by(slab_layer, x, y) %>%                            # Group by each horizontal pixel
    mutate(group = cur_group_id()) %>%                        # Create a grouping column for the summary scheme
    ungroup()
  
  return(scheme)
} 

#' Return a slabR Scheme for a StrathE2E Summary
#'
#' Use this scheme to extract the values for the nearest depth column to a point location.
#' 
#' The function takes the spatial information for a whole NEMO-MEDUSA file and identifies the grid point closest
#' to a point location (the sf object). The returned sampling scheme contains every depth layer above the seafloor
#' with an equal weight. To retrieve the average depth column over a time period, use the "slabR" analysis with 
#' `NEMO_MEDUSA`. To Simply subset the depth column from model output, returning all values, use the "1D" analysis 
#' option.
#'
#' @param space a list of latitudes, longitudes, and depths for NEMO-MEDUSA model output, as returned by `get_spatial()`
#' @param bathymetry a dataframe containg the NEMO_MEDUSA bathymetry.
#' @param sf an sf point object defining the location to extract a representive water column for.
#' @return A dataframe containing a slabR scheme.
#' @family NEMO-MEDUSA spatial tools
#' @export

scheme_column <- function(space, bathymetry, sf = NULL){
  
  #space <- get_spatial(File, grid_W = F); bathymetry <- Bathymetry; sf <- location  
  
  grid <- reshape2::melt(space$nc_lon) %>% 
    rename(x = "Var1", y = "Var2", longitude = "value") %>% 
    cbind(latitude = reshape2::melt(space$nc_lat)$value) %>%  # Convert space object into a dataframe of lat/lons 
    mutate(depth = list(space$nc_depth))                      # Create a list column with the depths straddling our target depth.
  
  if(!is.null(sf)) {                                        # If an SF polygon was provided
    
    points <- st_as_sf(grid, coords = c("longitude", "latitude"), crs = 4326, remove = F) # Convert to SF points
    
    if(st_crs(points) != st_crs(sf)) points <- st_transform(points, st_crs(sf)) # Ensure points and sf share a projection
    
    grid <- st_join(sf, points, join = st_nearest_feature) %>% # Find nearest grid point
      st_drop_geometry() %>%                                  # Drop heavy geometry column
      drop_na()                                               # Drop points outside of the polygon
  }  
  
  grid <- left_join(grid, bathymetry) %>% 
    drop_na()
  
  scheme <- grid %>% 
    mutate(Max_depth = Bathymetry,
           Min_depth = 0,
           weight = map2(Min_depth, Max_depth, .f = calculate_depth_share, depths = space$nc_depth), # Calculate the weights to interpolate between layers
           layer = list(which(space$nc_depth %in% depth[[1]]))) %>% # Which array layers match the depths we're interpolating?
    unnest(c(layer, weight, depth)) %>%                       # Open up the layer indices, weights, and the real world depths
    filter(weight != 0) %>%                                   # Drop points below the seafloor
    group_by(layer) %>%                                       # Group by depth layer
    mutate(group = cur_group_id()) %>%                        # Create a grouping column for the summary scheme
    ungroup() %>%                                             
    mutate(weight = 1)                                        # Now we only have 1 point per group, turn off averaging in slabR
  
  return(scheme)
}

