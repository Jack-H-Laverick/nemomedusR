
#### Utilities to help with manipulating data ####

#' Calculate Water Layer Thicknesses Within a Vector
#'
#' This function calculates the thickness of the water layer around each point in a vector. These thicknesses
#' are needed to calculate weighted averages across a depth window.
#'
#' The function calculates the midpoints between values in a vector and an optional set of boundary depths.
#' The differences between these midpoints are calculated to return the thickness of a depth layer centered
#' on each point in the original vector. In the absence of boundary depths, the maximum and minimum depths are used.
#'
#' If a depth falls outside the target window it gets a thickness of 0. If all depths fall outside the window
#' the function returns an all 0 vector. If this is passed to `weighted.mean` the result will be NaN.
#'
#'
#' @param depths A numeric vector of increasing depths.
#' @param min_depth The shallowest depth in the depth window. (defaults to minimum depth in the vector)
#' @param max_depth The deepest depth in the depth window. (defaults to minimum depth in the vector)
#' @return A vector of water layer thicknesses to match the length of the original vector.
#' @family NEMO-MEDUSA spatial tools
#' @examples
#' # Get a vector of depths
#' depths <- seq(0, 100, by = 10)
#'
#' # Water layer thickness within the vector
#' calculate_depth_share(depths)
#'
#' # Water layer thickness using limits of a depth window
#' calculate_depth_share(depths, min_depth = 25, max_depth = 75)
#'
#' # Special case when the depth vector falls outside the target depth window
#' calculate_depth_share(depths, min_depth = 400, max_depth = 600)
#' @export
calculate_depth_share <- function(depths, min_depth = min(depths), max_depth = max(depths)) {

  contained <- depths <= max_depth & depths >= min_depth       # Which depths are within our window?

  if (all(!contained)) {

    depths <- rep(0, length(depths))                           # If none of the vector entries are within the depth window, return all 0s

  } else {

    weights <- diff(c(min_depth, RcppRoll::roll_mean(depths[contained], n = 2), max_depth)) # Calculate the midpoints and the difference between them (width either side of the original point)

    depths[contained] <- weights                               # Assign the weights to positions in the vector in our window
    depths[!contained] <- 0                                    # Positions outside the vector get 0 weight
  }
  return(depths)

}

#' Calculate the Weights for Linear Interpolation Between Two Depths
#'
#' This function calculates the distance between a shallower and deeper depth to a target depth. These values 
#' are then swapped so they can be used in a weighted average to linearly interpolate to the target depth.
#'
#' @param depths A numeric vector of a shallower and deeper depth.
#' @param target An intermediate depth we would like to interpolate to.
#' @return A vector of weights to match the vector of depths.
#' @family NEMO-MEDUSA spatial tools
#' @examples
#' calculate_proximity_weight(c(0, 100), target = 30)
#' 
#' calculate_proximity_weight(c(0, 100), target = 50)
#' 
#' calculate_proximity_weight(c(0, 100), target = 80)
#' @export
calculate_proximity_weight <- function(depths, target) rev(diff(c(depths[1], target, depths[2])))

#' Get Latitudes, Longitudes, & Depths From NEMO-MEDUSA Model Outputs 
#'
#' This function gets the latitudes, longitudes, and depths which define the spatial location of points in an array of NEMO-MEDUSA outputs.
#'
#' Each variable of interest in the netcdf file is imported, and then collected into a list. A T/F switch is included for
#' when the file is a grid_W type as the depth levels are different and stored differently.
#'
#' @param file The full name of a netcdf file.
#' @param grid_W Is the file a grid_W file, TRUE or FALSE. Needed as the depth variable is unique to this file type.
#' @return A list of three elements:
#' \itemize{
#'  \item{\emph{nc_lat -}}{ A matrix of latitudes which maps onto the first and second dimension of a NEMO-MEDUSA array.}
#'  \item{\emph{nc_lon -}}{ A matrix of longitudes which maps onto the first and second dimension of a NEMO-MEDUSA array.}
#'  \item{\emph{nc_depth -}}{ A vector of depths which match the third dimension of a NEMO-MEDUSA array.}
#'  }
#' @family NEMO-MEDUSA variable extractors
#' @export
get_spatial <- function(file, grid_W = F) {
  
  nc_raw <- ncdf4::nc_open(file)                       # Open up a netcdf file to see it's raw contents (var names)
  
  nc_lat <- ncdf4::ncvar_get(nc_raw, "nav_lat")        # Extract a matrix of all the latitudes
  nc_lon <- ncdf4::ncvar_get(nc_raw, "nav_lon")        # Extract a matrix of all the longitudes

  if(grid_W == F) {                                    # Extract a vector of depths
    nc_depth <- nc_raw$dim$deptht$vals} else {
    nc_depth <- nc_raw$dim$depthw$vals 
    }                 

  ncdf4::nc_close(nc_raw)                              # You must close an open netcdf file when finished to avoid data loss
  
  all <- list("nc_lat" = nc_lat, "nc_lon" = nc_lon, "nc_depth" = nc_depth)
  return(all)
}

#' Convert XY coordinates to a value index for a matrix
#'
#' This function converts a double index for a matrix (xy) into a single index for values. This allows vectorised subsetting 
#' of incomplete rows and columns.
#'
#' @param x vector of row numbers for target values.
#' @param y vector of column numbers for target values.
#' @param nrow The number of rows in the matrix
#' @return A vector of indices for values in a matrix.
#' @family NEMO-MEDUSA spatial tools
#' @export
xyindex_to_nindex <- function(x, y, nrow) {x + ((y-1)*nrow)}

#' Calculate the Domain Area per Grid Point
#'
#' This function takes an array of a variable, and an array of water thicknesses to perform a weighted average across depth. The depth
#' window to be averaged can be specified, so this function can be used to create both shallow and deep layers (or more for that matter).
#'
#' @param points A Simple Feature object og the grid points within the model domain.
#' @param area A Simple Feature object containing the model domain.
#' @return the `points` object is returned, but instead of points, the geometry column now contains polygons representing the area closest to each point. A column for the size of this area is also gained.
#' @family NEMO-MEDUSA spatial tools
#' @export
voronoi_grid <- function(points, area) {

  result <- purrr::map(1:nrow(area), ~{                            # For each polygon in area
    voronoi <- points %>%                                          # Take the grid points
      sf::st_geometry() %>%                                        # To get sfc from sf
      sf::st_union() %>%                                           # To get a sfc of MULTIPOINT type
      sf::st_voronoi(envelope = sf::st_geometry(area[.x,])) %>%    # Voronoi polygon for the area
      sf::st_collection_extract(type = "POLYGON") %>%              # A list of polygons
      sf::st_sf() %>%                                              # From list to sf object
      sf::st_join(points) %>%                                      # put names back
      sf::st_intersection(area[.x,]) %>%                           # Cut to shape of target area
      dplyr::mutate(Cell_area = units::drop_units(sf::st_area(.))) # Area of each polygon
  }) %>%
    dplyr::bind_rows() %>%                                         # Combine the results from each area
    sf::st_sf(geomc = .$geometry, crs = 4326)                      # Reinstate attributes of the geometry column

}
