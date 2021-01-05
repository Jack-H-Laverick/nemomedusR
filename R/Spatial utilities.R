
#### Utilities to help with manipulating data ####

#' Check Whether a Vector Contains Data
#'
#' This function provides a check for whether a vector of locations contains any data. This is useful
#' for checking whether a lat/lon pixel doesn't report any data for the deep zone. This is expected if the
#' sea floor at those coordinates is shallower than the maximum depth we set for the shallow zone.
#'
#' @param x A vector of values to check.
#' @return TRUE (if only NAs) or FALSE (if any entry is not NA).
#' @family NEMO-MEDUSA spatial tools
#' @export
empty <- function(x) all(is.na(x))

#' Calculate Water Layer Thicknesses Within an Array
#'
#' This function calculates the thickness of the water layer around each point in a NEMO-MEDUSA array. These thicknesses
#' are needed to calculate weighted averages across a depth window.
#'
#' The function starts by calculating the thickness of the depth window represented by the shallowest layer in the array.
#' This is the depth to halfway between the shallowest and the next layer, subtracting any depth to the top of the window, i.e.
#' If the top of the window is the sea surface, the first array layer is at 5m, and the second is at 10m, the water thickness about
#' the first point is 7.5m (mean(5, 10) - 0). The function then checks whether the depth of this layer of points is shallower than the depth
#' to the seafloor for any of the points. For these points which exist below the seafloor, thicknesses are replaced by the depth to the seafloor,
#' minus the depth to the top of the depth window.
#'
#' The function then populates the thicknesses of the mid-layers in a for-loop, working increasingly deeper. The thickness about a point is
#' calculated as the difference between the depth halfway to the point above, and the depth halfway to the point below. Checks and replacements
#' are performed, as explained above, to see whether any points are now beyond the limit of the bathymetry, or approaching the `top` or
#' `bottom` of the depth window. If they are, the claculations are performed again using the new max and min depths.
#'
#' The thicknesses for the deepest layer of points are calculated as the depth to the seafloor minus the mid-depth of the deepest two layers.
#'
#' If a weight is <= 0, it is replaced with NA as this indicates land.
#'
#' @param top The shallowest depth in the depth window.
#' @param bottom The deepest depth in the depth window.
#' @return An array of water layer thicknesses to match the dimensions of a NEMO-MEDUSA array.
#' @family NEMO-MEDUSA spatial tools
#' @export
get_weights <- function(top, bottom, bathymetry) {
  #top <- 200                                                                    # shallowest depth of the slice
  #bottom <- 1000                                                                # What's the bottom of the slice? for incorporating into a function

  weights <- array(NA, c(nrow(bathymetry), ncol(bathymetry),38))               # Initialise an array to hold weights
  first <- weights[,,1]
  first[] <- mean(Space$nc_depth[1:2]) - top %>% rep(times = length(bathymetry))# Water above the first midpoint minus the depth at the top of the slice
  marks <- bathymetry > mean(Space$nc_depth[1:2])                              # Which cells contain model outputs deeper than the seafloor?
  first[!marks] <- bathymetry[!marks] - top                                    # Replace these with the depth to sea floor
  weights[,,1] <- first

  weights[,,38] <- bathymetry - mean(Space$nc_depth[37:38])                    # The remaining water column thickness is the sea floor - the deepest midpoint.

  for (i in 2:37) {
    #i <- 23
    last_midpoint <- mean(Space$nc_depth[(i-1):i])                               # Find the mid depth to the layer above
    next_midpoint <- mean(Space$nc_depth[i:(i+1)])                               # Find the mid depth to the layer below

    if(top > last_midpoint) above <- top else above <- last_midpoint             # If the top of the slice is deeper than the previous midpoint, use the top of the slice
    if(bottom < next_midpoint) below <- bottom else below <- next_midpoint       # If the next midpoint is deeper than the bottom of the slice, use the bottom of the slice

    weights[,,i] <- below - above %>% rep(times = length(bathymetry))                         # Calculate layer thickness and repeat to fill the array

    marks <- bathymetry > below                                                  # Is the seafloor deeper than the bottom of the layer?
    weights[,,i][!marks] <- bathymetry[!marks] - above                           # If not, replace these with the depth to sea floor - the top of the water layer

  }                                                          # Roll through each matrix and calculate the water thickness using the next depth, bottom of the slice, or bathymetry, whichever is smaller
  no_weight <- weights[] <= 0; weights[no_weight] <- NA                        # Finally if a weight is <= 0 get NA

  return(weights)
}

#' Calculate Water Layer Thicknesses Within an Array (Velocities)
#'
#' This function is a variant of `get_weights` which applies to water velocities. Water velocities from NEMO-MEDUSA are on a different
#' vector of depths. This function calculates the thickness of the water layer around each point in a NEMO-MEDUSA array. These thicknesses
#' are needed to calculate weighted averages across a depth window.
#'
#' The function starts by calculating the thickness of the depth window represented by the shallowest layer in the array.
#' This is the depth to halfway between the shallowest and the next layer, subtracting any depth to the top of the window, i.e.
#' If the top of the window is the sea surface, the first array layer is at 5m, and the second is at 10m, the water thickness about
#' the first point is 7.5m (mean(5, 10) - 0). The function then checks whether the depth of this layer of points is shallower than the depth
#' to the seafloor for any of the points. For these points which exist below the seafloor, thicknesses are replaced by the depth to the seafloor,
#' minus the depth to the top of the depth window.
#'
#' The function then populates the thicknesses of the mid-layers in a for-loop, working increasingly deeper. The thickness about a point is
#' calculated as the difference between the depth halfway to the point above, and the depth halfway to the point below. Checks and replacements
#' are performed, as explained above, to see whether any points are now beyond the limit of the bathymetry, or approaching the `top` or
#' `bottom` of the depth window. If they are, the claculations are performed again using the new max and min depths.
#'
#' The thicknesses for the deepest layer of points are calculated as the depth to the seafloor minus the mid-depth of the deepest two layers.
#'
#' If a weight is <= 0, it is replaced with NA as this indicates land.
#'
#' @param top The shallowest depth in the depth window.
#' @param bottom The deepest depth in the depth window.
#' @return An array of water layer thicknesses to match the dimensions of a NEMO-MEDUSA array.
#' @family NEMO-MEDUSA spatial tools
#' @export
get_weights.W <- function(top, bottom, bathymetry) {

  weights <- array(NA, c(nrow(bathymetry),ncol(bathymetry),39))                  # Initialise an array to hold weights
  first <- weights[,,1]
  first[] <- mean(DepthsW[1:2]) - top %>% rep(times = length(bathymetry))        # Water above the first midpoint minus the depth at the top of the slice
  marks <- bathymetry > mean(DepthsW[1:2])                                       # Which cells contain model outputs deeper than the seafloor?
  first[!marks] <- bathymetry[!marks] - top                                      # Replace these with the depth to sea floor
  weights[,,1] <- first

  weights[,,39] <- bathymetry - mean(DepthsW[38:39])                             # The remaining water column thickness is the sea floor - the deepest midpoint.

  for (i in 2:38) {
    #i <- 23
    last_midpoint <- mean(DepthsW[(i-1):i])                                      # Find the mid depth to the layer above
    next_midpoint <- mean(DepthsW[i:(i+1)])                                      # Find the mid depth to the layer below

    if(top > last_midpoint) above <- top else above <- last_midpoint             # If the top of the slice is deeper than the previous midpoint, use the top of the slice
    if(bottom < next_midpoint) below <- bottom else below <- next_midpoint       # If the next midpoint is deeper than the bottom of the slice, use the bottom of the slice

    weights[,,i] <- below - above %>% rep(times = length(bathymetry))            # Calculate layer thickness and repeat to fill the array

    marks <- bathymetry > below                                                  # Is the seafloor deeper than the bottom of the layer?
    weights[,,i][!marks] <- bathymetry[!marks] - above                           # If not, replace these with the depth to sea floor - the top of the water layer

  }                                                          # Roll through each matrix and calculate the water thickness using the next depth, bottom of the slice, or bathymetry, whichever is smaller
  no_weight <- weights[] <= 0; weights[no_weight] <- NA                        # Finally if a weight is <= 0 get NA

  return(weights)
}

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

#' Get Latitudes, Longitudes, & Depths
#'
#' This function gets the latitudes, longitudes, and depths which define the spatial location of points in an array of NEMO-MEDUSA outputs.
#'
#' Each variable of interest in the netcdf file is imported, and then collected into a list.
#'
#' @param file The full name of a netcdf file.
#' @return A list of three elements:
#' \itemize{
#'  \item{\emph{nc_lat -}}{ A matrix of latitudes which maps onto the first and second dimension of a NEMO-MEDUSA array.}
#'  \item{\emph{nc_lon -}}{ A matrix of longitudes which maps onto the first and second dimension of a NEMO-MEDUSA array.}
#'  \item{\emph{nc_depth -}}{ A vector of depths which match the third dimension of a NEMO-MEDUSA array.}
#'  }
#' @family NEMO-MEDUSA variable extractors
#' @export
get_spatial <- function(file) {
  nc_raw <- nc_open(file)                              # Open up a netcdf file to see it's raw contents (var names)
  nc_lat <- ncvar_get(nc_raw, "nav_lat")               # Extract a matrix of all the latitudes
  nc_lon <- ncvar_get(nc_raw, "nav_lon")               # Extract a matrix of all the longitudes
  nc_depth <- ncvar_get(nc_raw, "deptht")              # Extract a matrix of depths
  nc_close(nc_raw)                                     # You must close an open netcdf file when finished to avoid data loss
  all <- list("nc_lat" = nc_lat, "nc_lon" = nc_lon, "nc_depth" = nc_depth)
  return(all)
}

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

#' Get Latitudes, Longitudes, & Depths
#'
#' This function gets the latitudes, longitudes, and depths which define the spatial location of points in an array of NEMO-MEDUSA outputs.
#'
#' Each variable of interest in the netcdf file is imported, and then collected into a list.
#'
#' @param file The full name of a netcdf file.
#' @return A list of three elements:
#' \itemize{
#'  \item{\emph{nc_lat -}}{ A matrix of latitudes which maps onto the first and second dimension of a NEMO-MEDUSA array.}
#'  \item{\emph{nc_lon -}}{ A matrix of longitudes which maps onto the first and second dimension of a NEMO-MEDUSA array.}
#'  \item{\emph{nc_depth -}}{ A vector of depths which match the third dimension of a NEMO-MEDUSA array.}
#'  }
#' @family NEMO-MEDUSA variable extractors
#' @export
get_spatial <- function(file) {
  nc_raw <- nc_open(file)                              # Open up a netcdf file to see it's raw contents (var names)
  nc_lat <- ncvar_get(nc_raw, "nav_lat")               # Extract a matrix of all the latitudes
  nc_lon <- ncvar_get(nc_raw, "nav_lon")               # Extract a matrix of all the longitudes
  nc_depth <- ncvar_get(nc_raw, "deptht")              # Extract a matrix of depths
  nc_close(nc_raw)                                     # You must close an open netcdf file when finished to avoid data loss
  all <- list("nc_lat" = nc_lat, "nc_lon" = nc_lon, "nc_depth" = nc_depth)
  return(all)
}

#' Convert a U-V velocity field to speed and direction in degrees
#'
#' This function takes a vector of u and v velocities and calculates the direction and speed of the combined movement.
#'
#'This function was lifted from the `Rsenal` package, where it was originally used to calculate wind speeds. All I've done
#'is built a wrapper which accounts for different conventions when describing wind and flow directions.
#'
#' @param u A vector of Zonal currents (from West to East).
#' @param v A vector of Meridional currents (from South to North).
#' @return a dataframe of two columns is returned. Speed contains the composite speed of both velocities on the same scale.
#' Direction is the resolved direction of the flow in degrees, 0 heads north, 90 East, 180 South, 270 West.
#' @family NEMO-MEDUSA spatial tools
#' @export
vectors_2_direction <- function (u, v) {
  u <- -u                                        # This function was built to use wind direction
  v <- -v                                        # Winds  are "opposite", people care about where wind comes from, not where it goes

  # Lovingly lifted from the "Rsenal" package

  degrees <- function(radians) 180 * radians/pi
  mathdegs <- degrees(atan2(v, u))
  wdcalc <- ifelse(mathdegs > 0, mathdegs, mathdegs + 360)
  uvDirection <- ifelse(wdcalc < 270, 270 - wdcalc, 270 - wdcalc + 360)
  uvSpeed <- sqrt(u^2 + v^2)
  return(cbind(uvDirection, uvSpeed))
  }

## used for Light and air temperature data which goes into NM, this data uses a different grid and has time stored differently

#' Get Indices to Use When Clipping netcdf Files at Import
#'
#' This function works out how much of a netcdf file to read, to capture the data between a given Lat-Lon window.
#'
#' The function reads in a vector for both latitudes and longitudes, and tests whether each entry is within the specified
#' window. The max and min position in these vectors where the condition == TRUE are taken to define the ends of the window
#' to import. The vectors of latitudes and longitudes between these limits are kept, so they can be added to the variables
#' of interest during extraction.
#'
#' @param file The full name of a netcdf file containing a longitude and latitude dimension.
#' @param w Degrees West to read from.
#' @param e Degrees East to read to.
#' @param s Degrees South to read from.
#' @param n Degrees North to read to.
#' @return A list of three elements:
#' \itemize{
#'  \item{\emph{Lats -}}{ A vector of latitudes from `s` to `n`.}
#'  \item{\emph{Lons -}}{ A vector of longitudes from `w` to `e`.}
#'  \item{\emph{Limits -}}{ A dataframe containing the index to start reading from (Lon_start, Lat_start)
#'  and the length of the vector to read (Lon_count, Lat_count.}
#'  }
#' @family NEMO-MEDUSA spatial tools
#' @export
Window <- function(file, w, e, s, n) {

  #file <- examples[1,]$File ; w = 0 ; e = 180 ; s = 0 ; n = 90

  raw <- ncdf4::nc_open(file)
  lon <- raw$dim$longitude$vals %>% between(w, e)
  W <- min(which(lon == TRUE))
  E <- max(which(lon == TRUE))

  lat <- raw$dim$latitude$vals %>% between(s, n)
  S <- min(which(lat == TRUE))
  N <- max(which(lat == TRUE))

  lons <- raw$dim$longitude$vals[W:E]
  lats <- raw$dim$latitude$vals[S:N]

  Limits <- data.frame("Lon_start" = W, "Lon_count" = E - W + 1, "Lat_start" = S, "Lat_count" = N - S + 1)

  Limits <- list(Lats = lats, Lons = lons, Limits = Limits)
  return(Limits)
}

#' Pull Coordinates from a Simple Feature Geometry Column
#'
#' This function takes an SF object, and adds two columns containing the coordinates in the geometry column.
#'
#' @param data An SF (Simple Feature) object.
#' @return The same object, now with two columns containing the coordinates in the geometry column.
#' @family NEMO-MEDUSA spatial tools
#' @export
sfc_as_cols <- function(x, names = c("x","y")) {
  stopifnot(inherits(x,"sf") && inherits(sf::st_geometry(x),"sfc_POINT"))
  ret <- sf::st_coordinates(x)
  ret <- as.data.frame(ret)
  stopifnot(length(names) == ncol(ret))
  x <- x[ , !names(x) %in% names]
  ret <- setNames(ret,names)
  dplyr::bind_cols(x,ret)
}

#' Reproject from Latitude and Longitude to Project CRS
#'
#' This function takes a dataframe containing a latitude and longitude column, and replaces them with an X and Y column of coordinates
#' in a new CRS.
#'
#'The function converts a dataframe into an SF object and reprojects into a new CRS. Two coordinate columns are extracted from the
#'geometry column using `sfc_as_cols`, before the geometry column is dropped.
#'
#' @param data A dataframe containing Longitude and Latitude.
#' @param crs The new Coordinate Reference System  to project to.
#' @return A dataframe, now with an x and y column specifying the coordinates for points in the projects Coordiante Reference System.
#' @family NEMO-MEDUSA spatial tools
#' @export
reproj <- function(data, crs) {

  data %>%
    sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% # Specify original projection (crs)
    sf::st_transform(crs = crs) %>%                                   # Transform to crs specified in region file
    sfc_as_cols() %>%                                                 # Extract geometry column for geom_segment to work
    sf::st_set_geometry(NULL)                                         # Chuck geometry column
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
