#' Summarise Across Depths in a NEMO-MEDUSA Array
#'
#' This function takes an array of a variable, and an array of water thicknesses to perform a weighted average across depth. The depth
#' window to be averaged can be specified, so this function can be used to create both shallow and deep layers (or more for that matter).
#'
#' @param data An array containing estimates of a variable.
#' @param depth A TRUE FALSE vector indicating the depth layers to be extracted.
#' @param weights An array containing water thicknesses.
#' @return A matrix containing the weighted averages of the variable across the depth window of interest.
#' @family NEMO-MEDUSA spatial tools
#' @export
stratify  <- function(data, depth, weights) {
  
  # data <- nc_zonal ; depth <- Deep_mark ; weights <- dw                   # testing
  .Deprecated("array_w_mean with scheme_strathE2E")
  
  new <- data[,,depth] * weights[,,depth]                                      # Select slice of depths to average, multiply values by the weights
  empties <- emptyRcpp(new)                                                 # Find pixels with all depths shown by NA (locations of fake 0s)
  
  new2 <- rowSums(new, dims = 2, na.rm = TRUE)                                 # fast (C) sum the weighted values at a pixel
  denominator <- rowSums(weights[,,depth], dims = 2, na.rm = TRUE)             # fast (C) sum the weights
  weighted_mean <- new2/denominator                                            # Divide by the sum of the weights
  weighted_mean[empties] <- NA                                                 # Sum replaces an all NA dimension with 0, overwrite these by position
  return(weighted_mean)
}

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
#empty <- function(x){ .Deprecated("emptyRcpp") ; all(is.na(x))}

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
# get_weights <- function(top, bottom, bathymetry) {
#   #top <- 200                                                                    # shallowest depth of the slice
#   #bottom <- 1000                                                                # What's the bottom of the slice? for incorporating into a function
#   
#   .Deprecated("calculate_depth_share with apply(MARGIN(1,2)) or similar")
#   
#   weights <- array(NA, c(nrow(bathymetry), ncol(bathymetry),38))               # Initialise an array to hold weights
#   first <- weights[,,1]
#   first[] <- mean(Space$nc_depth[1:2]) - top %>% rep(times = length(bathymetry))# Water above the first midpoint minus the depth at the top of the slice
#   marks <- bathymetry > mean(Space$nc_depth[1:2])                              # Which cells contain model outputs deeper than the seafloor?
#   first[!marks] <- bathymetry[!marks] - top                                    # Replace these with the depth to sea floor
#   weights[,,1] <- first
#   
#   weights[,,38] <- bathymetry - mean(Space$nc_depth[37:38])                    # The remaining water column thickness is the sea floor - the deepest midpoint.
#   
#   for (i in 2:37) {
#     #i <- 23
#     last_midpoint <- mean(Space$nc_depth[(i-1):i])                               # Find the mid depth to the layer above
#     next_midpoint <- mean(Space$nc_depth[i:(i+1)])                               # Find the mid depth to the layer below
#     
#     if(top > last_midpoint) above <- top else above <- last_midpoint             # If the top of the slice is deeper than the previous midpoint, use the top of the slice
#     if(bottom < next_midpoint) below <- bottom else below <- next_midpoint       # If the next midpoint is deeper than the bottom of the slice, use the bottom of the slice
#     
#     weights[,,i] <- below - above %>% rep(times = length(bathymetry))                         # Calculate layer thickness and repeat to fill the array
#     
#     marks <- bathymetry > below                                                  # Is the seafloor deeper than the bottom of the layer?
#     weights[,,i][!marks] <- bathymetry[!marks] - above                           # If not, replace these with the depth to sea floor - the top of the water layer
#     
#   }                                                          # Roll through each matrix and calculate the water thickness using the next depth, bottom of the slice, or bathymetry, whichever is smaller
#   no_weight <- weights[] <= 0; weights[no_weight] <- NA                        # Finally if a weight is <= 0 get NA
#   
#   return(weights)
# }

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
# get_weights.W <- function(top, bottom, bathymetry) {
#   
#   .Deprecated("calculate_depth_share with apply(MARGIN(1,2)) or similar")
#   
#   weights <- array(NA, c(nrow(bathymetry),ncol(bathymetry),39))                  # Initialise an array to hold weights
#   first <- weights[,,1]
#   first[] <- mean(DepthsW[1:2]) - top %>% rep(times = length(bathymetry))        # Water above the first midpoint minus the depth at the top of the slice
#   marks <- bathymetry > mean(DepthsW[1:2])                                       # Which cells contain model outputs deeper than the seafloor?
#   first[!marks] <- bathymetry[!marks] - top                                      # Replace these with the depth to sea floor
#   weights[,,1] <- first
#   
#   weights[,,39] <- bathymetry - mean(DepthsW[38:39])                             # The remaining water column thickness is the sea floor - the deepest midpoint.
#   
#   for (i in 2:38) {
#     #i <- 23
#     last_midpoint <- mean(DepthsW[(i-1):i])                                      # Find the mid depth to the layer above
#     next_midpoint <- mean(DepthsW[i:(i+1)])                                      # Find the mid depth to the layer below
#     
#     if(top > last_midpoint) above <- top else above <- last_midpoint             # If the top of the slice is deeper than the previous midpoint, use the top of the slice
#     if(bottom < next_midpoint) below <- bottom else below <- next_midpoint       # If the next midpoint is deeper than the bottom of the slice, use the bottom of the slice
#     
#     weights[,,i] <- below - above %>% rep(times = length(bathymetry))            # Calculate layer thickness and repeat to fill the array
#     
#     marks <- bathymetry > below                                                  # Is the seafloor deeper than the bottom of the layer?
#     weights[,,i][!marks] <- bathymetry[!marks] - above                           # If not, replace these with the depth to sea floor - the top of the water layer
#     
#   }                                                          # Roll through each matrix and calculate the water thickness using the next depth, bottom of the slice, or bathymetry, whichever is smaller
#   no_weight <- weights[] <= 0; weights[no_weight] <- NA                        # Finally if a weight is <= 0 get NA
#   
#   return(weights)
# }
