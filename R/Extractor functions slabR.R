
#### slabR ####

#' Extract summaries from NEMO-MEDUSA arrays using slabR
#'
#' These functions read in target variables from NEMO-MEDUSA model outputs and return weighted averages 
#' according to a summary scheme.
#' 
#' @details Each variable of interest in a netcdf file is imported, only reading within an x/y window specified 
#' with `start` and `count`. The values are then passed to `array_w_mean()` to summarise according to a scheme. 
#' 
#' Each model output file type contains a different set of variables, extracted by the relevant function variant:
#' 
#' | File type &nbsp; &nbsp; &nbsp; &nbsp;| Variables    |
#' |--------------|-------------|
#' | grid_T_      | Salinity, temperature, sea ice concentration.|
#' | grid_U_      | Zonal currents.|
#' | grid_V_      | Meridional currents.|
#' | grid_W_      | Vertical velocity, vertical eddy diffusivitiy.|
#' | icemod_      | Ice presence, ice thickness, snow thickness.|
#' | ptrc_T_      | DIN, phytoplankton nitrogen content.|
#' 
#' Some function variants have different arguments:
#' 
#' grid_T_ and icemod_ functions accept an `ice_scheme`. Cryosphere data is contained in matrices not arrays.
#'
#' grid_W_ expects it's own `scheme_w` as depth levels are different between these and other files. it also expects it's
#' own `start_w` and `count_w`. 
#'  
#' @md
#' @param path the path to the NEMO-MEDUSA model outputs.
#' @param file the name of a netcdf file containing the title variables.
#' @param scheme a summary scheme as expected by `array_w_mean()`.
#' @param start an optional vector of indices to start subsetting at. See ncdf4 documentation.
#' @param count an optional vector of steps to subset along. See ncdf4 documentation.
#' @param scheme_w as above but for grid_W files.
#' @param start_w as above but for grid_W files.
#' @param count_w as above but for grid_W files.
#' @param ice_scheme a separate summary scheme for the 2D ice variables.
#' @param ... soaks up unused function arguments passed by the wrapper functions handling file architecture.
#' @return A matrix with a column of group averages per title variable for a single day.
#' @family NEMO-MEDUSA variable extractors
#' @name extractors_slabR
NULL

#' @rdname extractors_slabR
#' @export
get_grid_T_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ice_scheme, ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_saline <- ncdf4::ncvar_get(nc_raw, "vosaline", start, count)      # Extract an array of salinities
  nc_temp <- ncdf4::ncvar_get(nc_raw, "votemper", start, count)        # Extract an array of temperatures
  nc_ice <- ncdf4::ncvar_get(nc_raw, "soicecov", start[-3], count[-3]) # Extract a matrix of ice fractions
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                # Bind as columns
    Salinity = array_w_mean(nc_saline, scheme),                                # Summarise salinity according to scheme
    Temperature = array_w_mean(nc_temp, scheme),                               # Summarise temperature according to scheme
    Ice_conc  = c(nc_ice[ice_scheme],                                          # Subset ice instead of average as it's 2D
                  rep(NA, max(scheme$group)-length(ice_scheme))))              # Introduce NAs to fill deep layer summaries
  
    return(all)
}

#' @rdname extractors_slabR
#' @export
get_ptrc_T_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_DIN <- ncdf4::ncvar_get(nc_raw, "DIN", start, count)      # Extract an array for the variable
  nc_DET <- ncdf4::ncvar_get(nc_raw, "DET", start, count)
  nc_PHD <- ncdf4::ncvar_get(nc_raw, "PHD", start, count)
  nc_PHN <- ncdf4::ncvar_get(nc_raw, "PHN", start, count)
  nc_Phyt <- nc_PHD + nc_PHN
  ncdf4::nc_close(nc_raw)                                                          # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                    # Bind as matrix
    DIN = array_w_mean(nc_DIN, scheme),                                            # summarise DIN according to scheme
    Detritus = array_w_mean(nc_DET, scheme),                                       # summarise Detritus according to scheme
    Phytoplankton = array_w_mean(nc_Phyt, scheme))                                 # summarise Phytoplankton according to scheme
    return(all)
}

#' @rdname extractors_slabR
#' @export
get_grid_W_slabR   <- function(path, file, scheme_w, start_w = c(1,1,1,1), count_w = c(-1,-1,-1,-1), ...) {
  
  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_vel <- ncdf4::ncvar_get(nc_raw, "vovecrtz", start_w, count_w) # Extract an array for the variable
  nc_dif <- ncdf4::ncvar_get(nc_raw, "votkeavt", start_w, count_w)
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss
  
  all <- cbind(                                                                # Bind as matrix
    Vertical_velocity = array_w_mean(nc_vel, scheme_w),                        # Summarise vertical velocity according to scheme
    Vertical_diffusivity = array_w_mean(nc_dif, scheme_w))                     # Summarise diffusivity according to scheme
  return(all)
}

#' @rdname extractors_slabR
#' @export
get_grid_V_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ...) {
  
  nc_raw <- ncdf4::nc_open(paste0(path, file))                               # Open up a netcdf file to see it's raw contents (var names)
  nc_merid <- ncdf4::ncvar_get(nc_raw, "vomecrty", start, count)         # Pull meridinal currents
  ncdf4::nc_close(nc_raw)                                                    # You must close an open netcdf file when finished to avoid data loss
  
  all <- cbind(Meridional = array_w_mean(nc_merid, scheme))                  # summarise meridional currents according to scheme
  return(all)
}

#' @rdname extractors_slabR
#' @export
get_grid_U_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ...) {
  
  nc_raw <- ncdf4::nc_open(paste0(path, file))                               # Open up a netcdf file to see it's raw contents (var names)
  nc_zonal <- ncdf4::ncvar_get(nc_raw, "vozocrtx", start, count) # Pull zonal current
  ncdf4::nc_close(nc_raw)                                                    # You must close an open netcdf file when finished to avoid data loss
  
  all <- cbind(Zonal = array_w_mean(nc_zonal, scheme))                       # summarise zonal currents according to scheme
  return(all)
}

#' @rdname extractors_slabR
#' @export
get_icemod_slabR <- function(path, file, scheme, start = c(1,1,1,1), count = c(-1,-1,-1,-1), ice_scheme, ...) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))           # Open up a netcdf file to see it's raw contents (var names)
  nc_Ice <- ncdf4::ncvar_get(nc_raw, "ice_pres", start[-3], count[-3])# Extract a matrix of ice presence
  nc_Ithick <- ncdf4::ncvar_get(nc_raw, "iicethic", start[-3], count[-3]) # Extract ice thicknesses
  nc_Sthick <- ncdf4::ncvar_get(nc_raw, "isnowthi", start[-3], count[-3]) # Extract snow thicknesses
  ncdf4::nc_close(nc_raw)                                # You must close an open netcdf file when finished to avoid data loss

    all <- cbind(                                                     # Bind as a matrix
    Ice_pres = c(nc_Ice[ice_scheme],                                  # Subset ice instead of average as it's 2D
                 rep(NA, max(scheme$group)-length(ice_scheme))),      # Introduce NAs to fill deep layer summaries
    Ice_Thickness = c(nc_Ithick[ice_scheme],                          # Subset ice instead of average as it's 2D
                      rep(NA, max(scheme$group)-length(ice_scheme))), 
    Snow_Thickness = c(nc_Sthick[ice_scheme],                         
                       rep(NA, max(scheme$group)-length(ice_scheme))))
    return(all)
}
