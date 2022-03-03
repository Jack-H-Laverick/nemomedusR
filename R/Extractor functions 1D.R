
#### 1D extractors for taking a water column  ####

#' Extract a 1D water column from NEMO-MEDUSA
#' 
#' These functions read in target variables from NEMO-MEDUSA model outputs and return a water column.
#'
#' @details Each variable of interest in the netcdf file is imported, only reading at a single specified pixel for all depths.
#'  
#' Each model output file type contains a different set of variables, extracted by the relevant function variant:
#' 
#' | File type &nbsp; &nbsp; &nbsp; &nbsp;| Variables    |
#' |--------------|-------------|
#' | grid_T_      | Salinity, temperature, sea ice concentration.|
#' | grid_U_      | Zonal currents.|
#' | grid_V_      | Meridional currents.|
#' | grid_W_      | Vertical velocity, vertical eddy diffusivitiy.|
#' | ptrc_T_      | DIN, phytoplankton nitrogen content, phytoplankton chlorophyll, <br> zooplankton nitrogen content, detritus, silicate.|
#' 
#' Some function variants have different arguments:
#'
#' grid_W_ expects it's own `scheme_w` as depth levels are different between these and other files.
#'  
#' @md
#' @param filename The name of a netcdf file containing the title variables.
#' @param date The numeric string representing the date contained by the data file.
#' @param start an optional vector of indices to start subsetting at. See ncdf4 documentation.
#' @param count an optional vector of steps to subset along. See ncdf4 documentation.
#' @param ... soaks up unused function arguments passed by the wrapper functions handling file architecture.
#' @return A dataframe containing estimates of all variables at all depths for a single point in the model grid
#' @family NEMO-MEDUSA variable extractors
#' @name extractors_1D
NULL

#' @rdname extractors_1D
#' @export
get_grid_T_1D   <- function(filename, date, start, count) {
  
  nc_raw <- ncdf4::nc_open(filename)                                # Open up a netcdf file to see it's raw contents (var names)
  nc_saline <- ncdf4::ncvar_get(nc_raw, "vosaline", start = start, count = count)          # Extract an array of salinities
  nc_temp <- ncdf4::ncvar_get(nc_raw, "votemper", start = start, count = count)            # Extract an array of temperatures
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss
  
  all <- cbind(                                                                # Bind as columns
    Salinity = nc_saline,
    Temperature = nc_temp,
    Date = date,
    Layer = seq_along(nc_saline))
  return(all)
}

#' @rdname extractors_1D
#' @export
get_ptrc_T_1D <- function(filename, date, start, count) {
  
  nc_raw <- ncdf4::nc_open(filename)                                # Open up a netcdf file to see it's raw contents (var names)
  nc_DIN <- ncdf4::ncvar_get(nc_raw, "DIN", start = start, count = count)   # Extract an array of Dissolved inorganic nitrogen
  nc_CHD <- ncdf4::ncvar_get(nc_raw, "CHD", start = start, count = count)   # Extract an array of CHD
  nc_CHN <- ncdf4::ncvar_get(nc_raw, "CHN", start = start, count = count)
  nc_DET <- ncdf4::ncvar_get(nc_raw, "DET", start = start, count = count)
  nc_PHD <- ncdf4::ncvar_get(nc_raw, "PHD", start = start, count = count)
  nc_PHN <- ncdf4::ncvar_get(nc_raw, "PHN", start = start, count = count)
  nc_SIL <- ncdf4::ncvar_get(nc_raw, "SIL", start = start, count = count)
  nc_ZME <- ncdf4::ncvar_get(nc_raw, "ZME", start = start, count = count)
  nc_ZMI <- ncdf4::ncvar_get(nc_raw, "ZMI", start = start, count = count)
  
  
  ncdf4::nc_close(nc_raw)                                                           # You must close an open netcdf file when finished to avoid data loss
  
  all <- cbind(                                                                # Bind as columns
    DIN = nc_DIN,
    Diatom_chl = nc_CHD,
    Non_diatom_chl = nc_CHN,
    Diatom_phytoplankton_N = nc_PHD,
    Non_diatom_phytoplankton_N = nc_PHN,
    Detritus = nc_DET,
    Silicate = nc_SIL,
    Microzooplankton_N = nc_ZMI,
    Mesozooplankton_N = nc_ZME,
    Date = date,
    Layer = seq_along(nc_DIN))
  return(all)
}

#' @rdname extractors_1D
#' @export
get_grid_W_1D   <- function(filename, date, start, count) {
  
  nc_raw <- ncdf4::nc_open(filename)                                # Open up a netcdf file to see it's raw contents (var names)
  nc_vel <- ncdf4::ncvar_get(nc_raw, "vovecrtz", start = start, count = count)          # Extract an array of velocities
  nc_dif <- ncdf4::ncvar_get(nc_raw, "votkeavt", start = start, count = count)           # Extract an array of diffusivity
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss
  
  all <- cbind(                                                                # Bind as columns
    Vertical_velocity = nc_vel,
    Vertical_diffusivity = nc_dif,
    Date = date, 
    Layer = seq_along(nc_vel))
  return(all)
}

#' @rdname extractors_1D
#' @export
get_grid_V_1D <- function(filename, date, start, count) {
  
  nc_raw <- ncdf4::nc_open(filename)                                # Open up a netcdf file to see it's raw contents (var names)
  nc_merid <- ncdf4::ncvar_get(nc_raw, "vomecrty", start = start, count = count)          # Extract an array of meridinal currents
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss
  
  all <- cbind(                                                                # Bind as columns
    Meridional_velocity = nc_merid,
    Date = date, 
    Layer = seq_along(nc_merid))
  return(all)
}

#' @rdname extractors_1D
#' @export
get_grid_U_1D <- function(filename, date, start, count) {
  
  nc_raw <- ncdf4::nc_open(filename)                                # Open up a netcdf file to see it's raw contents (var names)
  nc_zonal <- ncdf4::ncvar_get(nc_raw, "vozocrtx", start = start, count = count)          # Extract an array of Zonal Current
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss
  
  all <- cbind(                                                                # Bind as columns
    Zonal_velocity = nc_zonal,
    Date = date, 
    Layer = seq_along(nc_zonal))
  return(all)
}
