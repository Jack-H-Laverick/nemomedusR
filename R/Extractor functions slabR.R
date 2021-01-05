
#### NEMO - MEDUSA file extractors    ####

#' Get Salinity, Temperature & Sea Ice Concentration
#'
#' This function reads in the title variables from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
#' Salinity, temperature & sea ice concentration can be found in files labelled "grid_T_".
#'
#' Each variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
#' `count3D`. The values are then passed to `array_w_mean` to calculate two average matrices, one a weighted vertical average
#' of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#'
#' @param path The path to the NEMO-MEDUSA model outputs.
#' @param file The name of a netcdf file containing the title variables.
#' @param scheme
#' @param ice_scheme
#' @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for title variables
#' in the shallow and deep zone. The dataframe contains the data for a single day.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_grid_T_slabR <- function(path, file, scheme, ice_scheme, space, ...) {

  print(stringr::str_glue("{file} Extracting Salinity, Temperature, and Sea Ice concentration"))
  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_saline <- ncdf4::ncvar_get(nc_raw, "vosaline", space$start3D, space$count3D)      # Extract an array of salinities
  nc_temp <- ncdf4::ncvar_get(nc_raw, "votemper", space$start3D, space$count3D)        # Extract an array of temperatures
  nc_ice <- ncdf4::ncvar_get(nc_raw, "soicecov", space$start3D[-3], space$count3D[-3]) # Extract a matrix of ice fractions
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                # Bind as columns
    Salinity = array_w_mean(nc_saline, scheme),                                # Summarise salinity according to scheme
    Temperature = array_w_mean(nc_temp, scheme),                               # Summarise temperature according to scheme
    Ice_conc  = c(nc_ice[ice_scheme],                                          # Subset ice instead of average as it's 2D
                  rep(NA, max(scheme$voxel)-length(ice_scheme))))              # Introduce NAs to fill deep layer summaries
  
    return(all)
}

#' Get Dissolved Inorganic Nitrogen & Phytoplankton Nitrogen
#'
#' This function reads in the title variables from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
#' DIN & phytoplankton Nitrogen content can be found in files labelled "ptrc_T_".
#'
#' Each variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
#' `count3D`. The values are then passed to `array_w_mean` to calculate two average matrices, one a weighted vertical average
#' of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#'
#' @param path The path to the NEMO-MEDUSA model outputs.
#' @param file The name of a netcdf file containing the title variables.
#' @param scheme
#' @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for title variables
#' in the shallow and deep zone. The dataframe contains the data for a single day.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_ptrc_T_slabR <- function(path, file, scheme, space, ...) {

  print(stringr::str_glue("{file} Extracting Dissolved Inorganic Nitrogen and Chlorophyll"))
  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_DIN <- ncdf4::ncvar_get(nc_raw, "DIN", space$start3D, space$count3D)      # Extract an array for the variable
  nc_DET <- ncdf4::ncvar_get(nc_raw, "DET", space$start3D, space$count3D)
  nc_PHD <- ncdf4::ncvar_get(nc_raw, "PHD", space$start3D, space$count3D)
  nc_PHN <- ncdf4::ncvar_get(nc_raw, "PHN", space$start3D, space$count3D)
  nc_Phyt <- nc_PHD + nc_PHN ; rm(nc_PHD, nc_PHN)
  ncdf4::nc_close(nc_raw)                                                          # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                    # Bind as matrix
    DIN = array_w_mean(nc_DIN, scheme),                                            # summarise DIN according to scheme
    Detritus = array_w_mean(nc_DET, scheme),                                       # summarise Detritus according to scheme
    Phytoplankton = array_w_mean(nc_Phyt, scheme))                                 # summarise Phytoplankton according to scheme
    return(all)
}

#' Get Ice Presence & Thickness, & Snow Thickness
#'
#' This function reads in the title variables from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
#' Ice presence & thickness, & snow thickness can be found in files labelled "icemod_".
#'
#' Each variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
#' `count3D`. The values are then passed to `array_w_mean` to calculate two average matrices, one a weighted vertical average
#' of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#'
#' @param path The path to the NEMO-MEDUSA model outputs.
#' @param file The name of a netcdf file containing the title variables.
#' @param scheme
#' @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for title variables
#' in the shallow and deep zone. The dataframe contains the data for a single day.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_icemod_slabR <- function(path, file, scheme, ice_scheme, space, ...) {

  print(stringr::str_glue("{file} Extracting Ice presence, and Ice and Snow thickness"))
  nc_raw <- ncdf4::nc_open(paste0(path, file))           # Open up a netcdf file to see it's raw contents (var names)
  nc_Ice <- ncdf4::ncvar_get(nc_raw, "ice_pres", space$start3D[-3], space$count3D[-3])# Extract a matrix of ice presence
  nc_Ithick <- ncdf4::ncvar_get(nc_raw, "iicethic", space$start3D[-3], space$count3D[-3]) # Extract ice thicknesses
  nc_Sthick <- ncdf4::ncvar_get(nc_raw, "isnowthi", space$start3D[-3], space$count3D[-3]) # Extract snow thicknesses
  ncdf4::nc_close(nc_raw)                                # You must close an open netcdf file when finished to avoid data loss

    all <- cbind(                                                     # Bind as a matrix
    Ice_pres = c(nc_Ice[ice_scheme],                                  # Subset ice instead of average as it's 2D
                 rep(NA, max(scheme$voxel)-length(ice_scheme))),      # Introduce NAs to fill deep layer summaries
    Ice_Thickness = c(nc_Ithick[ice_scheme],                          # Subset ice instead of average as it's 2D
                      rep(NA, max(scheme$voxel)-length(ice_scheme))), 
    Snow_Thickness = c(nc_Sthick[ice_scheme],                         
                       rep(NA, max(scheme$voxel)-length(ice_scheme))))
    return(all)
}

#' Get Vertical Velocity and Vertical Eddy Diffusivity
#'
#' This function reads in the title variables from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
#' Vertical velocity and vertical eddy diffusivitiy can be found in files labelled "W".
#'
#' Each variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
#' `count3D`. The values are then passed to `array_w_mean` to calculate two average matrices, one a weighted vertical average
#' of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#'
#' @param path The path to the NEMO-MEDUSA model outputs.
#' @param file The name of a netcdf file containing the title variables.
#' @param scheme_w
#' @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for title variables
#' in the shallow and deep zone. The dataframe contains the data for a single day.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_grid_W_slabR   <- function(path, file, scheme_w, space, ...) {

  print(stringr::str_glue("{file} Extracting Vertical water movements"))
  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_vel <- ncdf4::ncvar_get(nc_raw, "vovecrtz", space$start3DW, space$count3DW) # Extract an array for the variable
  nc_dif <- ncdf4::ncvar_get(nc_raw, "votkeavt", space$start3DW, space$count3DW)
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                # Bind as matrix
    Vertical_velocity = array_w_mean(nc_vel, scheme_w),                        # Summarise vertical velocity according to scheme
    Vertical_diffusivity = array_w_mean(nc_dif, scheme_w))                     # Summarise diffusivity according to scheme
  return(all)
}

#' Get Meridional Currents
#'
#' This function reads in the title variable from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
#' Meridional currents can be found in files labelled "grid_V_".
#'
#' Each variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
#' `count3D`. The values are then passed to `array_w_mean` to calculate two average matrices, one a weighted vertical average
#' of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#'
#' @param path The path to the NEMO-MEDUSA model outputs.
#' @param file The name of a netcdf file containing the title variable.
#' @param scheme
#' @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for the title variable
#' in the shallow and deep zone. The dataframe contains the data for a single day.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_grid_V_slabR <- function(path, file, scheme, space, ...) {

  print(stringr::str_glue("{file} Extracting Meridional currents"))
  nc_raw <- ncdf4::nc_open(paste0(path, file))                               # Open up a netcdf file to see it's raw contents (var names)
  nc_merid <- ncdf4::ncvar_get(nc_raw, "vomecrty", space$start3D, space$count3D)         # Pull meridinal currents
  ncdf4::nc_close(nc_raw)                                                    # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(Meridional = array_w_mean(nc_merid, scheme))                  # summarise meridional currents according to scheme
  return(all)
}

#' Get Zonal Currents
#'
#' This function reads in the title variable from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
#' Zonal currents can be found in files labelled "grid_U_".
#'
#' The variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
#' `count3D`. The values are then passed to `array_w_mean` to calculate two average matrices, one a weighted vertical average
#' of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#'
#' @param path The path to the NEMO-MEDUSA model outputs.
#' @param file The name of a netcdf file containing the title variables.
#' @param scheme
#' @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for the title variable
#' in the shallow and deep zone. The dataframe contains the data for a single day.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_grid_U_slabR <- function(path, file, scheme, space, ...) {

  print(stringr::str_glue("{file} Extracting Zonal currents"))
  nc_raw <- ncdf4::nc_open(paste0(path, file))                               # Open up a netcdf file to see it's raw contents (var names)
  nc_zonal <- ncdf4::ncvar_get(nc_raw, "vozocrtx", space$start3D, space$count3D) # Pull zonal current
  ncdf4::nc_close(nc_raw)                                                    # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(Zonal = array_w_mean(nc_zonal, scheme))                       # summarise zonal currents according to scheme
  return(all)
}
