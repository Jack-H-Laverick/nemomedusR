
#### NEMO - MEDUSA file extractors    ####

# Get Salinity, Temperature & Sea Ice Concentration
#
# This function reads in the title variables from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
# Salinity, temperature & sea ice concentration can be found in files labelled "grid_T_".
#
# Each variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
# `count3D`. The values are then passed to `stratify` to calculate two average matrices, one a weighted vertical average
# of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#
# @param path The path to the NEMO-MEDUSA model outputs.
# @param file The name of a netcdf file containing the title variables.
# @param start3D
# @param count3D
# @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for title variables
# in the shallow and deep zone. The dataframe contains the data for a single day.
# @family NEMO-MEDUSA variable extractors


#' @title DEPRECATED. Extract a StrathE2E summary from NEMO-MEDUSA 
#' 
#' These functions read in target variables from NEMO-MEDUSA model outputs and return a summarised slab.
#' @name extractors_SE2E
NULL

#' @rdname extractors_SE2E
#' @export
get_grid_T_StrathE2E <- function(path, file, space) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_saline <- ncdf4::ncvar_get(nc_raw, "vosaline", space$start3D, space$count3D)          # Extract an array of salinities
  nc_temp <- ncdf4::ncvar_get(nc_raw, "votemper", space$start3D, space$count3D)            # Extract an array of temperatures
  nc_ice <- ncdf4::ncvar_get(nc_raw, "soicecov", space$start3D[-3], space$count3D[-3])                               # Extract a matrix of ice fractions
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                # Bind as columns
    Ice_conc = c(as.numeric(nc_ice),                                           # shallow ice
                 rep(NA, length(nc_ice))),                                     # NAs for ice in deep layer
    Salinity = c(as.numeric(stratify(nc_saline, space$shallow, space$s.weights)), # Collapse shallow salinity into 2D and convert to long format
                as.numeric(stratify(nc_saline, space$deep, space$d.weights))),    # And deep salinity as one column
    Temperature = c(as.numeric(stratify(nc_temp, space$shallow, space$s.weights)),# Collapse shallow temperatures into 2D and convert to long format
                    as.numeric(stratify(nc_temp, space$deep, space$d.weights))))  # And deep temperature as one column
    return(all)
}

# Get Dissolved Inorganic Nitrogen & Phytoplankton Nitrogen
#
# This function reads in the title variables from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
# DIN & phytoplankton Nitrogen content can be found in files labelled "ptrc_T_".
#
# Each variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
# `count3D`. The values are then passed to `stratify` to calculate two average matrices, one a weighted vertical average
# of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#
# @param path The path to the NEMO-MEDUSA model outputs.
# @param file The name of a netcdf file containing the title variables.
# @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for title variables
# in the shallow and deep zone. The dataframe contains the data for a single day.
# @family NEMO-MEDUSA variable extractors
# @export

#' @rdname extractors_SE2E
#' @export
get_ptrc_T_StrathE2E <- function(path, file, space) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_DIN <- ncdf4::ncvar_get(nc_raw, "DIN", space$start3D, space$count3D)      # Extract an array for the variable
  nc_DET <- ncdf4::ncvar_get(nc_raw, "DET", space$start3D, space$count3D)
  nc_PHD <- ncdf4::ncvar_get(nc_raw, "PHD", space$start3D, space$count3D)
  nc_PHN <- ncdf4::ncvar_get(nc_raw, "PHN", space$start3D, space$count3D)
  nc_Phyt <- nc_PHD + nc_PHN ; rm(nc_PHD, nc_PHN)
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                # Bind as matrix
    DIN = c(as.numeric(stratify(nc_DIN, space$shallow, space$s.weights)),      # Collapse shallow DIN into 2D and convert to long format
            as.numeric(stratify(nc_DIN, space$deep, space$d.weights))),        # and deep as one column
    Detritus = c(as.numeric(stratify(nc_DET, space$shallow, space$s.weights)), # Collapse shallow Detritus into 2D and convert to long format
      as.numeric(stratify(nc_DET, space$deep, space$d.weights))),              # and deep as one column
    Phytoplankton = c(as.numeric(stratify(nc_Phyt, space$shallow, space$s.weights)), # Collapse shallow chlorophyll into 2D and convert to long format
                    as.numeric(stratify(nc_Phyt, space$deep, space$d.weights))))   # and deep as one column
    return(all)
}

# Get Ice Pesence & Thickness, & Snow Thickness
#
# This function reads in the title variables from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
# Ice presence & thickness, & snow thickness can be found in files labelled "icemod_".
#
# Each variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
# `count3D`. The values are then passed to `stratify` to calculate two average matrices, one a weighted vertical average
# of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#
# @param path The path to the NEMO-MEDUSA model outputs.
# @param file The name of a netcdf file containing the title variables.
# @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for title variables
# in the shallow and deep zone. The dataframe contains the data for a single day.
# @family NEMO-MEDUSA variable extractors
# @export

#' @rdname extractors_SE2E
#' @export
get_icemod_StrathE2E <- function(path, file, space) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))           # Open up a netcdf file to see it's raw contents (var names)
  nc_Ice <- ncdf4::ncvar_get(nc_raw, "ice_pres", space$start3D[-3], space$count3D[-3])# Extract a matrix of ice presence
  nc_Ithick <- ncdf4::ncvar_get(nc_raw, "iicethic", space$start3D[-3], space$count3D[-3]) # Extract ice thicknesses
  nc_Sthick <- ncdf4::ncvar_get(nc_raw, "isnowthi", space$start3D[-3], space$count3D[-3]) # Extract snow thicknesses
  ncdf4::nc_close(nc_raw)                                # You must close an open netcdf file when finished to avoid data loss

  length <- length(nc_Ice)                               # How many NAs will we need to fill?

  all <- cbind(                                          # Bind as a matrix
    Ice_pres = c(as.numeric(nc_Ice),                     # Flatten ice values
                 rep(NA, length)),                       # And repeat NAs for the deep layer to match the length of other get_functions
    Ice_Thickness = c(as.numeric(nc_Ithick),
                      rep(NA, length)),
    Snow_Thickness = c(as.numeric(nc_Sthick),
                       rep(NA, length)))
    return(all)
}

# Get Vertical Velocity and Vertical Eddy Diffusivity
#
# This function reads in the title variables from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
# Vertical velocity and vertical eddy diffusivitiy can be found in files labelled "W".
#
# Each variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
# `count3D`. The values are then passed to `stratify` to calculate two average matrices, one a weighted vertical average
# of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#
# @param path The path to the NEMO-MEDUSA model outputs.
# @param file The name of a netcdf file containing the title variables.
# @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for title variables
# in the shallow and deep zone. The dataframe contains the data for a single day.
# @family NEMO-MEDUSA variable extractors
# @export

#' @rdname extractors_SE2E
#' @export
get_grid_W_StrathE2E   <- function(path, file, space) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_vel <- ncdf4::ncvar_get(nc_raw, "vovecrtz", space$start3DW, space$count3DW) # Extract an array for the variable
  nc_dif <- ncdf4::ncvar_get(nc_raw, "votkeavt", space$start3DW, space$count3DW)
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                # Bind as matrix
    Vertical_velocity = c(as.numeric(stratify(nc_vel, space$shallow_W, space$s.weights_W)), # Collapse shallow DIN into 2D and convert to long format
                          as.numeric(stratify(nc_vel, space$deep_W, space$d.weights_W))),   # And deep as one column
    Vertical_diffusivity = c(as.numeric(stratify(nc_dif, space$shallow_W, space$s.weights_W)), # Collapse shallow chlorophyll into 2D and convert to long format
                             as.numeric(stratify(nc_dif, space$deep_W, space$d.weights_W))))   # And deep as one column
  return(all)
}

# Get Meridional Currents
#
# This function reads in the title variable from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
# Meridional currents can be found in files labelled "grid_V_".
#
# Each variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
# `count3D`. The values are then passed to `stratify` to calculate two average matrices, one a weighted vertical average
# of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#
# @param path The path to the NEMO-MEDUSA model outputs.
# @param file The name of a netcdf file containing the title variable.
# @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for the title variable
# in the shallow and deep zone. The dataframe contains the data for a single day.
# @family NEMO-MEDUSA variable extractors
# @export

#' @rdname extractors_SE2E
#' @export
get_grid_V_StrathE2E <- function(path, file, space) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                               # Open up a netcdf file to see it's raw contents (var names)
  nc_merid <- ncdf4::ncvar_get(nc_raw, "vomecrty", space$start3D, space$count3D)         # Pull meridinal currents
  ncdf4::nc_close(nc_raw)                                                    # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                              # Create 1 column matrix
    Meridional = c(as.numeric(stratify(nc_merid, space$shallow, space$s.weights)), # Collapse shallow meridional currents into 2D and convert to long format
                   as.numeric(stratify(nc_merid, space$deep, space$d.weights))))   # And deep in the same column
  return(all)
}

# Get Zonal Currents
#
# This function reads in the title variable from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
# Zonal currents can be found in files labelled "grid_U_".
#
# The variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
# `count3D`. The values are then passed to `stratify` to calculate two average matrices, one a weighted vertical average
# of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#
# @param path The path to the NEMO-MEDUSA model outputs.
# @param file The name of a netcdf file containing the title variables.
# @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for the title variable
# in the shallow and deep zone. The dataframe contains the data for a single day.
# @family NEMO-MEDUSA variable extractors
# @export

#' @rdname extractors_SE2E
#' @export
get_grid_U_StrathE2E <- function(path, file, space) {

  nc_raw <- ncdf4::nc_open(paste0(path, file))                               # Open up a netcdf file to see it's raw contents (var names)
  nc_zonal <- ncdf4::ncvar_get(nc_raw, "vozocrtx", space$start3D, space$count3D)         # Pull zonal current
  ncdf4::nc_close(nc_raw)                                                    # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                               # Create a matrix
    Zonal = c(as.numeric(stratify(nc_zonal, space$shallow, space$s.weights)), # Collapse shallow zonal currents into 2D and convert to long format
              as.numeric(stratify(nc_zonal, space$deep, space$d.weights))))   # And for deep as one column
  return(all)
}
