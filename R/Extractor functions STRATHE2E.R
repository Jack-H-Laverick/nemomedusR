
#### NEMO - MEDUSA file extractors    ####

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

  new <- data[,,depth] * weights[,,depth]                                      # Select slice of depths to average, multiply values by the weights
  empties <- emptyRcpp(new)                                                 # Find pixels with all depths shown by NA (locations of fake 0s)

  new2 <- rowSums(new, dims = 2, na.rm = TRUE)                                 # fast (C) sum the weighted values at a pixel
  denominator <- rowSums(weights[,,depth], dims = 2, na.rm = TRUE)             # fast (C) sum the weights
  weighted_mean <- new2/denominator                                            # Divide by the sum of the weights
  weighted_mean[empties] <- NA                                                 # Sum replaces an all NA dimension with 0, overwrite these by position
  return(weighted_mean)
}

#' Get Salinity, Temperature & Sea Ice Concentration
#'
#' This function reads in the title variables from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
#' Salinity, temperature & sea ice concentration can be found in files labelled "grid_T_".
#'
#' Each variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
#' `count3D`. The values are then passed to `stratify` to calculate two average matrices, one a weighted vertical average
#' of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#'
#' @param path The path to the NEMO-MEDUSA model outputs.
#' @param file The name of a netcdf file containing the title variables.
#' @param start3D
#' @param count3D
#' @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for title variables
#' in the shallow and deep zone. The dataframe contains the data for a single day.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_grid_T_StrathE2E <- function(path, file, space) {

  print(stringr::str_glue("{file} Extracting Salinity, Temperature, and Sea Ice concentration"))
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

#' Get Dissolved Inorganic Nitrogen & Chlorophyll a
#'
#' This function reads in the title variables from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
#' DIN & chlorphyll a can be found in files labelled "ptrc_T_".
#'
#' Each variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
#' `count3D`. The values are then passed to `stratify` to calculate two average matrices, one a weighted vertical average
#' of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#'
#' @param path The path to the NEMO-MEDUSA model outputs.
#' @param file The name of a netcdf file containing the title variables.
#' @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for title variables
#' in the shallow and deep zone. The dataframe contains the data for a single day.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_ptrc_T_StrathE2E <- function(path, file, space) {

  print(stringr::str_glue("{file} Extracting Dissolved Inorganic Nitrogen and Chlorophyll"))
  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_DIN <- ncdf4::ncvar_get(nc_raw, "DIN", space$start3D, space$count3D)      # Extract an array for the variable
  nc_DET <- ncdf4::ncvar_get(nc_raw, "DET", space$start3D, space$count3D)
  nc_CHD <- ncdf4::ncvar_get(nc_raw, "CHD", space$start3D, space$count3D)
  nc_CHN <- ncdf4::ncvar_get(nc_raw, "CHN", space$start3D, space$count3D)
  nc_Chl <- nc_CHD + nc_CHN ; rm(nc_CHD, nc_CHN)
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                                # Bind as matrix
    DIN = c(as.numeric(stratify(nc_DIN, space$shallow, space$s.weights)),      # Collapse shallow DIN into 2D and convert to long format
            as.numeric(stratify(nc_DIN, space$deep, space$d.weights))),        # and deep as one column
    Detritus = c(as.numeric(stratify(nc_DET, space$shallow, space$s.weights)), # Collapse shallow Detritus into 2D and convert to long format
      as.numeric(stratify(nc_DET, space$deep, space$d.weights))),              # and deep as one column
    Chlorophyll = c(as.numeric(stratify(nc_Chl, space$shallow, space$s.weights)), # Collapse shallow chlorophyll into 2D and convert to long format
                    as.numeric(stratify(nc_Chl, space$deep, space$d.weights))))   # and deep as one column
    return(all)
}

#' Get Ice Pesence & Thickness, & Snow Thickness
#'
#' This function reads in the title variables from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
#' Ice presence & thickness, & snow thickness can be found in files labelled "icemod_".
#'
#' Each variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
#' `count3D`. The values are then passed to `stratify` to calculate two average matrices, one a weighted vertical average
#' of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#'
#' @param path The path to the NEMO-MEDUSA model outputs.
#' @param file The name of a netcdf file containing the title variables.
#' @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for title variables
#' in the shallow and deep zone. The dataframe contains the data for a single day.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_icemod_StrathE2E <- function(path, file, space) {

  print(stringr::str_glue("{file} Extracting Ice presence, and Ice and Snow thickness"))
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

#' Get Vertical Velocity and Vertical Eddy Diffusivity
#'
#' This function reads in the title variables from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
#' Vertical velocity and vertical eddy diffusivitiy can be found in files labelled "W".
#'
#' Each variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
#' `count3D`. The values are then passed to `stratify` to calculate two average matrices, one a weighted vertical average
#' of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#'
#' @param path The path to the NEMO-MEDUSA model outputs.
#' @param file The name of a netcdf file containing the title variables.
#' @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for title variables
#' in the shallow and deep zone. The dataframe contains the data for a single day.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_grid_W_StrathE2E   <- function(path, file, space) {

  print(stringr::str_glue("{file} Extracting Vertical water movements"))
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

#' Get Meridional Currents
#'
#' This function reads in the title variable from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
#' Meridional currents can be found in files labelled "grid_V_".
#'
#' Each variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
#' `count3D`. The values are then passed to `stratify` to calculate two average matrices, one a weighted vertical average
#' of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#'
#' @param path The path to the NEMO-MEDUSA model outputs.
#' @param file The name of a netcdf file containing the title variable.
#' @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for the title variable
#' in the shallow and deep zone. The dataframe contains the data for a single day.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_grid_V_StrathE2E <- function(path, file, space) {

  print(stringr::str_glue("{file} Extracting Meridional currents"))
  nc_raw <- ncdf4::nc_open(paste0(path, file))                               # Open up a netcdf file to see it's raw contents (var names)
  nc_merid <- ncdf4::ncvar_get(nc_raw, "vomecrty", space$start3D, space$count3D)         # Pull meridinal currents
  ncdf4::nc_close(nc_raw)                                                    # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                              # Create 1 column matrix
    Meridional = c(as.numeric(stratify(nc_merid, space$shallow, space$s.weights)), # Collapse shallow meridional currents into 2D and convert to long format
                   as.numeric(stratify(nc_merid, space$deep, space$d.weights))))   # And deep in the same column
  return(all)
}

#' Get Zonal Currents
#'
#' This function reads in the title variable from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
#' Zonal currents can be found in files labelled "grid_U_".
#'
#' The variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
#' `count3D`. The values are then passed to `stratify` to calculate two average matrices, one a weighted vertical average
#' of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#'
#' @param path The path to the NEMO-MEDUSA model outputs.
#' @param file The name of a netcdf file containing the title variables.
#' @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for the title variable
#' in the shallow and deep zone. The dataframe contains the data for a single day.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_grid_U_StrathE2E <- function(path, file, space) {

  print(stringr::str_glue("{file} Extracting Zonal currents"))
  nc_raw <- ncdf4::nc_open(paste0(path, file))                               # Open up a netcdf file to see it's raw contents (var names)
  nc_zonal <- ncdf4::ncvar_get(nc_raw, "vozocrtx", space$start3D, space$count3D)         # Pull zonal current
  ncdf4::nc_close(nc_raw)                                                    # You must close an open netcdf file when finished to avoid data loss

  all <- cbind(                                                               # Create a matrix
    Zonal = c(as.numeric(stratify(nc_zonal, space$shallow, space$s.weights)), # Collapse shallow zonal currents into 2D and convert to long format
              as.numeric(stratify(nc_zonal, space$deep, space$d.weights))))   # And for deep as one column
  return(all)
}

#' Get Rivers
#'
#' This function reads river NEMO-MEDUSA river runoff files and reshapes for StrathE2E.
#'
#' The river runoff file for a single year is imported. Points are checked to see if they fall in the model domain.
#' The points which do are summed by months.
#'
#' NOTE! The river run off files provided to MiMeMo have estimates in open ocean. Runoff is channeled to cells in
#' the vicinity of river mouths so that single grid cells don't get the full amount. Of course river run off really
#' comes from the coast, so you may need to modify your domain polygon.
#'
#' @param File The name (with path) of a netcdf file containing river runoff data.
#' @param Year The year the file contains data for.
#' @param domain The polygon used to capture point estimates relevant to the model domain.
#' @return A dataframe containing the total river runoff into the model domain by month with a year column.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_rivers <- function(File, Year, domain) {

  data <- raster::raster(File, varname = "nav_lat") %>%
    raster::as.data.frame(xy = TRUE) %>%                                                              # Get a dataframe of xy positions and latitudes
    dplyr::left_join(raster::as.data.frame(raster::raster(File, varname = "nav_lon"), xy = TRUE)) %>% # Bind to the same for longitudes
    dplyr::left_join(raster::as.data.frame(raster::brick(File, varname = "sorunoff"), xy = TRUE)) %>% # Bind to the river run off values for all months
    sf::st_as_sf(coords = c("nav_lon", "nav_lat"), crs = 4326) %>%                                    # Convert to sf
    dplyr::select(-c(x,y)) %>%                                                                        # Drop cell positions used for binding
    tidyr::drop_na() %>%                                                                              # Drop empty pixels
    sf::st_transform(sf::st_crs(domain)) %>%                                                          # Transform points to the same crs as domain polygon
    sf::st_join(domain) %>%                                                                           # Check which points are in the domain
    sf::st_drop_geometry() %>%                                                                        # Drop sf formatting
    tidyr::drop_na() %>%                                                                              # Drop points outside the domain
    colSums()                                                                                         # Total all river runoff in a month

  names(data) <- c(1:12)                                                                              # Label each column by month

  result <- data.frame(Year = Year, Month = 1:12, Runoff = data)                                      # Create a dataframe of runoff by month and add year column

  return(result)
}

#' Condense Daily netcdf Files into a Month by Type
#'
#' This function takes the metadata for multiple netcdf files and creates a spatial summary for the month.
#' \strong{For one file type only}.
#'
#' The function bridges the step between extracting data from a netcdf file, and creating an average dataframe
#' for a whole month of NEMO-MEDUSA model outputs.
#'
#' Different file types require a different get function. This function takes a collection of netcdf files of
#' the same time from the same month, and passes them to the correct `get_*` for data extraction. The results
#' are passed back to this function, before the estimates from different days at the same lat-lon-depth
#' combination are averaged to get a single number for the month.
#'
#' Creating intermediate monthly objects allows the terabytes of data to be reduced without encountering memory issues.
#' Also, working on independent monthly packets of data means we can parallelise any data processing for speed.
#'
#' @param data A dataframe containing the metadata of multiple netcdf files which share a type.
#' @param ... Additional arguments passed to the relevant `get_*` function.
#' @return The function returns a dataframe containing the monthly average shalllow and deep spatial grids for
#' variables of interest.
#' @family NEMO-MEDUSA variable extractors
#' @export
type_in_month_old <- function(data, ...) {

  Type <- data[1,3]                                 # Pull type from the file

get <- eval(parse(text = paste0("get_", Type)))     # Change the extracting function based on file contents

  Month.type <- purrr::map2(data$Path, data$File, get, ...) %>% # Summarise the contents of each file
    abind::abind(along = 3) %>%                     # Stack matrices behind each other
    rowMeans(dims = 2)                              # Quickly take the mean for all variables in the month

  return(Month.type)
}

#' Condense Daily netcdf Files into a Monthly Summary
#'
#' This function takes the metadata for multiple netcdf files. It then creates a common spatial summary for all
#' the variables of interest for a month.
#'
#' The function takes the metadata for a collection of files which contain data from the same month. The files
#' are split into data packets which share the same file type, before being passed to `type_in_month` to be summarised.
#' `type_in_month` reduces the NEMO-MEDUSA model outputs from large arrays to effectively two matrices. The summaries for
#' each file type are returned to this function and get bound into a single dataframe. Points outside the project window are
#' removed before saving the dataframe for the month in "./Objects/Months/".
#'
#' Creating intermediate monthly objects allows the terabytes of data to be reduced without encountering memory issues.
#' Also, working on independent monthly packets of data means we can parallelise any data processing for speed.
#'
#' @param data A dataframe containing the metadata of multiple netcdf files from a common month.
#' @param crop
#' @param ... Additional arguments to be passed to get_* functions.
#' @return The function returns a dataframe containing the monthly average shalllow and deep spatial grids for
#' \strong{all} the variables of interest in NEMO-MEDUSA outputs.
#' @family NEMO-MEDUSA variable extractors
#' @export
whole_month_old <- function(data, crop, grid, ...) {

  Month <- data[1,5] ; Year <- data[1,4]                                    # Pull date

  Month <- split(data, f = list(data$Type)) %>%                             # Split out the files for this month by type, so they can be averaged together
    purrr::map(type_in_month, ...) %>%                                      # Pull a whole month of data from a single file type
    do.call(cbind, .) %>%                                                   # Join together all the data packets
    cbind(grid) %>%                                                         # Add coordinates and depth labels
    filter(Shore_dist > 0) %>%                                              # Drop points on land
    mutate(Year = as.integer(Year),                                         # Add time
           Month =  as.integer(Month)) %>%
    dplyr::right_join(crop) %>%    # b)                                     # Cut out rows outside of plotting window
    saveRDS(., file = paste("./Objects/Months/NM", Month, Year, "rds", sep = "."))    # save out a data object for one whole month
}

#' Get Surface Irradiance & Air Temperature
#'
#' This function reads either title variable from a NEMO-MEDUSA model \strong{DRIVERS} file, and returns a monthly time series.
#'
#' The appropriate variable in the netcdf file is imported according to the `Type` parameter, only reading within an
#' x/y window specified in `Space`. The function then drops points outside the model domain before constructing the monthly
#' time series.
#'
#' Unlike other NEMO-MEDUSA related get_* functions, this function constructs the monthly time series directly. A wrapper function
#' to handle time isn't required, as each netcdf file for driving data contains 360 day steps for a single year. `stratify` also
#' isn't called, as these variables have no depth dimension.
#'
#' @param File The location of the netcdf file containing NEMO-MEDUSA driving data.
#' @param Type The variable contained in the netcdf file either "T150" (air temperature) or "SWF" (surface irradiance).
#' @param Year The year the necdf file contains data for.
#' @return A dataframe containing a monthly time series within a year of either average air temperature or surface
#' irradiance. Air temperature is also split by shore zone.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_air <- function(File, Type, Year) {

  #File <- Airtemp_files$File[1] ; Type <- Airtemp_files$Type[1] ; Year <- Airtemp_files$Year[1] # test
  if(Type == "SWF") months <- Light_months                                   # Get the file holding the months
  if(Type == "T150") months <- Airtemp_months                                # For the time steps of this data

  nc_raw <- ncdf4::nc_open(File)                                             # Open up a netcdf file to see it's raw contents (var names)
  nc_var <- ncdf4::ncvar_get(nc_raw, Type, c(Space$Limits$Lon_start, Space$Limits$Lat_start, 1),  # Extract the variable of interest
                      c(Space$Limits$Lon_count, Space$Limits$Lat_count, -1)) # cropped to window, with all time steps
  ncdf4::nc_close(nc_raw)                                                           # You must close an open netcdf file when finished to avoid data loss

  Data <- as.data.frame.table(nc_var, responseName = "Measured") %>%         # Reshape array as dataframe
    dplyr::rename(Longitude = Var1, Latitude = Var2, Time_step = Var3) %>%   # Name the columns
    dplyr::mutate(Longitude = rep(rep(Space$Lons,                            # Replace the factor levels with dimension values
                               times = length(unique(Latitude))), times = length(unique(Time_step))),
           Latitude = rep(rep(Space$Lats,
                              each = length(unique(Longitude))), times = length(unique(Time_step))),
           Time_step = rep(1:length(unique(Time_step)),
                           each = length(unique(Latitude)) * length(unique(Longitude)))) %>%
    dplyr::right_join(domains_mask) %>%                                      # Crop to domain
    dplyr::left_join(months) %>%                                             # Assign a month to each time step
    dplyr::mutate(Year = Year,                                               # Attach Year
           Type = Type)                                                      # Attach variable name

  if(Type == "SWF") Data <- dplyr::group_by(Data, Month, Year, Type)         # We don't need to bother accounting for shore in light data
  if(Type == "T150") Data <- dplyr::group_by(Data, Month, Year, Type, Shore) # We care about shore for temperature, retain interesting columns

  Summary <- dplyr::summarise(Data, Measured = stats::weighted.mean(Measured, Cell_area)) # Average by time step.

  return(Summary)
}

#' Get Surface Irradiance & Air Temperature
#'
#' This function reads either title variable from a NEMO-MEDUSA model \strong{DRIVERS} file, and returns a monthly time series.
#'
#' The appropriate variable in the netcdf file is imported according to the `Type` parameter, only reading within an
#' x/y window specified in `Space`. The function then drops points outside the model domain before constructing the monthly
#' time series.
#'
#' Unlike other NEMO-MEDUSA related get_* functions, this function constructs the monthly time series directly. A wrapper function
#' to handle time isn't required, as each netcdf file for driving data contains 360 day steps for a single year. `stratify` also
#' isn't called, as these variables have no depth dimension.
#'
#' @param File The location of the netcdf file containing NEMO-MEDUSA driving data.
#' @param Type The variable contained in the netcdf file either "T150" (air temperature) or "SWF" (surface irradiance).
#' @param Year The year the necdf file contains data for.
#' @return A dataframe containing a monthly time series within a year of either average air temperature or surface
#' irradiance. Air temperature is also split by shore zone.
#' @importFrom data.table := as.data.table setnames
#' @family NEMO-MEDUSA variable extractors
#' @export
get_air_dt <- function(File, Type, Year) {

  #File <- all_files$File[1] ; Type <- all_files$Type[1] ; Year <- all_files$Year[1] # test
  if(Type == "SWF") months <- Light_months                                 # Get the file holding the months
  if(Type == "T150") months <- Airtemp_months                              # For the time steps of this data

nc_raw <- ncdf4::nc_open(File)                                             # Open up a netcdf file to see it's raw contents (var names)
nc_var <- ncdf4::ncvar_get(nc_raw, Type, c(Space$Limits$Lon_start, Space$Limits$Lat_start, 1),  # Extract the variable of interest
                    c(Space$Limits$Lon_count, Space$Limits$Lat_count, -1)) # cropped to window, with all time steps
ncdf4::nc_close(nc_raw)                                                    # You must close an open netcdf file when finished to avoid data loss

DT <- as.data.table(nc_var, value.name = "Measured") %>%                   # Pull array
  setnames(old = c("V1", "V2", "V3"), new = c("Longitude", "Latitude", "Time_step")) %>% # Name the columns
  .[, ':='(Longitude = Space$Lons[Longitude],                              # read ':=' as mutate
           Latitude = Space$Lats[Latitude],                                # Replace the factor levels with dimension values
           Month = months[Time_step, "Month"],                             # Assign months to time steps
           Year = Year,                                                    # Add year
           Type = Type)]                                                   # Add variable name
setkey(DT, Longitude, Latitude)                                            # Set key columns for speedy joins

DT <- DT[domains_mask] #%>%
  #merge(domains_mask, all.y = TRUE)                                        # Crop to domain

## Variable specific summaries

if(Type == "SWF") Data <- DT[,by = .(Month, Year, Type),                   # We don't need to bother accounting for shore in light data
                             .(Measured = weighted.mean(Measured, Cell_area))]  # Average by time ste, weighted by cell area
if(Type == "T150") Data <- DT[,by= .(Month, Year, Type, Shore),            # We care about shore for temperature, retain interesting columns
                             .(Measured = weighted.mean(Measured, Cell_area))]  # Average by time ste, weighted by cell area
return(Data)
}
