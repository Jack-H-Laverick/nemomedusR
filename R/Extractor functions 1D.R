
#### NEMO - MEDUSA file extractors for taking a water column  ####

#' Get Salinity & Temperature
#'
#' This function reads in the title variables from NEMO-MEDUSA model outputs and returns a water column.
#'
#' Each variable of interest in the netcdf file is imported, only reading at a single specified pixel for all depths.
#'
#' @param filename The name of a netcdf file containing the title variables.
#' @param date The numeric string representing the date contained by the data file.
#' @param x The x index in the grid to sample at.
#' @param y The y index in the grid to sample at.
#' @return A dataframe containing estimates of all variables at all depths for a single point in the model grid
#' @family NEMO-MEDUSA variable extractors
#' @export
get_grid_T_1D   <- function(filename, date, x, y) {
  
  print(stringr::str_glue("{filename} Extracting Salinity, Temperature"))
  nc_raw <- ncdf4::nc_open(filename)                                # Open up a netcdf file to see it's raw contents (var names)
  nc_saline <- ncdf4::ncvar_get(nc_raw, "vosaline", start=c(x,y,1,1), count=c(1,1,-1,-1))          # Extract an array of salinities
  nc_temp <- ncdf4::ncvar_get(nc_raw, "votemper", start=c(x,y,1,1), count=c(1,1,-1,-1))            # Extract an array of temperatures
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss
  
  all <- data.frame(                                                                # Bind as columns
    Salinity = nc_saline,
    Temperature = nc_temp,
    Date = date) %>% 
    rownames_to_column(var="Depth")
  return(all)
}

#' Get Dissolved Inorganic Nitrogen & Chlorophyll a
#'
#' This function reads in the title variables from NEMO-MEDUSA model outputs and returns a water column.
#'
#' Each variable of interest in the netcdf file is imported, only reading at a single specified pixel for all depths.
#'
#' @param filename The name of a netcdf file containing the title variables.
#' @param date The numeric string representing the date contained by the data file.
#' @param x The x index in the grid to sample at.
#' @param y The y index in the grid to sample at.
#' @return A dataframe containing estimates of all variables at all depths for a single point in the model grid
#' @family NEMO-MEDUSA variable extractors
#' @export
get_ptrc_T_1D <- function(filename, date, x, y) {
  
  print(stringr::str_glue("{filename} Extracting Dissolved Inorganic Nitrogen, Chlorophyll"))
  nc_raw <- ncdf4::nc_open(filename)                                # Open up a netcdf file to see it's raw contents (var names)
  nc_DIN <- ncdf4::ncvar_get(nc_raw, "DIN", start=c(x,y,1,1), count=c(1,1,-1,-1))          # Extract an array of Dissolved inorganic nitrogen
  nc_CHD <- ncdf4::ncvar_get(nc_raw, "CHD", start=c(x,y,1,1), count=c(1,1,-1,-1)) # Extract an array of CHD
  nc_CHN <- ncdf4::ncvar_get(nc_raw, "CHN", start=c(x,y,1,1), count=c(1,1,-1,-1))
  nc_Chl<-nc_CHD+nc_CHN
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss
  
  all <- data.frame(                                                                # Bind as columns
    DIN = nc.din,
    Chlorophyll = nc_Chl,
    Date = date) %>% 
    rownames_to_column(var="Depth")
  return(all)
}

#' Get Vertical Velocity and Vertical Eddy Diffusivity
#'
#' This function reads in the title variables from NEMO-MEDUSA model outputs and returns a water column.
#'
#' Each variable of interest in the netcdf file is imported, only reading at a single specified pixel for all depths.
#'
#' @param filename The name of a netcdf file containing the title variables.
#' @param date The numeric string representing the date contained by the data file.
#' @param x The x index in the grid to sample at.
#' @param y The y index in the grid to sample at.
#' @return A dataframe containing estimates of all variables at all depths for a single point in the model grid
#' @family NEMO-MEDUSA variable extractors
#' @export
get_grid_W_1D   <- function(filename, date, x, y) {
  
  print(stringr::str_glue("{filename} Extracting Vertical Velocity, Vertical Eddy Diffusivity"))
  nc_raw <- ncdf4::nc_open(filename)                                # Open up a netcdf file to see it's raw contents (var names)
  nc_vel <- ncdf4::ncvar_get(nc_raw, "vovercrtz", start=c(x,y,1,1), count=c(1,1,-1,-1))          # Extract an array of velocities
  nc_dif <- ncdf4::ncvar_get(nc_raw, "votkeavt", start=c(x,y,1,1), count=c(1,1,-1,-1))           # Extract an array of diffusivity
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss
  
  all <- data.frame(                                                                # Bind as columns
    vel = nc_vel,
    Chlorophyll = nc_dif,
    Date = date) %>% 
    rownames_to_column(var="Depth")
  return(all)
}


#' Get Meridional Currents
#'
#' This function reads in the title variables from NEMO-MEDUSA model outputs and returns a water column.
#'
#' Each variable of interest in the netcdf file is imported, only reading at a single specified pixel for all depths.
#'
#' @param filename The name of a netcdf file containing the title variables.
#' @param date The numeric string representing the date contained by the data file.
#' @param x The x index in the grid to sample at.
#' @param y The y index in the grid to sample at.
#' @return A dataframe containing estimates of all variables at all depths for a single point in the model grid
#' @family NEMO-MEDUSA variable extractors
#' @export
get_grid_V_1D <- function(filename, date, x, y) {
  
  print(stringr::str_glue("{filename} Extracting Meridonal Currents"))
  nc_raw <- ncdf4::nc_open(filename)                                # Open up a netcdf file to see it's raw contents (var names)
  nc_merid <- ncdf4::ncvar_get(nc_raw, "vomecrty", start=c(x,y,1,1), count=c(1,1,-1,-1))          # Extract an array of meridinal currents
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss
  
  all <- data.frame(                                                                # Bind as columns
    merid = nc_merid,
    Date = date) %>% 
    rownames_to_column(var="Depth")
  return(all)
}

#' Get Zonal Currents
#'
#' This function reads in the title variables from NEMO-MEDUSA model outputs and returns a water column.
#'
#' Each variable of interest in the netcdf file is imported, only reading at a single specified pixel for all depths.
#'
#' @param filename The name of a netcdf file containing the title variables.
#' @param date The numeric string representing the date contained by the data file.
#' @param x The x index in the grid to sample at.
#' @param y The y index in the grid to sample at.
#' @return A dataframe containing estimates of all variables at all depths for a single point in the model grid
#' @family NEMO-MEDUSA variable extractors
#' @export
get_grid_U_1D <- function(filename, date, x, y) {
  
  print(stringr::str_glue("{filename} Extracting Zonal Currents"))
  nc_raw <- ncdf4::nc_open(filename)                                # Open up a netcdf file to see it's raw contents (var names)
  nc_DIN <- ncdf4::ncvar_get(nc_raw, "vozocrtx", start=c(x,y,1,1), count=c(1,1,-1,-1))          # Extract an array of Zonal Current
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss
  
  all <- data.frame(                                                                # Bind as columns
    zonal=nc_zonal,
    Date = date) %>% 
    rownames_to_column(var="Depth")
  return(all)
}

#' Condense Daily netcdf Files into a Month by Type
#'
#' This function takes the metadata for multiple netcdf files and creates a spatial summary for the month.
#' \strong{For one file type only}.
#'
#' The function bridges the step between extracting data from a netcdf file, and creating a dataframe of all variables
#' for a whole month of NEMO-MEDUSA model outputs.
#'
#' Different file types require a different get function. This function takes a collection of netcdf files of
#' the same time from the same month, and passes them to the correct `get_*` for data extraction. The results
#' are passed back to this function.
#' 
#' Depending on the analysis chosen, further operations may be performed. 
#' 
#' - 1D, no further steps.
#' - StrathE2E, estimates from different days at the same lat-lon-depth combination are averaged to get a single 
#' number for the month.
#'
#' Creating intermediate monthly objects allows the terabytes of data to be reduced without encountering memory issues.
#' Also, working on independent monthly packets of data means we can parallelise any data processing for speed.
#'
#' @param data A dataframe containing the metadata of multiple netcdf files which share a type.
#' @param analysis A text string dictating the type of summary to create. Currently supported summaries include "StrathE2E" & "1D".
#' @param ... Additional arguments passed to the relevant `get_*` function.
#' @return The function returns a dataframe containing the monthly average shalllow and deep spatial grids for
#' variables of interest.
#' @family NEMO-MEDUSA variable extractors
#' @export
type_in_month <- function(data, analysis, ...) {

Type <- data[1,3]                                               # Pull type from the file
  
get <- match.fun(paste0("get_", Type, analysis))                # Change the extracting function based on file contents

if(analysis == "StrathE2E") {

  
  Month.type <- purrr::map2(data$Path, data$File, get, ...) %>% # Summarise the contents of each file
  abind::abind(along = 3) %>%                                   # Stack matrices behind each other
  rowMeans(dims = 2)                                            # Quickly take the mean for all variables in the month
}

if(analysis == "1D") {

    Month.type <- purrr::map2(data$filename, data$date, get, ...) %>% # Extract the contents of each file
    data.table::rbindlist()                                     # Fast row bind for the data frames
} 

if(!analysis %in% c("StrathE2E", "1D")) {
  
  print("Whoops, I haven't coded that analysis yet. Did you want a StrathE2E or a 1D summary?")

  }

  return(Month.type)
}

#' Condense Daily netcdf Files into a Monthly Summary
#'
#' This function takes the metadata for multiple netcdf files. It then creates a common spatial summary for all
#' the variables of interest for a month.
#'
#' The function takes the metadata for a collection of files which contain data from the same month. The files
#' are split into data packets which share the same file type, before being passed to `type_in_month` to be summarised.
#' `type_in_month` reduces the NEMO-MEDUSA model outputs from large arrays to summaries depending on the analysis chosen.
#' The summaries for each file type are returned to this function and get bound into a single dataframe. 
#' 
#' StrathE2E analysis returns effectively two matrices. Points outside the project window are
#' removed before saving the dataframe for the month in "./Objects/Months/".
#' 
#' 1D analysis returns a vector of estimates at one pixel across all depths.
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
whole_month <- function(data, analysis, out_dir = "./Objects/Months", crop = NULL, grid = NULL, ...) {

if(analysis == "StrathE2E") {
  Month <- data[1,5] ; Year <- data[1,4]                                    # Pull date

  Month <- split(data, f = list(data$Type)) %>%                             # Split out the files for this month by type, so they can be averaged together
    purrr::map(type_in_month, analysis = analysis, ...) %>%                 # Pull a whole month of data from a single file type
    do.call(cbind, .) %>%                                                   # Join together all the data packets
    cbind(grid) %>%                                                         # Add coordinates and depth labels
    filter(Shore_dist > 0) %>%                                              # Drop points on land
    mutate(Year = as.integer(Year),                                         # Add time
           Month =  as.integer(Month)) %>%
    dplyr::right_join(crop) %>%                                             # Cut out rows outside of plotting window
    saveRDS(., file = paste(out_dir, "/NM", Month, Year, "rds", sep = ".")) # save out a data object for one whole month
}
  
if(analysis == "1D") {

  Month <- split(data, f = list(data$Type)) %>%                             # Split out the files for this month by type, so they can be averaged together
    purrr::map(type_in_month, analysis = analysis, ...) %>%                 # Pull a whole month of data from a single file type
    do.call(cbind, .) %>%                                                   # Join together all the data packets
    mutate(Date = as.POSIXct(Date, format = c("%Y%M%D"))) %>%               # Add time
    saveRDS(., file = stringr::str_glue("{out_dir}/NM.{lubridate::month(.$Date)}.{lubridate::year(.$Date)}.rds")) # save out a data object for one whole month
}
  
if(!analysis %in% c("StrathE2E", "1D")) {
  print("Whoops, I haven't coded that analysis yet. Did you want a StrathE2E or 1D summary?")
  }
}
