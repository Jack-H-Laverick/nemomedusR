#' Perform Temporal Summary Operations Within `NEMO-MEDUSA` by Model Output File Type
#'
#' This function takes the metadata for multiple netcdf files, passes it to be summarised spatially according to 
#' `analysis`, and then performs a temporal operation, either binding the time steps together, or averaging them.
#'
#' The function bridges the step between extracting data from a netcdf file, and creating a summary of all variables
#' for a chunk of NEMO-MEDUSA model outputs.
#'
#' Different file types require a different get function. This function takes a collection of netcdf files of
#' the same time from the same month, and passes them to the correct `get_*` for data extraction. The results
#' are passed back to this function for handling the time dimension.
#' 
#' Depending on the analysis chosen, operations performed are: 
#' 
#' - 1D, no further steps.
#' - StrathE2E or slabR, estimates from different days at the same lat-lon-depth combination are averaged to get a single 
#' number for the target time step.
#'
#' Creating intermediate summary objects allows the terabytes of data to be reduced without encountering memory issues.
#' Also, working on independent packets of data means we can parallelise any data processing for speed.
#'
#' @param data A dataframe containing the metadata of multiple netcdf files.
#' @param analysis A text string dictating the type of summary to create. Currently supported summaries include "StrathE2E" & "1D".
#' @param ... Additional arguments passed to the relevant `get_*` function.
#' @return The function returns a dataframe if concatenating, or a matrix if averaging, containing the spatial summaries
#' of NEMO-MEDUSA model output.
#' @export

temporal_operation <- function(data, analysis, ...) {
  
  if(!analysis %in% c("StrathE2E", "1D", "slabR")) {
    stop("Whoops, I haven't coded that analysis yet. Did you want a StrathE2E, slabR, or a 1D summary?")}
  
  Type <- data[1,3]                                               # Pull type from the file
  get <- match.fun(paste0("get_", Type, analysis))                # Change the extracting function based on file contents
  
  if(analysis %in% c("StrathE2E", "slabR")) {

    temporal_summary <- purrr::map2(data$Path, data$File, get, ...) %>% # Summarise the contents of each file
      abind::abind(along = 3) %>%                                 # Stack matrices behind each other
      rowMeans(dims = 2)                                          # Quickly take the mean for all variables in the month
  }
  
  if(analysis == "1D") {
    
    temporal_summary <- purrr::map2_df(data$filename, data$date, get, ...) #%>% # Extract the contents of each file
    #data.table::rbindlist()                                      # Fast row bind for the data frames
  } 

  return(temporal_summary)
}

#' Condense NEMO MEDUSA Model Outputs into a Summary file
#'
#' This function takes the metadata for multiple netcdf files. It then creates a single summary file for all
#' variables according to an analysis argument.
#'
#' The function takes the metadata for a collection of files which contain data from the same month. The files
#' are split into data packets which share the same file type, before being passed to `temporal_operation` to be summarised.
#' `temporal_operation` reduces the NEMO-MEDUSA model outputs from large arrays to summaries depending on the analysis chosen.
#' The summaries for each file type are returned to this function and get bound into a single dataframe. 
#' 
#' StrathE2E analysis returns effectively two matrices. Points outside the project window are
#' removed before saving the dataframe for the month in "./Objects/Months/".
#' 
#' slabR uses a C++ routine to summarise arrays according to indices and a grouping scheme 
#' before saving the dataframe for the month in "./Objects/Months/".
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
NEMO_MEDUSA <- function(data, analysis, out_dir = "./Objects/Months", crop = NULL, grid = NULL, ...) {

  if(!analysis %in% c("StrathE2E", "1D", "slabR")) {
    stop("Whoops, I haven't coded that analysis yet. Did you want a StrathE2E, slabR, or 1D summary?")}

  if(analysis == "StrathE2E") {
    Month <- data[1,5] ; Year <- data[1,4]                                    # Pull date
    
    print("Consider using slabR with scheme_strathE2E instead of StrahE2E")
    
    Month <- split(data, f = list(data$Type)) %>%                             # Split out the files for this month by type, so they can be averaged together
      purrr::map(temporal_operation, analysis = analysis, ...) %>%            # Pull a whole month of data from a single file type
      do.call(cbind, .) %>%                                                   # Join together all the data packets
      cbind(grid) %>%                                                         # Add coordinates and depth labels
      filter(Bathymetry != 0) %>%                                             # Drop points on land
      mutate(Year = as.integer(Year),                                         # Add time
             Month =  as.integer(Month)) %>%
      dplyr::right_join(crop) %>%                                             # Cut out rows outside of plotting window
      saveRDS(file = paste0(out_dir, "/NM.", Month, ".", Year, ".rds"))       # save out a data object for one whole month
  }
  
  if(analysis == "slabR") {
    Month <- data[1,5] ; Year <- data[1,4]                                    # Pull date
    
    Month <- split(data, f = list(data$Type)) %>%                             # Split out the files for this month by type, so they can be averaged together
      purrr::map(temporal_operation, analysis = analysis, ...) %>%            # Pull a whole month of data from a single file type
      do.call(cbind, .) %>%                                                   # Join together all the data packets
      cbind(grid) %>%                                                         # Add coordinates and depth labels
      mutate(Year = as.integer(Year),                                         # Add time
             Month =  as.integer(Month)) %>%
      saveRDS(file = paste0(out_dir, "/NM.", Month, ".", Year, ".rds")) # save out a data object for one whole month
  }
  
  if(analysis == "1D") {
    
    Month <- split(data, f = list(data$Type)) %>%                             # Split out the files for this month by type, so they can be averaged together
      purrr::map(temporal_operation, analysis = analysis, ...) %>%            # Pull a whole month of data from a single file type
      purrr::reduce(left_join, by = c("Date", "Depth")) %>%                   # Join together all the data packets
      mutate(Month = stringr::str_sub(Date, start = 5, end = 6),
             Year = stringr::str_sub(Date, start = 1, end = 4)) %>% 
      saveRDS(file = paste0(out_dir, "/NM.", .$Month[1], ".", .$Year[1], ".rds")) # save out a data object for one whole month
  }
}
