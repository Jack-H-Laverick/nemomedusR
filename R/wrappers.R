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
#' @param analysis A text string dictating the type of summary to create. Deprecated, supports "StrathE2E", "1D", & "slabR" for back compatibility.
#' @param time_op A text string dictating the type of summary to create. Currently supported summaries include "mean" & "collect". When using time_op all spatial summaries are controlled by slabR schemes.
#' @param ... Additional arguments passed to the relevant `get_*` function.
#' @return The function returns a dataframe if concatenating, or a matrix if averaging, containing the spatial summaries
#' of NEMO-MEDUSA model output.
#' @export

temporal_operation <- function(data, analysis = NULL, time_op = NULL, ...) {
  
  if(length(analysis > 0)) {  
  if(!analysis %in% c("StrathE2E", "1D", "slabR") & !is.null(analysis)) {
    stop("Whoops, I haven't coded that analysis yet. Did you want a StrathE2E, slabR, or a 1D summary?")}
  
  if(!is.null(analysis) & !is.null(time_op)) {
    stop("Use either 'analysis' or 'time_op' not both.")}

    get <- match.fun(paste0("get_", data[1,"Type"], analysis))                # Change the extracting function based on file contents
  
  if(analysis %in% c("StrathE2E", "slabR")) {

    temporal_summary <- purrr::map2(data$Path, data$File, get, ...) %>% # Summarise the contents of each file
      abind::abind(along = 3) %>%                                 # Stack matrices behind each other
      rowMeans(dims = 2)                                          # Quickly take the mean for all variables in the month
  }
  
  if(analysis == "1D") {
    
    temporal_summary <- purrr::map2(data$filename, data$date, get, ...) %>% # Extract the contents of each file
      purrr::map(as.data.frame) %>%                                 
      data.table::rbindlist()                                               # Fast row bind for the data frames
  }} 

  if(length(time_op > 0)) {  
    
    get <- match.fun(paste0("get_", data[1,"Type"], "slabR"))               # Change the extracting function based on file contents
    
  if(time_op == "mean") {
    
    temporal_summary <- purrr::map2(data$Path, data$File, get, ...) %>% # Summarise the contents of each file
      abind::abind(along = 3) %>%                                 # Stack matrices behind each other
      rowMeans(dims = 2)                                          # Quickly take the mean for all variables in the month
  }
  
  if(time_op == "collect") {
    
    temporal_summary <- purrr::map2(data$Path, data$File, get, ...) %>% # Extract the contents of each file
      purrr::map(as.data.frame) %>%                                 
      purrr::map2(data$Date, cbind) %>%                                     # Add the specific date to each file
      data.table::rbindlist()                                               # Fast row bind for the data frames
  }}
  
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
#' @param data a dataframe containing the metadata of multiple netcdf files from a common month.
#' @param analysis A text string dictating the type of summary to create. Deprecated, supports "StrathE2E", "1D", & "slabR" for back compatibility.
#' @param time_op A text string dictating the type of summary to create. Currently supported summaries include "mean" & "collect". When using time_op all spatial summaries are controlled by slabR schemes.
#' @param crop on the way out with the StrathE2E analysis type.
#' @param out_dir a filepath to the directory to save summarised files in.
#' @param summary a dataframe of metadata to bind to summarised rows. The summary should be arranged by increasing 
#' group number from the summary scheme.
#' @param ... Additional arguments to be passed to get_* functions.
#' @return The function returns a dataframe containing the monthly average shalllow and deep spatial grids for
#' \strong{all} the variables of interest in NEMO-MEDUSA outputs.
#' @family NEMO-MEDUSA variable extractors
#' @export
NEMO_MEDUSA <- function(data, analysis = NULL, time_op = NULL, out_dir = "./Objects/Months", crop = NULL, summary = NULL, ...) {

  if(length(analysis > 0)) {  
  if(!analysis %in% c("StrathE2E", "1D", "slabR") & !is.null(analysis)) {
    stop("Whoops, I haven't coded that analysis yet. Did you want a StrathE2E, slabR, or 1D summary?")}

  if(!is.null(analysis) & !is.null(time_op)) {
    stop("Use either 'analysis' or 'time_op' not both.")}
  
  if(analysis == "StrathE2E") {
    Month <- data[1,5] ; Year <- data[1,4]                                    # Pull date
    
    print("Consider using slabR with scheme_strathE2E instead of StrathE2E")
    
    timestep <- split(data, f = list(data$Type)) %>%                          # Split out the files for this month by type, so they can be averaged together
      purrr::map(temporal_operation, analysis = analysis, ...) %>%            # Pull a whole month of data from a single file type
      do.call(cbind, .) %>%                                                   # Join together all the data packets
      cbind(summary) %>%                                                      # Add coordinates and depth labels
      filter(Bathymetry != 0) %>%                                             # Drop points on land
      mutate(Year = as.integer(Year),                                         # Add time
             Month =  as.integer(Month)) %>%
      dplyr::right_join(crop) %>%                                             # Cut out rows outside of plotting window
      saveRDS(file = paste0(out_dir, "/NM.", Month, ".", Year, ".rds"))       # save out a data object for one whole month
  }
  
  if(analysis == "slabR") {
    Month <- data[1,5] ; Year <- data[1,4]                                    # Pull date
    
    timestep <- split(data, f = list(data$Type)) %>%                          # Split out the files for this month by type, so they can be averaged together
      purrr::map(temporal_operation, analysis = analysis, ...)                # Pull a whole month of data from a single file type
      
    if(length(timestep) > 1) {timestep <- do.call(cbind, timestep)} else {    # Join together all the data packets if there are multiple file types
      timestep <- timestep[[1]] }                                             # If there's only one file type expose the output as a matrix

    timestep <- as.data.frame(timestep) %>% 
      cbind(summary) %>%                                                      # Add coordinates and depth labels
      mutate(Year = as.integer(Year),                                         # Add time
             Month =  as.integer(Month)) %>%
      saveRDS(file = paste0(out_dir, "/NM.", Month, ".", Year, ".rds")) # save out a data object for one whole month
  }
  
  if(analysis == "1D") {
    
    timestep <- split(data, f = list(data$Type)) %>%                          # Split out the files for this month by type, so they can be averaged together
      purrr::map(temporal_operation, analysis = analysis, ...) %>%            # Pull a whole month of data from a single file type
      purrr::reduce(~{                                                        # Join together all the data packets
        data.table::setDT(.x, key = c("Date", "Layer"))                       # Try to speed it up with DT keys
        data.table::setDT(.y, key = c("Date", "Layer"))  
        .y[.x]}) %>%                                              
      mutate(Month = stringr::str_sub(Date, start = 5, end = 6),
             Year = stringr::str_sub(Date, start = 1, end = 4)) %>% 
      saveRDS(file = paste0(out_dir, "/NM.", .$Month[1], ".", .$Year[1], ".rds")) # save out a data object for one whole month
  }}
  
  if(length(time_op > 0)) {  
    
  if(time_op == "collect") {
    Month <- data[1,5] ; Year <- data[1,4]                                    # Pull date
    
    timestep <- split(data, f = list(data$Type)) %>%                          # Split out the files for this month by type, so they can be averaged together
      purrr::map(temporal_operation, time_op = "collect", ...)                # Pull a whole month of data from a single file type
    
    if(length(timestep) > 1) {timestep <- do.call(cbind, timestep)} else {    # Join together all the data packets if there are multiple file types
      timestep <- timestep[[1]] }                                             # If there's only one file type expose the output as a matrix
    
    timestep <- as.data.frame(timestep) %>% 
      cbind(summary) %>%                                                      # Add coordinates and depth labels
      mutate(Year = as.integer(Year),                                         # Add time
             Month =  as.integer(Month)) %>%
      saveRDS(file = paste0(out_dir, "/NM.", Month, ".", Year, ".rds")) # save out a data object for one whole month
  }}
  
}
