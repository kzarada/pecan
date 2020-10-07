##' @title Script to launch NOAA download and temporal downscaling
##' @return None
##' @param site_id, site identification number 
##' @param lat_list, vector of latitudes that correspond to site codes
##' @param lon_list, vector of longitudes that correspond to site codes
##' @param output_directory, directory where the model output will be save
##' @param downscale, logical specifying whether to downscale from 6-hr to 1-hr
##' @param overwrite, logical stating to overwrite any existing output_file

##' @export
##' @export
##'
##' @author Quinn Thomas ; modified by Katie Zarada
##'
##'


download.NOAA_GEFS <- function(site_id,
                               lat.in,
                               lon.in,
                               outfolder,
                               forecast_time = "0",
                               start_date= Sys.Date(),
                               downscale = TRUE,
                               overwrite = FALSE){
  
  model_name <- "NOAAGEFS_6hr"
  model_name_ds <-"NOAAGEFS_1hr" #Downscaled NOAA GEFS
  model_name_raw <- "NOAAGEFS_raw"
  
  PEcAn.logger::logger.info(paste0("Downloading GEFS for site ", site_id, " for ", start_date))
  
  PEcAn.logger::logger.info(paste0("Overwrite existing files: ", overwrite))
 
    
  PEcAn.data.atmosphere::noaa_grid_download(lat_list = lat.in,
                                            lon_list = lon.in,
                                            forecast_time = forecast_time,
                                            forecast_date = start_date,
                                            model_name_raw = model_name_raw,
                                            num_cores = num_cores,
                                            output_directory = outfolder)
    
  PEcAn.data.atmosphere::process_gridded_noaa_download(lat_list = lat.in,
                                                      lon_list = lon.in,
                                                      site_id = site_id,
                                                      downscale = downscale,
                                                      overwrite = overwrite,
                                                      model_name = model_name,
                                                      model_name_ds = model_name_ds,
                                                      model_name_raw = model_name_raw,
                                                      num_cores = num_cores,
                                                      output_directory = output_folder)
  
}