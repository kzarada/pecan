# ----------------------------------------------------------------------
#------------------------------------------ Load required libraries-----
# ----------------------------------------------------------------------
library("PEcAn.all")
library("PEcAn.utils")
library("RCurl")
library("REddyProc")
library("tidyverse")
library("furrr")
library("R.utils")
library("dynutils")
plan(multiprocess)

# ----------------------------------------------------------------------------------------------
#------------------------------------------ That's all we need xml path and the out folder -----
# ----------------------------------------------------------------------------------------------
outputPath <- "/fs/data3/kzarada/ouput"
nodata <- FALSE
xmlTempName <-"gefs.sipnet.template.xml"
restart <-TRUE
days.obs <- 3  #how many of observed data to include -- not including today

setwd(outputPath)
#------------------------------------------------------------------------------------------------
#------------------------------------------ sourcing the required tools -------------------------
#------------------------------------------------------------------------------------------------
c(
  'Utils.R',
  'download_WCr.R',
  "gapfill_WCr.R",
  'prep.data.assim.R'
) %>% walk( ~ source(
  system.file("WillowCreek",
              .x,
              package = "PEcAn.assim.sequential")
))
#reading xml
settings <- read.settings("/fs/data3/kzarada/pecan/modules/assim.sequential/inst/WillowCreek/gefs.sipnet.template.xml")

#connecting to DB
con <-try(PEcAn.DB::db.open(settings$database$bety), silent = TRUE)

#------------------------------------------------------------------------------------------------
#------------------------------------------ Preparing the pecan xml -----------------------------
#------------------------------------------------------------------------------------------------
#--------------------------- Finding old sims
all.previous.sims <- list.dirs(outputPath, recursive = F)
all.previous.sims <- all.previous.sims[-c(which(all.previous.sims == "/fs/data3/kzarada/ouput/LAI_Reanalysis"),which(all.previous.sims == "/fs/data3/kzarada/ouput/Reanalysis"))] 
if (length(all.previous.sims) > 0 & !inherits(con, "try-error")) {
  
  tryCatch({
    # Looking through all the old simulations and find the most recent
    all.previous.sims <- all.previous.sims %>%
      map(~ list.files(path = file.path(.x, "SDA"))) %>%
      setNames(all.previous.sims) %>%
      discard( ~ !"sda.output.Rdata" %in% .x) # I'm throwing out the ones that they did not have a SDA output
    
    last.sim <-
      names(all.previous.sims) %>%  
      map_chr( ~ strsplit(.x, "_")[[1]][2]) %>%
      map_dfr(~ db.query(
        query = paste("SELECT started_at FROM workflows WHERE id =", .x),
        con = con
      ) %>% 
        mutate(ID=.x)) %>%
      mutate(started_at = as.Date(started_at)) %>%
      arrange(desc(started_at), desc(ID)) %>%
      head(1)
    # pulling the date and the path to the last SDA
    restart.path <-grep(last.sim$ID, names(all.previous.sims), value = T)
    sda.start <- last.sim$started_at
  },
  error = function(e) {
    restart.path <- NULL
    sda.start <- Sys.Date() - 9
    PEcAn.logger::logger.warn(paste0("There was a problem with finding the last successfull SDA.",conditionMessage(e)))
  })
  
  # if there was no older sims
  if (is.na(sda.start))
    sda.start <- Sys.Date() - 9
}

sda.end <- Sys.Date()
#-----------------------------------------------------------------------------------------------
#------------------------------------------ Download met and flux ------------------------------
#-----------------------------------------------------------------------------------------------
#Fluxes
if(!exists('prep.data'))
  prep.data <- prep.data.assim(
    sda.start - 90,# it needs at least 90 days for gap filling 
    sda.end,
    numvals = 100,
    vars = c("NEE", "LE"),
    data.len = days.obs * 24  
  ) 
obs.raw <-prep.data$rawobs
prep.data<-prep.data$obs


# if there is infinte value then take it out - here we want to remove any that just have one NA in the observed data 
prep.data <- prep.data %>% 
  map(function(day.data){
    #cheking the mean
    nan.mean <- which(is.infinite(day.data$means) | is.nan(day.data$means) | is.na(day.data$means))
    if ( length(nan.mean)>0 ) {
      
      day.data$means <- day.data$means[-nan.mean]
      day.data$covs <- day.data$covs[-nan.mean, -nan.mean] %>%
        as.matrix() %>%
        `colnames <-`(c(colnames(day.data$covs)[-nan.mean]))
    }
    day.data
  })


# Changing LE to Qle which is what sipnet expects
prep.data <- prep.data %>%
  map(function(day.data) {
    names(day.data$means)[names(day.data$means) == "LE"] <- "Qle"
    dimnames(day.data$covs) <- dimnames(day.data$covs) %>%
      map(function(name) {
        name[name == "LE"] <- "Qle"
        name
      })
    
    day.data
  })



# Finding the right end and start date
met.start <- obs.raw$Date%>% head(1) %>% lubridate::floor_date(unit = "day")
met.end <- met.start + lubridate::days(16)

# Download MODIS LAI Data 
tryCatch({
  lai <- call_MODIS(outfolder = '/fs/data3/kzarada/NEFI/MODIS/', 
                  start_date = paste0(lubridate::year(met.start), strftime(met.start, format = "%j")),
                  end_date = paste0(lubridate::year(met.end), strftime(met.end, format = "%j")), 
                  lat = 45.805925,
                  lon = -90.07961, 
                  size = 0, 
                  product = "MOD15A2H", 
                  band = "Lai_500m", 
                  band_qc = "FparLai_QC", 
                  band_sd = "LaiStdDev_500m",
                  siteID = NULL, 
                  package_method = "MODISTools", 
                  QC_filter = FALSE,
                  progress = TRUE)}, 
  error = function(e) {
                  lai <- NA
                  PEcAn.logger::logger.warn(paste0("MODIS Data not available for these dates",conditionMessage(e)))
                }
    )
if(!exists('lai')){lai = NULL}
  
#pad Observed Data to match met data 

date <-
  seq(
    from = lubridate::with_tz(as.POSIXct(met.start, format = "%Y-%m-%d"), tz = "UTC") + lubridate::days(1),
    to = lubridate::with_tz(as.POSIXct(met.end - lubridate::days(1), format = "%Y-%m-%d"), tz = "UTC"),
    by = "6 hour"
  )
pad.prep <- obs.raw %>%
  tidyr::complete(Date = seq(
    from = lubridate::with_tz(as.POSIXct(met.start, format = "%Y-%m-%d"), tz = "UTC") + lubridate::days(1),
    to = lubridate::with_tz(as.POSIXct(met.end - lubridate::days(1), format = "%Y-%m-%d"), tz = "UTC"),
    by = "6 hour"
  )) %>%
  mutate(means = NA, covs = NA) %>%
  dplyr::select(Date, means, covs) %>%
  dynutils::tibble_as_list()

names(pad.prep) <-date

#create the data type to match the other data 
pad.cov <- matrix(data = c(NA, NA, NA, NA, NA, NA, NA, NA, NA ), nrow = 3, ncol = 3, dimnames = list(c("NEE", "Qle", "LAI"), c("NEE", "Qle", "LAI")))
pad.means = c(NA, NA, NA)
names(pad.means) <- c("NEE", "Qle", "LAI")

#cycle through and populate the list 

pad <- pad.prep %>% 
          map(function(day.data){
            day.data$means <- pad.means
            day.data$covs <- pad.cov
            day.data
          })

#remove any obs that are in pad and prep.data

  index <- which(names(pad) %in% names(prep.data))
  if(length(index) > 1){pad <- pad[-index]}


#Add in LAI info to prep.data 

if(is.null(lai)){index <- rep(FALSE, length(names(prep.data)))}else{
  index <- as.Date(names(prep.data)) %in% as.Date(lai$calendar_date)
}


for(i in 1:length(index)){
  
  if(index[i]){
    lai.date <- which(as.Date(lai$calendar_date) == as.Date(prep.data[[i]]$Date))
    LAI <- c(0,0)
    prep.data[[i]]$means <- c(prep.data[[i]]$means, lai$data[lai.date])
    prep.data[[i]]$covs <- rbind(cbind(prep.data[[i]]$covs, c(0, 0)), c(0,0, 0))
    #lai$sd[lai.date]
    
  }else{prep.data[[i]]$means <- c(prep.data[[i]]$means, NA)
  prep.data[[i]]$covs <- rbind(cbind(prep.data[[i]]$covs, c(NA,NA)), c(NA,NA,NA))}
  
  names(prep.data[[i]]$means) <- c("NEE", "Qle", "LAI")
  rownames(prep.data[[i]]$covs) <- c("NEE", "Qle", "LAI")
  colnames(prep.data[[i]]$covs) <- c("NEE", "Qle", "LAI")
  
  
}
#add onto end of prep.data list 

prep.data = c(prep.data, pad)


# This line is what makes the SDA to run daily  ***** IMPORTANT CODE OVER HERE
prep.data<-prep.data %>%
  discard(~lubridate::hour(.x$Date)!=0)


obs.mean <- prep.data %>% map('means') %>% setNames(names(prep.data))
obs.cov <- prep.data %>% map('covs') %>% setNames(names(prep.data))









#-----------------------------------------------------------------------------------------------
#------------------------------------------ Fixing the settings --------------------------------
#-----------------------------------------------------------------------------------------------
#Using the found dates to run - this will help to download mets
settings$run$start.date <- as.character(met.start)
settings$run$end.date <- as.character(last(date))
settings$run$site$met.start <- as.character(met.start)
settings$run$site$met.end <- as.character(met.end)
#info
settings$info$date <- paste0(format(Sys.time(), "%Y/%m/%d %H:%M:%S"), " +0000")
# --------------------------------------------------------------------------------------------------
#---------------------------------------------- PEcAn Workflow -------------------------------------
# --------------------------------------------------------------------------------------------------
#Update/fix/check settings. Will only run the first time it's called, unless force=TRUE
settings <- PEcAn.settings::prepare.settings(settings, force=FALSE)
setwd(settings$outdir)
ggsave(
  file.path(settings$outdir, "Obs_plot.pdf"),
  ploting_fluxes(obs.raw) ,
  width = 16,
  height = 9
)

#Write pecan.CHECKED.xml
PEcAn.settings::write.settings(settings, outputfile = "pecan.CHECKED.xml")
# start from scratch if no continue is passed in
statusFile <- file.path(settings$outdir, "STATUS")
if (length(which(commandArgs() == "--continue")) == 0 && file.exists(statusFile)) {
  file.remove(statusFile)
}
# Do conversions
settings <- PEcAn.workflow::do_conversions(settings, T, T, T)

# Query the trait database for data and priors
if (PEcAn.utils::status.check("TRAIT") == 0) {
  PEcAn.utils::status.start("TRAIT")
  settings <- PEcAn.workflow::runModule.get.trait.data(settings)
  PEcAn.settings::write.settings(settings, outputfile = 'pecan.TRAIT.xml')
  PEcAn.utils::status.end()
} else if (file.exists(file.path(settings$outdir, 'pecan.TRAIT.xml'))) {
  settings <-
    PEcAn.settings::read.settings(file.path(settings$outdir, 'pecan.TRAIT.xml'))
}
# Run the PEcAn meta.analysis
if (!is.null(settings$meta.analysis)) {
  if (PEcAn.utils::status.check("META") == 0) {
    PEcAn.utils::status.start("META")
    PEcAn.MA::runModule.run.meta.analysis(settings)
    PEcAn.utils::status.end()
  }
}
#sample from parameters used for both sensitivity analysis and Ens
get.parameter.samples(settings, ens.sample.method = settings$ensemble$samplingspace$parameters$method)
# Setting dates in assimilation tags - This will help with preprocess split in SDA code
settings$state.data.assimilation$start.date <-as.character(met.start)
settings$state.data.assimilation$end.date <-as.character(met.end - lubridate::hms("06:00:00"))

if (nodata) {
  obs.mean <- obs.mean %>% map(function(x)
    return(NA))
  obs.cov <- obs.cov %>% map(function(x)
    return(NA))
}

# --------------------------------------------------------------------------------------------------
#--------------------------------- Restart -------------------------------------
# --------------------------------------------------------------------------------------------------

#@Hamze - should we add a if statement here for the times that we don't want to copy the path?
# @Hamze: Yes if restart == TRUE 
if(restart == TRUE){
  if(!dir.exists("SDA")) dir.create("SDA",showWarnings = F)

  #Update the SDA Output to just have last time step 
  temp<- new.env()
  load(file.path(restart.path, "SDA", "sda.output.Rdata"), envir = temp)
  temp <- as.list(temp)
  
  #we want ANALYSIS, FORECAST, and enkf.parms to match up with how many days obs data we have
  # +2 for days.obs since today is not included in the number. So we want to keep today and any other obs data 
  if(length(temp$ANALYSIS) > 1){
    for(i in rev((days.obs + 2):length(temp$ANALYSIS))){ 
    temp$ANALYSIS[[i]] <- NULL
  }
 
  for(i in rev((days.obs + 2):length(temp$FORECAST))){
    temp$FORECAST[[i]] <- NULL
  } 
  
  
  for(i in rev((days.obs + 2):length(temp$enkf.params))){
    temp$enkf.params[[i]] <- NULL
  } 
  }

  temp$t = 1 
  
  #change inputs path to match sampling met paths 
  
  for(i in 1: length(temp$inputs$ids)){
    
    temp$inputs$samples[i] <- settings$run$inputs$met$path[temp$inputs$ids[i]]
    
  }
  
  temp1<- new.env()
  list2env(temp, envir = temp1)
  save(list = c("ANALYSIS", 'FORECAST', "enkf.params", "ensemble.id", "ensemble.samples", 'inputs', 'new.params', 'new.state', 'run.id', 'site.locs', 't', 'Viz.output', 'X'),
       envir = temp1, 
       file = file.path(settings$outdir, "SDA", "sda.output.Rdata"))  

  
  
  temp.out <- new.env()
  load(file.path(restart.path, "SDA", 'outconfig.Rdata'), envir = temp.out)
  temp.out <- as.list(temp.out)
  temp.out$outconfig$samples <- NULL
  
  temp.out1 <- new.env()
  list2env(temp.out, envir = temp.out1)
  save(list = c('outconfig'), 
       envir = temp.out1, 
       file = file.path(settings$outdir, "SDA", "outconfig.Rdata"))



#copy over run and out folders 
  
  if(!dir.exists("run")) dir.create("run",showWarnings = F)
  copyDirectory(from = file.path(restart.path, "run/"), 
                to = file.path(settings$outdir, "run/"))
  if(!dir.exists("out")) dir.create("out",showWarnings = F)
  copyDirectory(from = file.path(restart.path, "out/"), 
                to = file.path(settings$outdir, "out/"))
} #restart == TRUE
 # --------------------------------------------------------------------------------------------------
#--------------------------------- Run state data assimilation -------------------------------------
# --------------------------------------------------------------------------------------------------

if(restart == FALSE) unlink(c('run','out','SDA'), recursive = T)

if ('state.data.assimilation' %in% names(settings)) {
  if (PEcAn.utils::status.check("SDA") == 0) {
    PEcAn.utils::status.start("SDA")
    PEcAn.assim.sequential::sda.enkf(
      settings, 
      restart=restart,
      Q=0,
      obs.mean = obs.mean,
      obs.cov = obs.cov,
      control = list(
        trace = TRUE,
        interactivePlot =FALSE,
        TimeseriesPlot =TRUE,
        BiasPlot =FALSE,
        debug =FALSE,
        pause=FALSE
      )
    )
    PEcAn.utils::status.end()
  }
}

  