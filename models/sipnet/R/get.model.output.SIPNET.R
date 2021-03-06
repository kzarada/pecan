#-------------------------------------------------------------------------------
# Copyright (c) 2012 University of Illinois, NCSA.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the 
# University of Illinois/NCSA Open Source License
# which accompanies this distribution, and is available at
# http://opensource.ncsa.illinois.edu/license.html
#-------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------#
##' Function to retrieve model output from local or remote server
##'
##' @name get.model.output.SIPNET
##' @title Retrieve model output from local or remote server
##' 
##' @import PEcAn.utils
##' @export
get.model.output.SIPNET <- function(settings) {
  
  model <- "SIPNET"
  
  ### Get model output on the localhost
  if (settings$host$name == "localhost") {
    get.results(settings)
  } else {
    
    ## model output is on a remote host
    remoteScript <- paste(settings$outdir, "PEcAn.functions.R", sep = "")
    
    ### Make a copy of required functions and place in file PEcAn.functions.R
    dump(c(paste("model2netcdf", model, sep = "."), "get.run.id", "read.ensemble.output", "read.sa.output", 
      "read.output", "get.results"), file = remoteScript)
    
    ### Add execution of get.results to the end of the PEcAn.functions.R file This will execute all the
    ### code needed to extract output on remote host
    cat("get.results()", file = remoteScript, append = TRUE)
    
    ### Copy required PEcAn.functions.R to remote host
    rsync("-outi", remoteScript, paste(settings$host$name, ":", settings$host$outdir, sep = ""))
    
    ### Run script on remote host
    system(paste("ssh -T", settings$host$name, "'", "cd", settings$host$outdir, "; R --vanilla < PEcAn.functions.R'"))
    
    ### Get PEcAn output from remote host
    rsync("-outi", from = paste(settings$host$name, ":", settings$host$outdir, "output.Rdata", 
      sep = ""), to = settings$outdir)
  }  ### End of if/else
  
} # get.model.output.SIPNET
