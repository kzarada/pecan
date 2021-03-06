context("met2model")

outfolder <- tempfile()
setup(dir.create(outfolder, showWarnings = FALSE))
teardown(unlink(outfolder, recursive = TRUE))

test_that("Met conversion runs without error", {
  skip(paste0(
    "BIOCRO met2model is currently broken. ",
    "See issue #2274 (https://github.com/PecanProject/pecan/issues/2274)."
  ))
  nc_path <- system.file("test-data", "CRUNCEP.2000.nc",
                         package = "PEcAn.utils")
  in.path <- dirname(nc_path)
  in.prefix <- "CRUNCEP"
  start_date <- "2000-01-01"
  end_date <- "2000-12-31"
  result <- met2model.BIOCRO(in.path, in.prefix, outfolder,
                             lat = 45.25,
                             lon = -84.75,
                             start_date = start_date,
                             end_date = end_date)
  expect_s3_class(result, "data.frame")
  expect_true(file.exists(result[["file"]][[1]]))
})
