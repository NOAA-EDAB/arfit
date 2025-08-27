#' Create Lazy data for OISS bottom temp
#'
#' from ecodata::seasonal_oisst_anom
#'

create_surface_temp_oisst_data <- function() {
  dataSet <- readRDS(
    here::here("data-raw", "surface_temp_oisst.rds")
  )

  surface_temp_oisst <- dataSet
  usethis::use_data(surface_temp_oisst, overwrite = TRUE)
}
