#' Create Lazy data for survey bottom temp
#'
#' from ecodata::bottom_temp_insitu
#'

create_bottom_temp_survey_data <- function() {
  dataSet <- readRDS(
    here::here("data-raw", "bottom_temp_survey.rds")
  )

  bottom_temp_survey <- dataSet |>
    dplyr::filter(Time > 2013)
  usethis::use_data(bottom_temp_survey, overwrite = TRUE)
}
