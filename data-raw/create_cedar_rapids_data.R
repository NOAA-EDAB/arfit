#' Create Lazy data for cedar rapids
#'
#'
#'

library(magrittr)
create_cedar_rapids_data <- function() {

  dataSet <- readr::read_table(here::here("data-raw","usgs_05464500_cedar_rapids.txt"),
                             skip = 35,
                             col_names = TRUE)

  # pick out 2 columns and rename them
  cedar_rapids <- dataSet %>%
    dplyr::select(year_nu,mean_va) %>%
    dplyr::rename(year = year_nu, riverflow=mean_va)


  usethis::use_data(cedar_rapids, overwrite = TRUE)

}
