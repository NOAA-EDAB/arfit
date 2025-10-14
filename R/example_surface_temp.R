#' Run test on example OISST data set
#'
#' Performs the short term trend analysis on bundled data set
#'
#' @param epu character string. Name of the EPU to test ("GB", "GOM")
#' @param nBootSims Numeric scalar. Number of bootstrap samples to perform
#'
#'
#'@export

example_surface_temp <- function(
  epu = "GB",
  nBootSims = 999
) {
  dataSet <- arfit::surface_temp_oisst |>
    dplyr::filter(EPU == epu) |>
    dplyr::select(Time, Value) |>
    dplyr::rename(x = Time, y = Value)

  res <- arfit::fit_real_data(dataSet, nBootSims = nBootSims, printFig = F)

  return(res)
}
