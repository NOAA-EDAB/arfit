#' Run test on bottom temperature example data set
#'
#' Performs the short term trend analysis on bundled data set
#'
#' @param epu character string. Name of the EPU to test ("GB", "GOM")
#' @param nBootSims Numeric scalar. Number of bootstrap samples to perform
#'
#'
#'@export

example_bottom_temp <- function(
  epu = "GB",
  nBootSims = 999
) {
  dataSet <- arfit::bottom_temp_survey |>
    dplyr::filter(EPU == epu) |>
    dplyr::select(Time, Value) |>
    dplyr::rename(x = Time, y = Value)

  res <- arfit::fit_real_data(dataSet, nBootSims = nBootSims, printFig = F)

  return(res)
}
