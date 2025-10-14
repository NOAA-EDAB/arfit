#' OISST temperature data
#'
#' Most recent 10 years of data from National Oceanographic and Atmospheric Administration optimum interpolation sea surface temperature high resolution data set (NOAA OISST V2)
#'
#'
#' @format A data frame
#'
#' \describe{
#'   \item{Time}{year}
#'   \item{Var}{name of the variable}
#'   \item{Value}{value of the variable}
#'   \item{EPU}{Ecological Production Unit in which Var was measured/calculated}
#'   \item{Units}{units measures}
#' }
#'
#' @source \url{https://noaa-edab.github.io/ecodata/index.html}
#'
#' @family data
#'
#' @section Data:
#'
#' The Data were extracted from ecodata v6.0.1 from the indicator \emph{seasonal_oisst_anom}
#'
"surface_temp_oisst"
