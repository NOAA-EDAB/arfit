#' Bottom temperature data from NEFSC surveys
#'
#' Most recent 10 years of data from the bottom temperature index which incorporates near-bottom temperature measurements collected on Northeast Fisheries Science Center (NEFSC) surveys between 1977-present. Early measurements were made using surface bucket samples, mechanical bathythermographs and expendable bathythermograph probes, but by 1991 the CTD – an acronym for conductivity temperature and depth – became standard equipment on all NEFSC surveys. Near-bottom refers to the deepest observation at each station that falls within 10 m of the reported water depth. Observations encompass the entire continental shelf area extending from Cape Hatteras, NC to Nova Scotia, Canada, inclusive of the Gulf of Maine and Georges Bank.
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
#' The Data were extracted from ecodata v6.0.1 from the indicator \emph{bottom_temp_insitu}
#'
"bottom_temp_survey"
