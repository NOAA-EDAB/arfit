#' Checks data set for issues
#'
#' Check column names, missing data, equally spaced time units
#' If NAs are found at the start or the end of the time series they are removed
#' since this wont affect the fitting, it will just use fewer data points.
#' If NAs are present in the middle of the time series, the location is returned
#' and is used in the fitting procedure
#'
#' @param dataSet Data frame. Passed internally
#'
#' @return A list
#' \item{dataSet}{The "cleaned" version of the data set}
#' \item{missingValues}{Numeric vector. The location of NAs}
#'
#' @export

check_data_validation <- function(dataSet) {

  # x and y names
  if (!all(c("x","y") %in% names(dataSet))) {
    stop("Please specify which columns are x and y. Two of the columns need to have names x and y")
  } else {
    # in case more than 2 columns of data
    dataSet <- dataSet |> dplyr::select(x,y)
  }

  # check for equally spaced time series
    spacing <- diff(dataSet$x)
  if(!all(spacing == spacing[1])){
    stop("Time (x field) is not equally spaced. You can not apply this test.
         You will need to explicitly add the missing Time value(s) using NA for the response(s)")
  } else {
    # We could "fix" the time series and add NAs
  }


  # recursive incase blocks of NA at start or end of data
  clean <- F
  while(!clean) {
    clean <- T
    # check for missing data
    if(head(is.na(dataSet$y),1) ) {
      # NA in first position. Drop it.
      # Treat the timeseries as length n-1
      dataSet <- dataSet |>
        dplyr::slice_tail(n = -1)
      clean <- F
    }

    if (tail(is.na(dataSet$y),1)) {
      # NA in last position. Drop it
      # Treat the timeseries as length n-1
      dataSet <- dataSet |>
        dplyr::slice_head(n = -1)
      clean <- F
    }

  }


  # check for other missing data
  if(any(is.na(dataSet$y))){
    missingValues <- which(is.na(dataSet$y))
    warning("You have some missing data in your response variable (y)")
  } else {
    missingValues <- NULL
  }

  return(list(dataSet=dataSet,missingValues=missingValues))

}
