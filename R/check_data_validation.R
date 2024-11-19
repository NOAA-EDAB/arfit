#' Checks data set for misspecification
#'
#' Check column names, missing data, consecutive equally spaced time units
#'
#'

check_data_validation <- function(dataSet) {

  # x and yd names
  if (!all(c("x","y") %in% names(dataSet))) {
    stop("Please specify which columns are x and y. Two of the columns need to have names x and y")
  } else {
    # in case more than 2 columns of data
    dataSet <- dataSet |> dplyr::select(x,y)
  }

  # check for equally spaced time series
  spacing <- diff(dataSet$x)
  if(!all(spacing == spacing[1])){
    # Create continuous time series and add NAs
    x <- data.frame(x = seq(from=min(dataSet$x), to=max(dataSet$x),by = min(spacing)))
    dataSet <- x |>
      dplyr::left_join(dataSet, by = "x")

    warning("Time (x field) is not equally spaced. The missing Time value(s) were added using NA for the response(s)")
  } else {
    # All equally spaced
    # maybe make x = 1:n?
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
