#' Checks data set for misspecification
#'
#' Check column names, missing data, consecutive equally spaced time units
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
    stop("Time (x field) is not equally spaced. You can not apply this test")
  } else {
    # maybe make x = 1:n?
  }

  # check for missing data
  if(any(is.na(dataSet$y))){
    # can not proceed if NA is first or last position
    if(head(is.na(dataSet$y),1) | tail(is.na(dataSet$y),1)) {
      stop("Missing data in your response variable at the begining or end of the time series is not currently permitted")
    } else {
      missingValues <- which(is.na(dataSet$y))
      warning("You have some missing data in your response variable (y).")
    }
  } else {
    missingValues <- NULL
  }

  return(list(dataSet=dataSet,missingValues=missingValues))

}
