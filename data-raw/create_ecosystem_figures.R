## create ecosystem figures

# data are part of package
library(arfit)

png("bottom_temp.png", width = 900, height = 600)
par(mfrow = c(1, 2), mai = c(1, 1.5, 0, 0), oma = c(0, 0, 1, 1))
for (var in c("GB", "GOM")) {
  dataSet <- bottom_temp_survey |>
    dplyr::filter(EPU == var) |>
    dplyr::select(Time, Value) |>
    dplyr::rename(x = Time, y = Value)
  res <- arfit::fit_real_data(dataSet, nBootSims = 999, printFig = F)

  message(paste0("bottom temp in ", var))
  print(res$pValue)

  nT <- nrow(res$data)

  ylab = ""
  plot(
    res$data$x,
    res$data$y,
    type = "l",
    xlab = "",
    xaxt = "n",
    yaxt = "n",
    ylab = ylab,
    lwd = 2,
    ylim = c(-0.5, 1.7)
  )
  if (var == "GB") {
    mtext(
      expression("Bottom Temperature Anomaly (\u00B0 C)"),
      side = 2,
      line = 3,
      cex = 2.5
    )
    text(2014.2, 1.7, "a)", cex = 2)
  } else {
    text(2014.2, 1.7, "b)", cex = 2)
  }

  lines(res$data$x, rep(res$null$betaEst, nT), col = "black", lty = 2, lwd = 2)
  #  lines(dataSet$x,alt$betaEst[1]+alt$betaEst[2]*c(1:nT),col="black",lty=3,lwd=2)
  lines(
    res$data$x,
    res$alt$betaEst[1] + res$alt$betaEst[2] * res$data$x,
    col = "black",
    lty = 3,
    lwd = 2
  )
  axis(
    side = 1,
    at = c(2014:2023),
    labels = as.character(2014:2023),
    cex.lab = 2.5,
    cex.axis = 2
  )
  axis(
    side = 2,
    at = seq(-.5, 1.5, .5),
    labels = as.character(seq(-.5, 1.5, .5)),
    pas = 2,
    cex.lab = 2.5,
    cex.axis = 2
  )
}
dev.off()

png("surface_temp.png", width = 900, height = 600)
par(mfrow = c(1, 2), mai = c(1, 1.5, 0, 0), oma = c(0, 0, 1, 1))
for (var in c("GB", "GOM")) {
  dataSet <- surface_temp_oisst |>
    dplyr::filter(EPU == var) |>
    dplyr::select(Time, Value) |>
    dplyr::rename(x = Time, y = Value)
  res <- arfit::fit_real_data(dataSet, nBootSims = 999, printFig = F)

  message(paste0("Surface temp in ", var))
  print(res$pValue)

  nT <- nrow(res$data)

  plot(
    res$data$x,
    res$data$y,
    type = "l",
    xlab = "",
    xaxt = "n",
    ylab = "",
    yaxt = "n",
    lwd = 2,
    ylim = c(-0.5, 1.7)
  )
  if (var == "GB") {
    mtext(
      expression("Surface Temperature Anomaly (\u00B0 C)"),
      side = 2,
      line = 3,
      cex = 2.5
    )
    text(2014.2, 1.7, "c)", cex = 2)
  } else {
    text(2014.2, 1.7, "d)", cex = 2)
  }

  lines(res$data$x, rep(res$null$betaEst, nT), col = "black", lty = 2, lwd = 2)
  #  lines(dataSet$x,alt$betaEst[1]+alt$betaEst[2]*c(1:nT),col="black",lty=3,lwd=2)
  lines(
    res$data$x,
    res$alt$betaEst[1] + res$alt$betaEst[2] * res$data$x,
    col = "black",
    lty = 3,
    lwd = 2
  )
  axis(
    side = 1,
    at = c(2014:2023),
    labels = as.character(2014:2023),
    cex.lab = 2.5,
    cex.axis = 2
  )
  axis(
    side = 2,
    at = seq(-.5, 1.5, .5),
    labels = as.character(seq(-.5, 1.5, .5)),
    pas = 2,
    cex.lab = 2.5,
    cex.axis = 2
  )
}
dev.off()
