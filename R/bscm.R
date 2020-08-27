#' Spatial interpolation using Bratseth's Succesive Correction Method
#'
#'
#' This function performs a statistical interpolation using Bratseth's
#' Successive Correction Method.  It is designed specifically for marine data,
#' and it uses bottom depth to calculate the similarity among the interpolation
#' points.
#'
#' @param fo numeric vector with the observations of the variable to be
#'  interpolated.
#' @param xo numeric vector with the longitude of the observations.
#' @param yo numeric vector with the latitude of the observations.
#' @param dypobs numeric vector with the bottom depth at the location of the
#'  observations.
#' @param x numeric vector with the longitude values of the interpolation grid.
#' @param y numeric vector with the latitude values of the interpolation grid.
#' @param dypxy numeric matrix with the bottom depth at the nodes of the
#'  interpolation grid.
#' @param radius numeric.  This is the topography parameter as described in
#' Skagseth abd Mork (2012). Defaults to 300m.
#' @param obserror numeric. The target observation error. Defaults to 0.05
#' @param b numeric vector.  Defaults to seq(200, 100, -10)
#' @param fb numeric matrix. Background field (not yet implemented)
#'
#'
#' @return Named list with the following elements:

#' \itemize{
#'   \item f: The values of the interpolated variable on the grid nodes.
#'   \item E: The error variance of the estimate (from 0 to 1)
#'   \item MSEobs: The Mean Square Error on the observation points in each iteration.
#'   \item ndist_ao: The number of observations within each influence value.
#'}
#'
#'@section Details:
#'Naturally, the numeric vectors \code{fo}, \code{xo}, \code{yo} and dypobs must have the same
#'length.  The length of vectors \code{x} and \code{y} should match the number of columns
#'and rows in matrix \code{dypxy}.
#'
#'
#'
#' @examples
#' sum(1:10)

bscm <- function(fo, xo, yo, dypobs = NULL, x, y, dypxy = NULL,
                 radius = 300, obserror = 0.05, fb = NULL, b=NULL) {

  ## Check input

  if (is.null(b)) b <- c(seq(200, 100, -10), rep(100, 9)) # default value

  nx <- length(x)
  ny <- length(y)

  AntObs <- length(xo)
  print(paste(AntObs, "points used"))

  ## Parameters
  d2r <- pi / 180
  JorRad2 <- 6370 * 6370
  cos_br <- d2r * cos(mean(yo) * d2r)

  niter <- length(b) # Number of iterations

  ## Calculate distances between observations
  ObsRad2 <- matrix(NA, ncol = AntObs, nrow = AntObs)

  for (j in 1:AntObs) {

    dr <- ((xo[j] - xo) * cos_br) ^ 2 + ((yo[j] - yo) * d2r) ^ 2

    if (!is.null(dypobs) & !is.null(dypxy)) {
      dr_h <- (abs(radius * (dypobs[j] - dypobs) / (dypobs[j] + dypobs))) ^ 2
      dr_h[is.na(dr_h)] <- 0
      ObsRad2[, j] <- JorRad2 * dr + dr_h
    }
    else {
      ObsRad2[, j] <- JorRad2 * dr
    }

  }

  ## Asuming I don't have a background field
  fb <- matrix(mean(fo), ncol = nx, nrow = ny)
  fbo <- mean(fo)
  fgo <- rep(mean(fo), times = length(xo))

  MSEobs <- rep(NA, times = niter + 1)
  MSEobs[1] <- (sum((fgo-fo)^2)) / AntObs


  ## Calculate distances between observation and predictions

  AnaRad2 <- matrix(NA, nrow = nx * ny, ncol = AntObs)
  ndist_ao <- matrix(NA, nrow = ny, ncol = nx)
  fbr <- rep(NA, length = nx * ny)

  for (k in 1:ny){
    for (i in 1:nx){
      r <- (k-1)*nx + i
      drx <- ((x[i] - xo) * cos_br) ^ 2 + ((y[k] - yo) * d2r) ^ 2
      drxh <- (abs(radius * (dypxy[k, i] - dypobs) / (dypxy[k, i] + dypobs))) ^ 2
      drxh <- tidyr::replace_na(drxh, 0)
      AnaRad2[r, ] <- JorRad2 * drx + drxh

      ndist_ao[k, i] <- sum(sqrt(AnaRad2[r, ]) < dplyr::last(b))
      fbr[r] <- fb[k, i]
    }}


  fr <- fbr
  db <- c(1, diff(b))


  ## Start iterations

  for (n in 1:niter){

    b2 <- 1 / (b[n] * b[n])

    if (db [n] != 0) {

      ObsCov <- exp(-b2 * ObsRad2)

      diag(ObsCov) <- diag(ObsCov) + obserror

      Wj <- rep(NA, AntObs)
      M <- rep(NA, AntObs)

      for (j in 1:AntObs) {
        for (i in 1:AntObs) {
          nev <- ObsCov[j, j] * ObsCov[i, i]
          Wj[i] <- ObsCov[j, i] / sqrt(nev)
        }
        M[j] <- ObsCov[j, j] * sum(Wj)
      }

    }

    a <- matrix(NA, ncol = AntObs, nrow = AntObs)
    ax <- matrix(NA, nrow = nx * ny, ncol = AntObs)

    for (j in 1:AntObs) {

      a[, j]  <- ObsCov[, j] / M[j]
      ax[, j] <- exp(-b2 * AnaRad2[, j]) / M[j]
    }

    fr  <-  fr + ax %*% (fo - fgo)
    fgo <-  fgo + a %*% (fo - fgo)

    MSEobs[n + 1] <- (sum((fgo - fo)^2)) / AntObs

  }


  ##%%%% Beregner Feil relativ til covariansen ****
  b2 <- 1 / (last(b) * last(b))

  ainv <- MASS::ginv(a)

  p <- ax %*% ainv

  er <- rep(NA, nx * ny)
  for (k in 1:ny) {
    for (i in 1:nx) {
      r <- (k - 1) * nx + i
      er[r] <- 1 - sum(M * ax[r, ] * p[r, ])

    }
  }

  ## Sjekker feilen

  inderr <- which(er < 0)
  if (min(er) < (-0.1)){
    print("WARNING: er<0 in the following locations:")
    print(paste(as.character(inderr), collapse = ", "))
    print(paste("min(er) =", round(min(er), 4)))
  }

  if(any(inderr)){
    er[inderr] <- 0
  }

  print("Mean square error of observations:")
  print(round(MSEobs, 5))

  f <- matrix(NA, ncol = nx, nrow = ny)
  E <- matrix(NA, ncol = nx, nrow = ny)

  ## Create grid
  for (k in 1:ny) {

    f[k, ] <- fr[(1 + (k - 1) * nx):(k * nx)]
    E[k, ] <- er[(1 + (k - 1) * nx):(k * nx)]

  }

  return(list(f = f, E = E, MSEobs = MSEobs, ndist_ao = ndist_ao))

}
