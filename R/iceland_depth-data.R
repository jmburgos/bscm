#' Bottom depth around Iceland
#'
#' A grid of bottom depth values around Iceland.  The grid has 89 longitude
#' values between -27.8 and -10.2 degrees, and 69 latitude values between
#' 61.1 and 67.9 degrees.  Bottom depth values were obtained from GEBCO
#' (General Bathymetric Chart of the Oceans).  Locations on land are coded
#' as NA.
#'
#' @format A named list with three components:
#' \describe{
#'   \item{lon}{longitude, in decimall degrees}
#'   \item{lat}{latitude, in decimal degrees}
#'   \item{depth}{bottom depth, in meters}
#' }
#' @source \url{https://www.gebco.net/}
"iceland_depth"
