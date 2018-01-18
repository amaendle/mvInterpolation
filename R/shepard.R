#' Get initial search-radius
#'
#' This function determines a search radius such that on average 7 points lie
#' within the ball with the radius.
#'
#' @param dat A matrix of the points for which the function is known
#' @param d Dimension of the coordinate space
#' @return The radius as a  numeric
#' @export
getR <- function(dat,d) {
  N <- dim(data)[1] #Number of data points
  #use package geometry for the convex hull
  Avol <- geometry::convhulln(dat[,1:d], options="FA")$vol
  #Avol <- 10^d  #alternatively use smpl@bbox
  return((7*Avol/(N*pi^(d/2)))^(1/d))
}
