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

#' Get weights for interpolation by inverse distance weighting
#'
#' Interpolation weights for inverse distance weighting.
#'
#' @param dat A matrix of the points for which the function is known
#' @param xinp Coordinates, where interpolation has to be computed
#' @return The weights for the IDW interpolation
#' @export
si_get <- function(data, xinp) {
  #d <- dim(smpl@coords)[2] #4 #Anzahl der Koordinaten  aus SPDF structure Auslesen
  d <- length(xinp) #4 #Anzahl der Koordinaten  aus Länge von xinp auslesen
  if (is(data, "SpatialPointsDataFrame")==TRUE) {
    data <- as.matrix(data@data)
  }
  N <- dim(data)[1] #Anzahl der Punkte
  iniR <- getR(d=d, N=N)
  #Abstände vom Interpolationspunkt:
  di <- as.matrix(dist(rbind(xinp,data[,1:d])))[-1,1]
  #für welche Indices sind die Di näher am Interpolationspunkt als der Initialradius?
  cp <- which(di<iniR ,arr.ind=TRUE)
  ncp <- length(cp)
  # wie viele Punkte benutze ich also? min 4 und max 10:
  n <- max(4,min(10,ncp))   # minimale und maximale Interpolationspunktzahl, evtl anders wählen
  #supscripts ij for increasing distance (die sortierten)
  ij <- sort(di, index.return=TRUE)$ix
  # Finaler Radius r:                     ########### evtl extra 0/1-Gewichtungsfaktor?, wenig hilfreich wegen fallunterscheidung
  if (ncp<=4) { radius<-di[ij][4+1]
  } else {
    if (ncp<=10) {radius<-iniR
    } else {
      radius<-di[ij][10+1]
    }
  }
  # Gewichte berechnen
  #si: (for N data points)
  si <- numeric(N)
  si <- si+(di>0)*(di<=radius/3)*(1/di)
  si <- si+(di>radius/3)*(di<=radius)*(27/(4*radius))*(((di/radius)-1)^2)
  return(si)  # Anmerkung: sum(si) nicht 1
}
