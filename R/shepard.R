#' Generate R documentation from inline comments.
#'
#' Roxygen2 allows you to write documentation in comment blocks co-located
#' with code.
#'
#' The only function Otherwise refer to the vignettes to see
#' how to format the documentation.
#'
#' @author Andreas Mändle, \email{maendle@@uni-bremen.de}
#' @references Shepard, Donald (1968). "A two-dimensional interpolation function for irregularly-spaced data". Proceedings of the 1968 ACM National Conference. pp. 517–524. doi:10.1145/800186.810616.
#' @keywords multivariate interpolation
"_PACKAGE"
#> [1] "_PACKAGE"


#' Get initial search-radius
#'
#' This function determines a search radius such that on average 7 points lie
#' within the ball with the radius.
#'
#' @param data A matrix of the points for which the function is known
#' @param d Dimension of the coordinate space
#' @return The radius as a  numeric
#' @examples
#' coords <- matrix(runif(10*4)*10, ncol=4)
#' getR(coords)
#' @export
getR <- function(data,d=NULL) {
  N <- dim(data)[1] #Number of data points
  if (is.null(d)) d<-dim(data)[2]
  #use package geometry for the convex hull
  Avol <- geometry::convhulln(data[,1:d], options="FA")$vol
  #Avol <- 10^d  #alternatively use smpl@bbox
  return((7*Avol/(N*pi^(d/2)))^(1/d))
}

#' Get weights for interpolation by inverse distance weighting
#'
#' Interpolation weights for inverse distance weighting.
#'
#' @param data A matrix of the points for which the function is known
#' @param xinp Coordinates, where interpolation has to be computed
#' @param p additional parameter as in the Shepard paper
#' @return The weights for the IDW interpolation
#' @examples
#' mydata <- matrix(runif(10*4)*10, ncol=4)
#' mydata <- cbind(mydata,abs(apply(mydata,1,sum)-3),abs(apply(mydata,1,prod)-4))
#' di_get(mydata, c(4,4,4,4))
#' @export
di_get <- function(data, xinp, p=2) {
  d <- length(xinp) # get number of coordinates from xinp
  if (is(data, "SpatialPointsDataFrame")==TRUE)
    data <- as.matrix(data@data)
  N <- dim(data)[1] #number of points

  # weights based on all the distances from xinp:
  di <- as.matrix(dist(rbind(xinp,data[,1:d])))[-1,1]
  #final weights to be used
  return(1/(di^p))
}

#' Get weights for interpolation by inverse distance weighting (IDW) for nearest neighbors
#'
#' Interpolation weights for inverse distance weighting of the nearest neighbors (4 to 10 neighbors, depending on the distance)
#'
#' @param data A matrix of the points for which the function is known
#' @param xinp Coordinates, where interpolation has to be computed
#' @param rad Function which gives the initial radius for the observations considered in the interpolation
#' @return The weights for the IDW interpolation based on the nearest neighbors
#' @examples
#' mydata <- matrix(runif(10*4)*10, ncol=4)
#' mydata <- cbind(mydata,abs(apply(mydata,1,sum)-3),abs(apply(mydata,1,prod)-4))
#' si_get(mydata, c(4,4,4,4))
#' @export
si_get <- function(data, xinp, rad=getR) {
  d <- length(xinp)  #Anzahl der Koordinaten  aus Länge von xinp auslesen
  if (is(data, "SpatialPointsDataFrame")==TRUE) {
    data <- as.matrix(data@data)
  }
  N <- dim(data)[1] #number of points
  if (is.numeric(rad)) {iniR <- rad
    } else {iniR <- rad(dat=data, d=d)}
  #distances from interpolation point
  di <- as.matrix(dist(rbind(xinp,data[,1:d])))[-1,1]
  #für welche Indices sind die Di näher am Interpolationspunkt als der Initialradius?
  cp <- which(di<iniR ,arr.ind=TRUE)
  ncp <- length(cp)
  # how many points should be used for interpolation? min 4 and max 10:
  n <- max(4,min(10,ncp))   # minimale und maximale Interpolationspunktzahl, evtl anders wählen
  #subscripts ij for increasing distance (die sortierten)
  ij <- sort(di, index.return=TRUE)$ix
  # final radius r:                     ########### evtl extra 0/1-Gewichtungsfaktor?, wenig hilfreich wegen fallunterscheidung
  if (ncp<=4) { radius<-di[ij][4+1]
  } else {
    if (ncp<=10) {radius<-iniR
    } else {
      radius<-di[ij][10+1]
    }
  }
  # Compute weights si (for N data points):
  si <- numeric(N)
  si <- si+(di>0)*(di<=radius/3)*(1/di)
  si <- si+(di>radius/3)*(di<=radius)*(27/(4*radius))*(((di/radius)-1)^2)
  return(si)  # Anmerkung: sum(si) nicht 1
}

#' Get cosinus of pairwise angles between data points and interpolation point as fixed point
#'
#' This function is basically a helper function for mvInterpolation::ti_get and not intended for direct use.
#'
#' @param data A matrix of data points
#' @param xinp Coordinates of the fixed point of the angles; usually the point where interpolation has to be computed
#' @param distance Distances between xinp and data can be given as an optional argument, if they have been coputed before, such that a new computation is not necessary
#' @return The cosini of the pairwise angles between two data points with fixed point xinp
#' @examples
#' coords <- matrix(runif(10*4)*10, ncol=4)
#' cosangmx(coords,rep(4,4))
#' @export
cosangmx <- function(data, xinp, distance=NULL) {
  d <- length(xinp) #4 #Anzahl der Koordinaten  aus Länge von xinp auslesen
  if (is(data, "SpatialPointsDataFrame")==TRUE)
    data <- as.matrix(data@data)
  N <- dim(data)[1] #Anzahl der Punkte

  #wenn distance nicht übergeben wurde (eigentlich nur dann...)
  if (is.null(distance)) distance <- as.matrix(dist(rbind(xinp,data[,1:d])))[-1,1]

  # cosinus vom Winkel von DiPDj: (ergibt sich durch inneres Produkt)
  xj14 <- -data[,1:d]
  rownames(xj14) <- paste0("-Xj",as.character(1:dim(xj14)[1]))
  pxj14 <- sweep(xj14, 2, xinp, "+") #P-X_j
  rownames(pxj14) <- paste0("P",rownames(pxj14))

  cosmat <- matrix(nrow=N, ncol=N)
  colnames(cosmat) <- paste0("i=", as.character(1:N))
  rownames(cosmat) <- paste0("j=", as.character(1:N))
  for (i in 1:N) {
    pxi <- (xinp-data[i,1:d])
    #rownames(pxi) <- "P-Xi"
    pxipxj14 <- sweep(pxj14,2, pxi, "*")
    #rownames(pxipxj14) <- paste0("(P-xi)(",rownames(pxipxj14),")")
    cosmat[i,] <- apply(pxipxj14,1,sum)/(distance[i]*distance[1:N])
    #names(cosangi) <- paste0("cosangi", as.character(1:dim(xj14)[1]))
  }
  return(cosmat)
}

#' Get weights for interpolation by angular distance weighting (ADW)
#'
#' Interpolation weights for angular distance weighting of the nearest neighbors (4 to 10 neighbors, depending on the distance)
#'
#' @param data A matrix of the points for which the function is known
#' @param xinp Coordinates, where interpolation has to be computed
#' @param rad Function which gives the initial radius for the observations considered in the interpolation
#' @return The weights for the ADW interpolation based on the nearest neighbors
#' @examples
#' mydata <- matrix(runif(10*4)*10, ncol=4)
#' mydata <- cbind(mydata,abs(apply(mydata,1,sum)-3),abs(apply(mydata,1,prod)-4))
#' ti_get(mydata, c(4,4,4,4))
#' @export
ti_get <- function(data, xinp, rad=getR) {
  d <- length(xinp) #Anzahl der Koordinaten  aus xinp auslesen
  if (is(data, "SpatialPointsDataFrame")==TRUE) data <- as.matrix(data@data)
  N <- dim(data)[1] #Anzahl der Punkte

  #iniR <- getR(d=d, N=N)
  if (is.numeric(rad)) {iniR <- rad
  } else {iniR <- rad(dat=data, d=d)}

  distance <- as.matrix(dist(rbind(xinp,data[,1:d])))[-1,1]
  #wieviele Di sind näher am Interpolationspunkt als der Initialradius?
  ncp <- length(which(distance<iniR ,arr.ind=TRUE))
  # wie viele Punkte benutze ich also? min 4 und max 10:
  n <- max(4,min(10,ncp))
  #subscripts ij for increasing distance (die sortierten)
  ij <- sort(distance, index.return=TRUE)$ix

  #cosmat_all <- cosangmx(data=data,xinp=xinp)
  cosmat <- cosangmx(data=data[ij[1:n],],xinp=xinp)  #arbeite hier NUR mit data im betrachteten Radius
  si<-si_get(data=data, xinp=xinp) # übergebe hier ALLE Daten
  si<-si[ij[1:n], drop=FALSE]
  ti <- apply( si*(1-cosmat), 2, sum )/sum(si)
  wi <- si^2*(1+ti)       ## Achtung, ^2 evtl anzupassen
  # Nullen wieder einfügen
  ti_all <- numeric(N)
  ti_all[ij[1:n]] <- wi[1:n]
  return(ti_all)
}

#' Get drift based on slopes
#'
#' The drift is needed for interpolation by angular distance weighting extended by slopes.
#' Interpolation weights are based on the nearest neighbors (4 to 10 neighbors, depending on the distance).
#'
#' @param data A matrix of the points for which the function is known
#' @param xinp Coordinates, where interpolation has to be computed
#' @param rad Function which gives the initial radius for the observations considered in the interpolation
#' @param vparam Scalar, parameter for the slopes
#' @return The weights for the ADW interpolation considering slopes and based on the nearest neighbors
#' @examples
#' mydata <- matrix(runif(10*4)*10, ncol=4)
#' mydata <- cbind(mydata,abs(apply(mydata,1,sum)-3),abs(apply(mydata,1,prod)-4))
#' slopezi_get(mydata, c(4,4,4,4))
#' @export
slopezi_get <- function(data, xinp, rad=getR, vparam=0.1) {
  d <- length(xinp) #Anzahl der Koordinaten  aus xinp auslesen   ## brauchen wir nur deswegen xinp???
  if (is(data, "SpatialPointsDataFrame")==TRUE) {data <- as.matrix(data@data)
  }else{data<-as.matrix(as.data.frame(data))} #strange but as it seems necessary, cf: https://stackoverflow.com/questions/15084803/length-of-dimnames-2-not-equal-to-array-extent-on-one-of-two-very-similar
  N <- dim(data)[1] #Anzahl der Punkte
  R<-(dim(data)[2]-d)

  distance_g <- as.matrix(dist(rbind(xinp,data[,1:d])))[-1,1]

  if (is.numeric(rad)) {iniR <- rad
  } else {iniR <- rad(dat=data, d=d)}    #Anfangsradius wird nur hier berechnet, nicht für die um 1 Element verkleinerten Subsets der Datengrundlage
  #wir müssen *jeden* Datenpunkt anschauen
  result <- array(dim=c(N,d,R))#alt:#matrix(nrow=N,ncol=d)
  for (i in 1:N) {
    #Anzahl der Nachbarn im Suchradius
    #Abstände von Xi
    distance <- as.matrix(dist(data[,1:d]))[i,]
    #für welche Indices (außer i) sind die Di näher am Interpolationspunkt als der Initialradius?
    cp <- setdiff(which(distance<iniR ,arr.ind=TRUE),i)
    ncp <- length(cp)
    # wie viele Punkte benutze ich also? min 4 und max 10:
    n <- max(4,min(10,ncp))
    # betrachtete Punkte
    #supscripts ij for increasing distance (die sortierten)
    ij <- sort(distance, index.return=TRUE)$ix
    #wirklich betrachtete ids   # Achtung: im Folgenden darf ich den Index i ja nicht verwenden...
    ciprime2_id <- setdiff(ij,i)[1:n]   #alternativ: einfach erstes Element weglassen und dann 1:n
    #wirklich betrachtete Daten
    ciprime2 <- data[ciprime2_id,]
    rownames(ciprime2)<-paste0("j=",ciprime2_id)
    PiDiffPjs <- sweep(ciprime2,2,data[i,],"-")
  ####  if (is.null(colnames(PiDiffPjs))) colnames(PiDiffPjs)[1:d] <- paste0("C",1:d)
    colnames(PiDiffPjs) <- paste0(colnames(PiDiffPjs),"_j-",colnames(PiDiffPjs),"_i")
    wi<-ti_get(data=data,xinp=xinp)[ciprime2_id]
    prodimzaehler <- array(dim=c(n,d,R))
    for (r in 1:R) {
      prodimzaehler[,,r] <- sweep(PiDiffPjs[,1:d],1,PiDiffPjs[,d+r],"*")
    }
    dimnames(prodimzaehler) = dimnames(PiDiffPjs[,1:d])
    dimnames(prodimzaehler)[[3]] <- paste0("response",1:R)
    zaehler <- apply(prodimzaehler*wi / (distance[ciprime2_id])^2,c(2,3),sum)
    #teste: bei zwei gleichen Responses kommt das gleiche raus?
    nenner <- sum(wi)
    result[i,,] <- zaehler/nenner
    #^^vielleicht ist das nciht so gedacht. wenn wi=0 dann ist auch slope ausgelöscht? nur für entsprechendes j,insgesamt nicht
  }
##  colnames(result)<- paste0("Slope",1:d)
##  rownames(result) <- paste0("i=",1:N)
  dimnames(result)[[3]] <- paste0("response",1:R)

  #define v
  zaehl <- vparam*apply(data[,-(1:d), drop=FALSE],2,function(x) diff(range(x)))
  nenn <-  apply(apply(result^2,c(1,3),sum),2,max)^(1/2)
  v <- (zaehl/nenn)

  # Ai *  x-xi | y-yi | ... | nicht mehr: z-zi
  rs <- sweep(result,c(1,2),sweep(-data[,1:d],2,xinp,"+"),"*")
  return(apply(rs,c(1,3),sum)* sweep(1/sweep(matrix(distance_g,nrow=N,ncol=R),2,v,"+"),2,v,"*")  )
}

#' Do the multivaiate interpolation based on weights and drift functions
#'
#' The drift is needed for interpolation by angular distance weighting extended by slopes.
#' Interpolation weights include either all observations or are based on the nearest neighbors
#' (4 to 10 neighbors, depending on the distance).
#'
#' The interpoliation works for a multidimensional input space and also for multi-responses.
#'
#' @param data A matrix of the points for which the function is known
#' @param xinp Coordinates, where interpolation has to be computed
#' @param weights Vector of the weights for the interpolation
#' @param drift Matrix of the drift for the interpolation
#' @param p Scalar, parameter for the weights, currently not used
#' @param vparam Scalar, parameter for the slopes, currently not used
#' @return The interpolated value(s) at xinp
#' @examples
#' mydata <- matrix(runif(10*4)*10, ncol=4)
#' mydata <- cbind(mydata,abs(apply(mydata,1,sum)-3),abs(apply(mydata,1,prod)-4))
#' inp <- rep(4,4)
#' #Inverse distance weighting (based on all observations)
#' wintpl(data=mydata, xinp=inp, weights=di_get(data=mydata, xinp=inp,2))
#' #Inverse distance weighting (based on nearest neighbors)
#' wintpl(data=mydata, xinp=inp, weights=si_get(data=mydata, xinp=inp)^2)
#' #Angular distance weighting
#' wintpl(data=mydata, xinp=inp, weights=ti_get(data=mydata, xinp=inp))
#' #Angular distance weighting under consideration of slopes
#' wintpl(data=mydata, xinp=inp, weights=ti_get(data=mydata, xinp=inp), drift=slopezi_get(mydata,inp))
#' @export
wintpl <- function(data, xinp, weights, drift=0, p=NULL, vparam=NULL) {
  d <- length(xinp) #read number of coordinates from xinp
  if (is(data, "SpatialPointsDataFrame")==TRUE)
    data <- as.matrix(data@data)
  N <- dim(data)[1] #number of points

  #If interpolation point is known don't interpolate
  r0 <- data[apply(data[,1:d], 1, function(t) identical(as.numeric(t), as.numeric(xinp))),-(1:d),drop=FALSE] #auch mit mehreren response-dimensionen #evtl maschinengenauigkeit hier aussteuern
  if (length(r0)>0) return(r0[1])
  #else do interpolation
  wi <- weights # si_get(data=smpl, xinp=inp)^2   #Achtung, Exponent hier oder in si_get evtl höher zu wählen!
  return(apply((data[,-(1:d),drop=FALSE]+drift)*wi,2,sum)/sum(wi))
  #seltsam: Rückgabetyp anders wenn Interpolationspunkt in Datensatz ist
}

