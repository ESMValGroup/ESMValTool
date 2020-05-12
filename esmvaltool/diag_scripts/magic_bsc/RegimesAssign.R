anom2regime <- function(ref, target, method = "distance", lat) { # nolint
  posdim <- which(names(dim(ref)) == "nclust")
  poslat <- which(names(dim(ref)) == "lat")
  poslon <- which(names(dim(ref)) == "lon")

  nclust <- dim(ref)[posdim]

  if (all(dim(ref)[-posdim] != dim(target))) {
    stop("The target should have the same dimensions [lat,lon] that
         the reference ")
  }
  if (is.null(names(dim(ref))) | is.null(names(dim(target)))) {
    stop(
      "The arrays should include dimensions names ref[nclust,lat,lon]
      and target [lat,lon]"
    )
  }

  if (length(lat) != dim(ref)[poslat]) {
    stop("latitudes do not match with the maps")
  }

  # This dimensions are reorganized
  ref <- aperm(ref, c(posdim, poslat, poslon))
  target <-
    aperm(target, c(which(names(dim(
      target
    )) == "lat"), which(names(dim(
      target
    )) == "lon")))

  # weights are defined
  latWeights <-
    InsertDim(sqrt(cos(lat * pi / 180)), 2, dim(ref)[3]) # nolint


  rmsdiff <- function(x, y) {
    dims <- dim(x)
    ndims <- length(dims)
    if (ndims != 2 | ndims != length(dim(y))) {
      stop("x and y should be maps")
    }
    map_diff <- NA * x
    for (i in 1:dims[1]) {
      for (j in 1:dims[2]) {
        map_diff[i, j] <- (x[i, j] - y[i, j])^2
      }
    }
    rmsdiff <- sqrt(mean(map_diff, na.rm = TRUE))
    return(rmsdiff)
  }

  if (method == "ACC") {
    corr <- rep(NA, nclust)
    for (i in 1:nclust) {
      corr[i] <-
        ACC(
          InsertDim(InsertDim(
            # nolint
            InsertDim(ref[i, , ] * latWeights, 1, 1), 2, 1 # nolint
          ), 3, 1),
          InsertDim(InsertDim(
            # nolint
            InsertDim(target * latWeights, 1, 1), 2, 1 # nolint
          ), 3, 1)
        )$ACC[2]
    }
    assign <- which(corr == max(corr))
  }

  if (method == "distance") {
    rms <- rep(NA, nclust)
    for (i in 1:nclust) {
      rms[i] <-
        rmsdiff(ref[i, , ] * latWeights, target * latWeights) # nolint
    }
    assign <- which(rms == min(rms, na.rm = TRUE))
  }
  return(assign)
}

RegimesAssign <- function(var_ano, ref_maps, lats, # nolint
                          method = "distance") {
  posdim <- which(names(dim(ref_maps)) == "nclust")
  poslat <- which(names(dim(ref_maps)) == "lat")
  poslon <- which(names(dim(ref_maps)) == "lon")
  poslat_ano <- which(names(dim(var_ano)) == "lat")
  poslon_ano <- which(names(dim(var_ano)) == "lon")

  nclust <- dim(ref_maps)[posdim]
  nlat <- dim(ref_maps)[poslat]
  nlon <- dim(ref_maps)[poslon]


  if (is.null(names(dim(ref_maps))) |
    is.null(names(dim(var_ano)))) {
    stop(
      "The arrays should include dimensions names ref[nclust,lat,lon]
      and target [lat,lon]"
    )
  }

  if (length(lats) != dim(ref_maps)[poslat]) {
    stop("latitudes do not match with the maps")
  }
  print(str(var_ano))
  assign <-
    Apply(
      data = list(target = var_ano),
      margins = c((seq_along(dim(
        var_ano
      )))[-c(poslat_ano, poslon_ano)]),
      fun = "anom2regime",
      ref = ref_maps,
      lat = lats,
      method = method
    )

  if (poslat_ano < poslon_ano) {
    dim_order <- c(nlat, nlon)
  } else {
    dim_order <- c(nlon, nlat)
  }

  anom_array <-
    array(var_ano, dim = c(
      prod(dim(var_ano)[-c(poslat_ano, poslon_ano)]),
      dim_order
    ))

  rm(var_ano)

  index <- as.vector(assign$output1)
  recon <-
    Composite(var = aperm(anom_array, c(3, 2, 1)), occ = index)
  freqs <- rep(NA, nclust)
  for (n in 1:nclust) {
    freqs[n] <- (length(which(index == n)) / length(index)) * 100
  }
  output <-
    list(
      composite = recon$composite,
      pvalue = recon$pvalue,
      cluster = assign$output1,
      frequency = freqs
    )
  return(output)
}
