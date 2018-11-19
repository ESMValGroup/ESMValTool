AtomicWeatherRegime <- function(
  data, EOFS = TRUE, neofs = 30, threshold = NULL, lon = NULL, lat = NULL,
  ncenters = NULL, method = "kmeans", nstart = 30) {
  names(dim(data)) <- c("sdate", "ftime", "lat", "lon")
  sdate <- which(names(dim(data)) == "sdate")
  ftime <- which(names(dim(data)) == "ftime")
  nftimes <- dim(data)[ftime]
  nsdates <- dim(data)[sdate]
  lon2 <- which(names(dim(data)) == "lon")
  lat2 <- which(names(dim(data)) == "lat")
  data <- aperm(data, c(ftime, sdate, lat2, lon2))
  nlon <- dim(data)[lon2]
  nlat <- dim(data)[lat2]
  dim(data) <- c(nftimes *  nsdates, nlat, nlon)

  if (is.null(ncenters)) {
    stop("ncenters must be specified")
  }
  if (EOFS  == TRUE && is.null(lon)) {
    stop("longitudes must be specified")
  }
  if (EOFS  == TRUE && is.null(lat)) {
    stop("latitudes must be specified")
  }

  if (EOFS  == TRUE) {
    #data <- princomp(data[,subs], cor = FALSE)$scores
    dataPC <- EOF(
      data,
      lat = as.vector(lat),
      lon = as.vector(lon),
      neofs = neofs
    )
    if (is.null(threshold)){
      threshold <- sum(dataPC$var)
      cluster_input <- dataPC$PC
    } else {
      threshold <- threshold
      minPC <- head(as.numeric(which(cumsum(dataPC$var) > threshold)), 1)
      cluster_input <- dataPC$PC[, 1 : minPC]    }
  } else {
    cluster_input <- data
    latWeights <- InsertDim(
      InsertDim(cos(lat*pi/180), 1, nftimes*nsdates),
      3,
      nlon
    )
    cluster_input <- cluster_input * latWeights
    dim(cluster_input) <- c(nftimes * nsdates, nlat * nlon)
  }
  if (method == "kmeans") {
    result <- kmeans(cluster_input, centers = ncenters,
                     iter.max = 100, nstart = nstart, trace = FALSE)
    reconstructed <- array(0, c(ncenters, nlat, nlon))
    data <- aperm(data, c(2, 3, 1))
    reconstructed <- Composite(data, result$cluster)
    names(dim(reconstructed$composite)) <- c("lon", "lat", "cluster")
    cluster_timeseries <- list(lengths = c(), values = c())
    frequency <- persistence <- matrix(NA, nsdates, ncenters)
    for (i in 1 : nsdates) {
      occurences <- rle(
        result$cluster[((i * nftimes) + 1 - nftimes) : (i * nftimes)]
      )
      cluster_timeseries <- list(
        lengths = c(cluster_timeseries$lengths, occurences$lengths),
        values = c(cluster_timeseries$values, occurences$values)
      )
      for (j in 1 : ncenters) {
        total <- sum(occurences$lengths[occurences$values == j])
        frequency[i,j] <- (total /  nftimes) * 100
        persistence[i,j] <- mean(occurences$lengths[occurences$values == j])
      }
    }
    # for (i in 1 : ncenters) {
    #   for (j in 1 : nsdates) {
    #     timeseries_subset <- cluster_timeseries$lengths#[(1 +((nftimes -1) * j)) : (nftimes * j)]
    #     cluster_subset <- cluster_timeseries$values[(1 +((nftimes -1) * j)) : (nftimes * j)]
    #     print(timeseries_subset)
    #     print(cluster_subset)
    #     print((1 +((nftimes -1) * j)) : (nftimes * j))
    #     stop()
    #     frequency[i, j] <- (sum(timeseries_subset[cluster_subset == i]) / length(nftimes)) * 100
    #     persistence[i,j] <- mean(timeseries_subset[cluster_subset == i])
    #   }
    #   #frequency[i] <- (sum(cluster_timeseries$lengths[cluster_timeseries$values == i]) / length(result$cluster)) * 100
    # }
  } else {
    result <- hclust(dist(cluster_input), method = method)
    clusterCut <- cutree(result, ncenters)
    data <- aperm(data, c(3, 2, 1))
    result <- Composite(data, clusterCut)
  }
  if (method == "kmeans") {
    return(list(
      composite = reconstructed$composite,
      pvalue = reconstructed$pvalue,
      cluster = as.array(result$cluster),
      center = as.array(result$center),
      cluster_lengths = as.array(cluster_timeseries$lengths),
      cluster_values = as.array(cluster_timeseries$values),
      persistence = as.array(persistence), frequency = frequency))
  } else {
    return(list(
      composite = result$composite,
      pvalue = result$pvalue,
      cluster = as.array(clusterCut)
    ))
  }
}

WeatherRegime <- function(
  data, EOFS = TRUE, neofs = 30, threshold = NULL, lon = NULL, lat = NULL,
  ncenters = NULL, method = "kmeans", nstart = 30, iter.max = 100,
  ncores = NULL) {
  if (length(dim(data)) > 4) {
    sdate <- which(names(dim(data)) == "sdate")
    ftime <- which(names(dim(data)) == "ftime")
    lon_dim <- which(names(dim(data)) == "lon")
    lat_dim <- which(names(dim(data)) == "lat")
    dims <- c(1 : length(dim(data)))[-c(sdate, ftime, lon_dim, lat_dim)]
    data <- aperm(data, c(sdate, ftime, lat_dim, lon_dim, dims))
    margins <- 5 : length(dim(data))
    result <- Apply(
      data = list(data),
      margins = list(margins),
      AtomicFun = "AtomicWeatherRegime",
      EOFS = EOFS,
      neofs = neofs,
      threshold = threshold,
      lon = lon,
      lat = lat,
      ncenters = ncenters,
      method = method,
      ncores = ncores
    )
  } else {
    result <- AtomicWeatherRegime(
      data,
      EOFS = EOFS,
      neofs = neofs,
      threshold = threshold,
      lon = lon,
      lat = lat,
      ncenters = ncenters,
      method = method
    )
  }
}
