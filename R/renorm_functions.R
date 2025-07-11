## Compute the MAD
MAD.scorer <- function(x = NULL) {
  return(median(abs(diff(as.numeric(x[!is.na(x)])))))
}
## Compute the MALMM
MALMM.scorer <- function(x = NULL, smo = 399) {
  malmm <- median(abs(x - median(x, na.rm = TRUE)), na.rm = TRUE)
  return(malmm)
}
## Compute the SSAD on a running median
SSAD.scorer <- function(x = NULL, smo = 399) {
  rm.diff <- diff(as.numeric(runmed(x[!is.na(x)], smo)))
  ssad <- sum(abs(rm.diff[rm.diff != 0]))
  return(ssad)
}
## Compute the MARMMM on a running median
MARMMM.scorer <- function(x = NULL, smo = 399) {
  rmed <- as.numeric(runmed(x[!is.na(x)], smo))
  marmmm <- median(abs(rmed - median(x, na.rm = TRUE)), na.rm = TRUE)
  return(marmmm)
}

## L2R tf fit loop function
l2r.fitloop_bak <- function(l2rObj, tfd, smo = 399, method = "loess") {
  
  ### FITLOOP
  minitf <- tfd[,-c(1:4), drop = FALSE]
  tfheads <- colnames(minitf)
  b <- ncol(minitf)+1
  posfit <- c()
  
  tmsg(paste0("Init (", l2rObj$rm.mad, ")"))
  
  while ( (b != 1) & (ncol(minitf) != 0) ) {
    
    biggy <- list()
    biggy <- append(biggy, list(l2rObj))
    rmtest <- l2rObj$rm.mad
    for (z in 1:length(minitf)) {
      if (method == "loess") tempfit <- l2r.fit(biggy[[1]]$l2r, minitf[,z], smo)
      if (method == "pcs") tempfit <- l2r.pcs(biggy[[1]]$l2r, minitf[,z], smo)
      if (method == "spline") tempfit <- l2r.spline(biggy[[1]]$l2r, minitf[,z], smo)
      biggy <- append(biggy, list(tempfit))
      rmtest <- c(rmtest, tempfit$rm.mad)
      # message(paste0(z, " / ", tempfit$rm.mad))
    }
    b <- which.min(rmtest)
    if (b > 1) {
      tmsg(paste0(" Positive fit with ", tfheads[b-1], " (", min(rmtest), ")"))
      posfit <- c(posfit, tfheads[b-1])
      l2rObj <- biggy[[b]]
      minitf <- as.data.frame(minitf[,-c(b-1)])
      tfheads <- tfheads[-c(b-1)]
    }
  }
  return(list(l2r = l2rObj, pos = posfit))
}

# normloop.res <- list(
#     data = input.data, 
#     renorm = l2r.fitloop(
#       l2r = input.data$L2R,
#       chr = input.data$chr,
#       tfd = RNdata, 
#       smo = smo,
#       sex.chr = sex.chr)
#   )
  
## L2R tf fit loop function
l2r.fitloop <- function(l2r = NULL, chr = NULL, tfd = NULL, method = 'loess', sex.chr = NULL, smo = 399) {
  
  ### Init
  # tfheads <- colnames(tfd)
  b <- ncol(tfd)+1
  posfit <- c()
  
  ## Get sex chr pos
  if (!is.null(sex.chr)) {
    sex.idx <- chr %in% sex.chr
    if (any(sex.idx)) {
      auto.med <- median(l2r[!sex.idx], na.rm = TRUE)
      ks.diff <- sapply(sex.chr, function(k) {
        k.idx <- chr %in% k
        k.med <- median(l2r[k.idx], na.rm = TRUE)
        k.diff <- if (any(k.idx)) auto.med - k.med else 0
        return(k.diff)
      })
    }
  }
  
  ## Compute initial MARMMM
  top.marmmm <- MARMMM.scorer(x = l2r, smo = smo)
  tmsg(paste0("Init (", top.marmmm, ")"))
  
  while ( (b != 1) & (ncol(tfd) != 0) ) {
    
    biggy <- list(ori = l2r)
    marmmmz <- c(top.marmmm)
    for (z in colnames(tfd)) {
      l2r.tmp <- l2r
      nona <- !is.na(l2r)
      ## Regress track
      if (method == "loess") tempfit <- l2r.loess(l2r = l2r[nona], tf = tfd[[z]][nona])
      if (method == "spline") tempfit <- l2r.spline(l2r = l2r[nona], tf = tfd[[z]][nona])
      if (method == "pcs") tempfit <- l2r.pcs(l2r = l2r[nona], tf = tfd[[z]][nona])
      l2r.tmp[nona] <- tempfit
      tempfit <- l2r.tmp
      ## Correct sex chromosomes
      if(any(sex.idx)) {
        auto.med <- median(tempfit[!sex.idx], na.rm = TRUE)
        for (k in sex.chr) {
          k.idx <- chr == k
          if (any(k.idx)) {
            k.med <- median(tempfit[k.idx], na.rm = TRUE)
            tempfit[k.idx] <- tempfit[k.idx] + (auto.med - k.med) - ks.diff[[k]]
          }
        }
      }
      ## Center
      tempfit <- tempfit - median(tempfit, na.rm = TRUE)
      ## Compute new MARMMM
      z.marmmm <- MARMMM.scorer(x = tempfit, smo = smo)
      ## Save
      biggy[[z]] <- tempfit
      marmmmz <- c(marmmmz, z.marmmm)
      rm(tempfit, z.marmmm)
    }
    b <- which.min(marmmmz)
    if (b > 1) {
      tmsg(paste0(" Positive fit with ", names(biggy)[b], " (", marmmmz[b], ")"))
      
      posfit <- c(posfit, names(biggy)[b])
      l2r <- biggy[[b]]
      top.marmmm <- marmmmz[b]
      tfd <- tfd[,-c(b-1), drop = FALSE]
    }
  }
  return(list(l2r = l2r, pos = posfit))
}

## FONCTION LOESSFIT
l2r.loess <- function(l2r = NULL, tf = NULL) {
  ## Check NAs
  if (any(is.na(tf))) tf <- stats::approxfun(seq_along(tf), tf, rule = 2)(seq_along(tf))
  ## Fit regression
  l2fN <- limma::loessFit(y = l2r, x = tf)
  ## Apply regression
  l2N <- l2r-l2fN$fitted
  ## Center
  l2N <- l2N - median(l2N, na.rm = TRUE)
  return(l2N)
}

## FONCTION PERCENTILE SCALING
l2r.pcs <- function(l2r = NULL, tf = NULL) {
  ## Check NAs
  if (any(is.na(tf))) tf <- stats::approxfun(seq_along(tf), tf, rule = 2)(seq_along(tf))
  ## Apply regression
  l2N <- GCnorm.pcs(measures = l2r, gc = tf)
  ## Center
  l2N <- l2N - median(l2N, na.rm = TRUE)
  return(l2N)
}

## FONCTION SPLINE (cubic)
l2r.spline <- function(l2r = NULL, tf = NULL, df = 3) {
  ## Check NAs
  if (any(is.na(tf))) tf <- stats::approxfun(seq_along(tf), tf, rule = 2)(seq_along(tf))
  ## Fit regression
  l2fN <- stats::lm(l2r ~ splines::bs(tf, degree = df))
  ## Apply regression
  l2N <- stats::residuals(l2fN)
  ## Center
  l2N <- l2N - median(l2N, na.rm = TRUE)
  return(l2N)
}

GCnorm.pcs <- function(measures = NULL, gc = NULL) {
  # message("Performing GC PC-scaling ...")
  gcpc.idx <- sapply(sort(unique(gc)), function(x) { return(which(gc == x)) })
  names(gcpc.idx) <- sort(unique(gc))
  gcpc.value <- sapply(1:length(gcpc.idx), function(x) { return(measures[gcpc.idx[[x]]]) })
  gcpc.med <- sapply(gcpc.value, median)
  measures.new <- rep(NA, length(gc))
  for (x in 1:length(gcpc.med))  measures.new[gcpc.idx[[x]]] <- measures[gcpc.idx[[x]]] - gcpc.med[x]
  return(measures.new)
}

## Main renormalization function
renorm.go_bak <- function(input.data = NULL, renorm.data = NULL, track.type = "gc", sex.chr = NULL, smo = 399, arraytype = NULL, genome = NULL) {
  if (!is.null(renorm.data)) {
    # load(renorm.rda, envir = environment())
    rn.arraytype <- renorm.data$info$value[renorm.data$info$key == "array_type"]
    # rn.genome <- renorm.data$info$value[renorm.data$info$key == "genome-version"]
    # rn.track.type <- renorm.data$info$value[renorm.data$info$key == "track_type"]
    # if ((rn.track.type != track.type) | (rn.arraytype != arraytype) | (rn.genome != genome)) stop(tmsg(paste0("Provided renormalization pack is not as intended ! Expected [", track.type, ", ", arraytype, ", ", genome, "], got [", rn.track.type, ", ", rn.arraytype, ", ", rn.genome, "] !")), call. = FALSE)
  } else {
    RN.pkg.name <- "affy.CN.norm.data"
    if (!(RN.pkg.name %in% installed.packages())) stop(tmsg(paste0("Package ", RN.pkg.name, " not found !")), call. = FALSE)
    RN.file <- system.file(paste0("data/", arraytype, ".", genome, ".tracks.rda"), package = RN.pkg.name)
    if (RN.file == "") stop(tmsg(paste0("Could not find a normalization tracks data package for [", arraytype, ", ", genome, "] in package '", RN.pkg.name, "' ! Please build your own GC / Wave data pack with ", RN.pkg.name, "::affy.tracks.compute(), and submit it using the 'renorm.data' option.")), call. = FALSE)
    data(list = paste0(arraytype, ".", genome, ".tracks"), package = RN.pkg.name, envir = environment())
  }
  # print(str(RNdata))
  # print(str(rownames(input.data)))
  RNdata <- renorm.data$tracks[renorm.data$tracks$ProbeSetName %in% input.data$ProbeSetName,]
  input.data <- input.data[input.data$ProbeSetName %in% RNdata$ProbeSetName,]
  # print(str(input.data))
  if (!all(unique(RNdata$ProbeSetName == input.data$ProbeSetName))) stop(tmsg(paste0(track.type, " data and L2R data are not synched, or ordered differently !")), call. = FALSE)
  # ndata <- data.frame(chr = paste0("chr", input.data$chrs), start = input.data$pos, end = input.data$pos, name = rownames(input.data), RNdata[,-c(1:4), drop = FALSE], stringsAsFactors = FALSE)
  ndata <- data.frame(chr = input.data$chr, start = input.data$pos, end = input.data$pos, name = input.data$ProbeSetName, RNdata[,-c(1:4), drop = FALSE], stringsAsFactors = FALSE)
  rm(RNdata, renorm.data)
  # print(str(ndata))
  rm.diff <- diff(as.numeric(runmed(input.data$L2R[!is.na(input.data$L2R)], smo)))
  my.rm.mad <- sum(abs(rm.diff[rm.diff != 0]))
  # print(paste0("RMMAD ", my.rm.mad))
  # print(paste0(summary(my.rm.mad)))
  normloop.res <- list(data = input.data, renorm = l2r.fitloop(l2rObj = list(l2r=input.data$L2R, rm.mad = my.rm.mad), tfd = ndata, smo = smo))
  
  return(normloop.res)
  
  # input.data[[paste0("L2R.", pack.type)]] <- normloop.res$l2r$l2r + median(input.data$L2R, na.rm = TRUE)
  # input.data$L2R <- input.data[[paste0("L2R.", pack.type)]]
  # return(input.data)
  
}

## Main renormalization function
renorm.go <- function(input.data = NULL, renorm.data = NULL, track.type = 'gc', method = 'loess', sex.chr = sex.chr, smo = 399, arraytype = NULL, genome = NULL) {
  if (is.null(renorm.data)) {
    RN.pkg.name <- "affy.CN.norm.data"
    if (!(RN.pkg.name %in% installed.packages())) stop(tmsg(paste0("Package ", RN.pkg.name, " not found (required) !")), call. = FALSE)
    RN.file <- system.file(paste0("data/", arraytype, ".", genome, ".tracks.rda"), package = RN.pkg.name)
    if (RN.file == "") stop(tmsg(paste0("Could not find a normalization tracks data package for [", arraytype, ", ", genome, "] in package '", RN.pkg.name, "' ! Please build your own GC / Wave data pack with ", RN.pkg.name, "::affy.tracks.compute(), and submit it using the 'renorm.data' option.")), call. = FALSE)
    data(list = paste0(arraytype, ".", genome, ".tracks"), package = RN.pkg.name, envir = environment())
  }
  ## Synch checks
  RinI <- renorm.data$tracks$pos$ProbeSetName %in% input.data$ProbeSetName
  IinR <- input.data$ProbeSetName %in% renorm.data$tracks$pos$ProbeSetName
  if (!all(unique(renorm.data$tracks$pos$ProbeSetName[RinI] == input.data$ProbeSetName[IinR]))) stop(tmsg(paste0(track.type, " data and L2R data are not synched, or ordered differently !")), call. = FALSE)
  
  RNdata <- if (track.type == 'pos') data.frame(width = (renorm.data$tracks$pos$end - renorm.data$tracks$pos$start +1))[RinI,, drop = FALSE] else renorm.data$tracks[[track.type]][RinI,, drop = FALSE]
  input.data <- input.data[IinR,]
  
  # ndata <- data.frame(chr = input.data$chr, start = input.data$pos, end = input.data$pos, name = input.data$ProbeSetName, RNdata, stringsAsFactors = FALSE)
  # rm(RNdata, renorm.data)
  
  # rm.diff <- diff(as.numeric(runmed(input.data$L2R[!is.na(input.data$L2R)], smo)))
  # my.rm.mad <- sum(abs(rm.diff[rm.diff != 0]))
  # print(paste0("RMMAD ", my.rm.mad))
  # print(paste0(summary(my.rm.mad)))
  normloop.res <- list(
    data = input.data, 
    renorm = l2r.fitloop(
      l2r = input.data$L2R,
      chr = input.data$chr,
      tfd = RNdata,
      method = method,
      smo = smo,
      sex.chr = sex.chr)
  )
  
  return(normloop.res)
  
  # input.data[[paste0("L2R.", pack.type)]] <- normloop.res$l2r$l2r + median(input.data$L2R, na.rm = TRUE)
  # input.data$L2R <- input.data[[paste0("L2R.", pack.type)]]
  # return(input.data)
  
}

