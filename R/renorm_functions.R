## L2R tf fit loop function
l2r.fitloop <- function(l2rObj, tfd, smo = 399, method = "loess") {
  
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

## FONCTION LOESSFIT
l2r.fit <- function(l2r, tf, smo) {
  l2fN <- limma::loessFit(l2r, tf)
  l2N <- l2r-l2fN$fitted
  rm.diff <- diff(as.numeric(runmed(l2N[!is.na(l2N)], smo)))
  Nrmspread <- sum(abs(rm.diff[rm.diff != 0]))
  return(list(l2r=l2N, rm.mad=Nrmspread))
}

## FONCTION PERCENTILE SCALING
l2r.pcs <- function(l2r, tf, smo) {
  l2N <- GCnorm.pcs(l2r, tf)
  rm.diff <- diff(as.numeric(runmed(l2N[!is.na(l2N)], smo)))
  Nrmspread <- sum(abs(rm.diff[rm.diff != 0]))
  return(list(l2r=l2N, rm.mad=Nrmspread))
}

GCnorm.pcs <- function(measures=NULL, gc=NULL) {
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
renorm.go <- function(input.data = NULL, renorm.rda = NULL, track.type = "GC", smo = 399, arraytype = NULL, genome = NULL) {
  if (!is.null(renorm.rda)) {
    load(renorm.rda, envir = environment())
    rn.arraytype <- renorm.data$info$value[renorm.data$info$key == "array_type"]
    rn.genome <- renorm.data$info$value[renorm.data$info$key == "genome-version"]
    rn.track.type <- renorm.data$info$value[renorm.data$info$key == "track_type"]
    if ((rn.track.type != track.type) | (rn.arraytype != arraytype) | (rn.genome != genome)) stop(tmsg(paste0("Provided renormalization pack is not as intended ! Expected [", track.type, ", ", arraytype, ", ", genome, "], got [", rn.track.type, ", ", rn.arraytype, ", ", rn.genome, "] !")))
  } else {
    RN.pkg.name <- "affy.CN.norm.data"
    if (!(RN.pkg.name %in% installed.packages())) stop(tmsg(paste0("Package ", RN.pkg.name, " not found !")))
    RN.file <- system.file(paste0("data/", arraytype, ".", genome, ".", track.type, ".rda"), package = RN.pkg.name)
    if (RN.file == "") stop(tmsg(paste0("Could not find a ", track.type, " data package for [", arraytype, ", ", genome, "] in package '", RN.pkg.name, "' ! Please build your own Wave data pack with ", RN.pkg.name, "::affy.wave.compute() and submit it using the 'renorm.rda' option.")))
    data(list = paste0(arraytype, ".", genome, ".", track.type), package = RN.pkg.name, envir = environment())
  }
  # print(str(RNdata))
  # print(str(rownames(input.data)))
  RNdata <- renorm.data$tracks[renorm.data$tracks$ProbeSetName %in% input.data$ProbeSetName,]
  input.data <- input.data[input.data$ProbeSetName %in% RNdata$ProbeSetName,]
  # print(str(input.data))
  if (!all(unique(RNdata$ProbeSetName == input.data$ProbeSetName))) stop(tmsg(paste0(track.type, " data and L2R data are not synched, or ordered differently !")))
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

